#pragma once
#include <algorithm>
#include <functional>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

// Naive, explicit de Bruijn graph oracle
class NaiveDBG {
public:
  using Node = std::string;
  using Label = char;
  using Edge = std::tuple<Node, Label, Node>;

  static NaiveDBG FromKmers(size_t k, const std::vector<std::string> &kmers) {
    NaiveDBG g;
    g.k_ = k;
    for (const auto &km : kmers) {
      if (km.size() != k)
        continue;
      Node u = km.substr(0, k - 1);
      Node v = km.substr(1, k - 1);
      Label c = km[k - 1];
      g.edges_.insert({u, c, v});
      g.nodes_.insert(u);
      g.nodes_.insert(v);
      g.kmers_.insert(km);
    }
    return g;
  }

  static NaiveDBG Union(const NaiveDBG &a, const NaiveDBG &b) {
    NaiveDBG g;
    g.k_ = a.k_;
    g.edges_ = a.edges_;
    g.edges_.insert(b.edges_.begin(), b.edges_.end());
    g.nodes_ = a.nodes_;
    g.nodes_.insert(b.nodes_.begin(), b.nodes_.end());
    g.kmers_ = a.kmers_;
    g.kmers_.insert(b.kmers_.begin(), b.kmers_.end());
    return g;
  }

  size_t K() const { return k_; }

  const std::set<Node> &Nodes() const { return nodes_; }
  const std::set<std::string> &Kmers() const { return kmers_; }
  const std::multiset<Edge> &Edges() const { return edges_; }

  bool Equals(const NaiveDBG &other) const {
    return CanonicalDigest() == other.CanonicalDigest();
  }

  // Canonical digest (stable, deterministic)
  std::string CanonicalDigest() const {
    std::vector<std::string> edge_strings;
    edge_strings.reserve(edges_.size());

    for (const auto &e : edges_) {
      edge_strings.push_back(std::get<0>(e) + "|" +
                             std::string(1, std::get<1>(e)) + "|" +
                             std::get<2>(e));
    }

    std::sort(edge_strings.begin(), edge_strings.end());

    std::ostringstream canonical;
    canonical << "k=" << k_ << "\n";
    for (const auto &s : edge_strings) {
      canonical << s << "\n";
    }

    // Hash the canonical serialization
    return StableHash(canonical.str());
  }

private:
  // Simple stable hash; replace with SHA-256 / BLAKE3 if desired
  static std::string StableHash(const std::string &s) {
    std::hash<std::string> hasher;
    size_t h = hasher(s);
    return std::to_string(h);
  }

  size_t k_{0};
  std::set<Node> nodes_;
  std::set<std::string> kmers_;
  std::multiset<Edge> edges_;
};

// Read-only semantic adaptor for any DBG implementation
class DBGView {
public:
  using Node = std::string;
  using Label = char;
  using Edge = std::tuple<Node, Label, Node>;

  virtual ~DBGView() = default;

  virtual size_t K() const = 0;

  // Semantic enumeration only
  virtual std::vector<Node> EnumerateNodes() const = 0;
  virtual std::vector<Edge> EnumerateEdges() const = 0;

  // Optional but helpful
  virtual std::vector<std::string> EnumerateKmers() const {
    std::vector<std::string> km;
    for (auto &e : EnumerateEdges()) {
      km.push_back(std::get<0>(e) + std::string(1, std::get<1>(e)));
    }
    return km;
  }
};

// Canonical digest for any DBGView
inline std::string CanonicalDigest(const DBGView &view) {
  std::vector<std::string> edges;
  edges.reserve(view.EnumerateEdges().size());

  for (const auto &e : view.EnumerateEdges()) {
    edges.push_back(std::get<0>(e) + "|" + std::string(1, std::get<1>(e)) +
                    "|" + std::get<2>(e));
  }

  std::sort(edges.begin(), edges.end());

  std::ostringstream canonical;
  canonical << "k=" << view.K() << "\n";
  for (const auto &s : edges) {
    canonical << s << "\n";
  }

  std::hash<std::string> hasher;
  return std::to_string(hasher(canonical.str()));
}

inline NaiveDBG ToNaiveDBG(const DBGView &view) {
  return NaiveDBG::FromKmers(view.K(), view.EnumerateKmers());
}