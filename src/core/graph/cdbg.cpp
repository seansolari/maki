
#include "maki/core/graph/cdbg.hpp"

ColouredGraphFiles graphFiles(fs::path base) {
  return ColouredGraphFiles{.l = base / "succ.dat",
                            .lR = base / "succ-rank.dat",
                            .lS = base / "succ-select.dat",
                            .W = base / "edges.dat",
                            .archive = base / "archive.dat",
                            .meta = base / "graph.dat"};
}

uint8_t ColouredGraph::block(size_t i) const {
  if (F[1] > i)
    return 0u | 0b1000u;
  else if (F[2] > i)
    return 1u | 0b1000u;
  else if (F[3] > i)
    return 2u | 0b1000u;
  else if (F[4] > i)
    return 3u | 0b1000u;
  else
    return 4u | 0b1000u;
}

uint8_t ColouredGraph::edge(size_t i) const { return W[i] & 0b0111; }

Dna4Sequence ColouredGraph::kmer(size_t i) const {
  size_t remaining = k + 1;
  Dna4Sequence s;
  s.reserve(remaining);
  std::optional i_ = i;
  // check first edge, it may be $
  if (auto e_ = edge(*i_); e_ > 0)
    s.push_back(parsing::dna5ToDna4(e_));
  i_ = bwd(*i_);
  --remaining;
  // continue
  while ((i_) && (remaining)) {
    s.push_back(parsing::dna5ToDna4(edge(*i_)));
    i_ = bwd(*i_);
    --remaining;
  }
  auto rs = s | std::views::reverse;
  return Dna4Sequence(rs);
}

// WARNING: no bounds check on input
std::optional<size_t> ColouredGraph::fwd(size_t i) const {
  uint8_t e = edge(i);
  if (e == 0u)
    return std::nullopt;
  size_t rnk = W.rank(i + 1, e | 0b1000u);
  return std::make_optional(lS(lR(F[e]) + rnk));
}

// WARNING: no bounds check on input
std::optional<size_t> ColouredGraph::bwd(size_t i) const {
  uint8_t c = block(i);
  if ((c ^ 0b1000u) == 0u)
    return std::nullopt;
  size_t r1 = (lR(i) + 1u), r2 = lR(F[c ^ 0b1000u]);
  return W.select(r1 - r2, c);
}

std::optional<size_t> ColouredGraph::outgoing(size_t i, uint8_t e) const {
  IndexRange r = getNode(i);

  size_t j = WPred(r.end, e | 0b1000u);
  if (j >= r.begin) {
    return std::make_optional(j);
  } else {
    j = WPred(r.end, e & 0b0111u);
    if (j >= r.begin)
      return std::make_optional(j);
    else
      return std::nullopt;
  }
}
