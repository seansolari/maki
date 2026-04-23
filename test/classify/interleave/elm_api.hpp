#pragma once
#include <memory>
#include <string>
#include <vector>

#include "dbg_api.hpp"

// Every ELM merge implementation must satisfy this API
template <typename Impl> struct ELMTestAPI {
  using Graph = typename Impl::Graph;

  static Graph Build(size_t k, const std::vector<std::string> &kmers) {
    return Impl::BuildFromKmers(k, kmers);
  }

  static Graph Merge(const Graph &a, const Graph &b) {
    return Impl::Merge(a, b);
  }

  static std::unique_ptr<DBGView> View(const Graph &g) { return Impl::View(g); }

  // Optional: refinement tracing
  static bool HasRefinementTrace() { return Impl::HasRefinementTrace(); }

  static auto Trace(const Graph &a, const Graph &b) {
    return Impl::TraceMerge(a, b); // implementation-defined type
  }
};