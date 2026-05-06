#pragma once

#include "maki/build/graph/archive/sink_manager.hpp"
#include "maki/build/graph/interleave_buffers.hpp"
#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/core/seq/concepts.hpp"
#include "maki/core/utils/logging.hpp"
#include <cstddef>
#include <filesystem>

namespace dbg {

struct BuildOptions {
  // Algorithm parameters
  std::size_t kmer_size = 31;
  std::size_t suffix_size = 8;
  // Output parameters
  std::filesystem::path out;
  // Space parameters
  std::size_t pool_size = 16;
  std::size_t reserve_per_chunk = 0;
  std::size_t chunks() const;
};

// -----------------------------------------------------------------------------
// Intermediate data
// -----------------------------------------------------------------------------

template <class MainBufferType> struct Buffers {
  Buffers(std::size_t length, std::size_t width, std::size_t k,
          std::size_t keff)
      : kmers(length, width, k, keff), temp(length, width, k, keff),
        terminals(), b(), t() {}

  /**
   * Collect k-mers with given suffix `s_` from input sequences.
   */
  void collectKmers(const std::vector<const SequenceContainer *> &seqs_,
                    const std::vector<SuffixTable> &blocks_, ShortSuffix s_) {
    std::size_t size = blocks_.back()[s_];

    LOG_DEBUG() << "Collecting k-mers for suffix " << s_.toString()
                << " ; count = " << size;

    kmers.resize(size);
    temp.resize(size);

    kmers.fill(seqs_, blocks_, s_);

    if (size > 1) {
      kmers.sort(&temp);
    }

    adjacentDifference(kmers, b);
  }

  /**
   * Load terminals data.
   */
  void setTerminals(TerminalRange &&data_) {
    const std::size_t n = data_.size();

    if (n == 0) {
      LOG_DEBUG() << "No terminal nodes for this suffix";
    }

    terminals = std::move(data_);
    adjacentDifference(terminals, t);
  }

  MainBufferType kmers, temp;
  TerminalRange terminals;
  sdsl::int_vector<2> b, t;
};

template <class MainBufferType>
struct BufferMaker : public Factory<Buffers<MainBufferType>> {
  BufferMaker(std::size_t length_, std::size_t width_, std::size_t k_,
              std::size_t keff_)
      : Factory<Buffers<MainBufferType>>(), length(length_), width(width_),
        k(k_), keff(keff_) {}

  std::size_t length;
  std::size_t width;
  std::size_t k;
  std::size_t keff;

  inline std::unique_ptr<Buffers<MainBufferType>> obtain() {
    return Factory<Buffers<MainBufferType>>::obtain(length, width, k, keff);
  }
};

// -----------------------------------------------------------------------------
// Pipeline
// -----------------------------------------------------------------------------

template <class MainBufferType, class... Sinks> struct Suffixwise {
  Suffixwise(const std::vector<const SequenceContainer *> &seqs,
             const TerminalRange &terms,
             std::shared_ptr<std::vector<SuffixTable>> &&suffixPlan,
             uint8_t bufferValueSize, std::size_t k, std::size_t s,
             push_summary *str)
      : suffixSize(s), graphStructure(str), sequences(seqs), terminals(terms),
        suffixCounts(std::move(suffixPlan)),
        kmerBuffers(std::make_shared<BufferMaker<MainBufferType>>(
            suffixCounts->back().maxValue(), bufferValueSize, k, k - s)),
        bundles() {
    LOG_INFO() << "Initialised suffix-wise k-mer extractor"
               << " (k=" << k << ", suffix=" << s << ")";

    LOG_INFO() << "Suffix plan maximum value = "
               << suffixCounts->back().maxValue();

    LOG_INFO() << "Terminal ranges total = " << terminals.size();
  }

  // input data
  std::size_t suffixSize;
  push_summary *graphStructure;

  // input buffers
  const std::vector<const SequenceContainer *> &sequences;
  const TerminalRange &terminals;

  // shared auxilliary data
  std::shared_ptr<std::vector<SuffixTable>> suffixCounts;
  std::shared_ptr<BufferMaker<MainBufferType>> kmerBuffers;

  // output buffers
  std::shared_ptr<::BundlePool<Sinks...>> bundles;

  // Initialisation ---

  // Source of data bundles where analysed data is placed
  void setPool(std::shared_ptr<::BundlePool<Sinks...>> &p) { bundles = p; }

  // Obtain a data bundle and tag it
  std::unique_ptr<ChunkBundleT<Sinks...>> getBundle(uint64_t id) const {
    auto bnd = bundles->acquire();
    bnd->id = id;
    return bnd;
  }

  // Finalise chunks ---

  // Account freshly inserted nodes and edges per H1 block
  void pushRegionStructure(push_summary &tkn) const {
    for (std::size_t i = 0; i < 5; ++i)
      std::atomic_ref(graphStructure->F[i]) += tkn.F[i];
    for (std::size_t i = 0; i < 5; ++i)
      std::atomic_ref(graphStructure->C[i]) += tkn.C[i];
  }
};

} // namespace dbg
