
#pragma once

#include "maki/build/kmers/buffers/nt_encoding.hpp"
#include "maki/build/kmers/buffers/terminals.hpp"
#include "maki/core/seq/seq_concepts.hpp"
#include "maki/core/seq/seq_io.hpp"
#include "archive/sink_manager.hpp"

template <typename Container> struct KmerBuffers {
  template <class... Args>
  KmerBuffers(Args& ...args,
              const std::vector<const SequenceContainer*> &seqs,
              const std::vector<SuffixTable> &blocks)
    : data_(args...), temp_(args...), seqs_(seqs), blocks_(blocks) {}

    void collect(ShortSuffix s_) {
      std::size_t size = blocks_.back()[s_];
      data_.resize(size);
      temp_.resize(size);
      data_.fill(seqs_, blocks_, s_);
      if (size > 1)
        data_.sort(&temp_);
    }

    Container& data() { return data_; }

protected:
  Container data_, temp_;
  const std::vector<const SequenceContainer*> &seqs_;
  const std::vector<SuffixTable> &blocks_;
};

template <typename Container, class... Sinks> struct AnalyseSuffix {
  using Buffers = KmerBuffers<Container>;
  using Bundle = ChunkBundleT<Sinks...>;
  
  std::unique_ptr<Bundle> operator()(uint64_t suffixRank) const {
    auto sfx = ShortSuffix::fromRank(suffixRank);
    if (sfx.size() == s_) {
      auto buf = _getbuffer();
      buf->collect(sfx);
      return _flush(suffixRank, buf->data(), terminals_->endsWith(sfx), parsing::dna4ToDna5(sfx.msb()));
    } else {
      return _flush(suffixRank, terminals_.retrieve(sfx));
    }
  }

  std::unique_ptr<Buffers> _getbuffer() const {

  }

  std::unique_ptr<Bundle> _getbundle(uint64_t id) const {
    auto bnd = pool_->acquire();
    bnd->id = id;
    return bnd;
  }

  void _flush(uint64_t rnk, Container &kmers, TerminalRange tmls, uint64_t msb) const {
    auto bnd = _getbundle(rnk);
  }

  void _flush(uint64_t rnk, TerminalRange tmls) const {
    auto bnd = _getbundle(rnk);

  }

protected:
  const std::size_t s_;
  TerminalRange terminals_;
  std::shared_ptr<BundlePool<Sinks...>> pool_;
};
