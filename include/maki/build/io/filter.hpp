
// ---------------------------------------------------------------------------
// Parsing large FNA files (filters)
// ---------------------------------------------------------------------------

struct RestrictedSequenceFragment final : public SequenceContainer {
  RestrictedSequenceFragment(Dna4SequenceConstIter it, std::size_t begin,
                             std::size_t end, uint64_t id)
      : _it(it), _begin(begin), _end(end), _id(id) {}

  inline std::size_t size() const noexcept { return _end - _begin; }
  inline void setEndToTerminal() { _terminal = true; }
  inline bool endIsTerminal() const noexcept { return _terminal; }

  virtual std::size_t numTerminals(std::size_t k) const override final;
  virtual poly_input_range<SequenceFragment> terminals() const override final;
  virtual std::size_t numKmers(std::size_t k) const override final;
  virtual poly_input_range<SequenceFragment>
  fragments(std::size_t k) const override final;

protected:
  Dna4SequenceConstIter _it;
  std::size_t _begin;
  std::size_t _end;
  uint64_t _id;
  bool _terminal = false;
};

struct ChunkedDna4Genome {
  Dna4Genome genome;
  std::vector<RestrictedSequenceFragment> chunks;

  ChunkedDna4Genome() = default;
  ChunkedDna4Genome(ChunkedDna4Genome &&) = default;
  ChunkedDna4Genome &operator=(ChunkedDna4Genome &&) = default;
  ChunkedDna4Genome(const ChunkedDna4Genome &) = delete;
  ChunkedDna4Genome &operator=(const ChunkedDna4Genome &) = delete;

  void chunk(std::size_t granularity, std::size_t overlap);
};

std::size_t chunks(const std::vector<ChunkedDna4Genome> &genomes);

// ---------------------------------------------------------------------------
// API
// ---------------------------------------------------------------------------

ChunkedDna4Genome parseFilterFNA(const std::string &fastaFile,
                                 std::size_t granularity, std::size_t k);


std::size_t RestrictedSequenceFragment::numTerminals(std::size_t k) const {
  if (_begin == 0)
    return k;
  else
    return 0;
}

poly_input_range<SequenceFragment>
RestrictedSequenceFragment::terminals() const {
  if (_begin == 0)
    return poly_input_range<SequenceFragment>(
        SequenceFragment(_it + _begin, _it + _end, _id, _terminal));
  else
    return poly_input_range<SequenceFragment>::empty_range();
}

std::size_t RestrictedSequenceFragment::numKmers(std::size_t k) const {
  return size() - k + 1;
}

poly_input_range<SequenceFragment>
RestrictedSequenceFragment::fragments([[maybe_unused]] std::size_t k) const {
  return poly_input_range<SequenceFragment>(
      SequenceFragment(_it + _begin, _it + _end, _id, _terminal));
}

void ChunkedDna4Genome::chunk(std::size_t granularity, std::size_t overlap) {
  if (!chunks.empty())
    chunks.clear();

  std::size_t ntChunkSize =
      std::max(genome.medianContigSize() / granularity, overlap);

  for (const Dna4Contig &ctg : genome.contigs) {
    for (const SequenceVector *drn : {&ctg.sequences[0], &ctg.sequences[1]}) {
      for (const AnnotatedSequence &rec : *drn) {
        assert(rec.annotations.size() == 0 &&
               "chunking scheme can only be applied to unannotated sequences");
        assert(rec.sequence.size() >= overlap &&
               "sequence size is less than requested overlap");
        std::size_t i = overlap;
        std::size_t j;
        while (i < rec.sequence.size()) {
          j = std::min(rec.sequence.size(), i + ntChunkSize);
          chunks.emplace_back(rec.sequence.cbegin(), i - overlap, j,
                              rec.nullFeatureId);
          i = j;
        }
        chunks.back().setEndToTerminal();
      }
    }
  }
}

std::size_t chunks(const std::vector<ChunkedDna4Genome> &genomes) {
  return std::accumulate(genomes.begin(), genomes.end(), (std::size_t)0,
                         [](std::size_t rtot, const ChunkedDna4Genome &g) {
                           return rtot + g.chunks.size();
                         });
}

ChunkedDna4Genome parseFilterFNA(const std::string &fastaFile,
                                 std::size_t granularity, std::size_t k) {
  // base input file stream that reads bytes
  zstr::ifstream zis(fastaFile);
  std::string fastaData(std::istreambuf_iterator<char>(zis), {});

  // process stream data
  ChunkedDna4Genome obj;
  std::istringstream fs(fastaData);
  parseFastaStream(obj.genome, fs, k);

  // create genome slices
  obj.chunk(granularity, k);

  return obj;
}
