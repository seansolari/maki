
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
