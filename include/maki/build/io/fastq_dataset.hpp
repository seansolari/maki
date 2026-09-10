#pragma once

#include "maki/core/seq/concepts.hpp"
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include <algorithm>
#include <cctype>
#include <concepts>
#include <cstddef>
#include <filesystem>
#include <iterator>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace reads {

enum class FastqLayout { unpaired, paired_files, interleaved };

class FastqDatasetError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

/**
 * One contiguous A/C/G/T region from an input read.
 *
 * sequence[i] and qualities[i] always refer to the same base.
 */
struct ReadFragment {
  using sequence_type = seqan3::bitpacked_sequence<seqan3::dna4>;

  using quality_type = seqan3::bitpacked_sequence<seqan3::phred42>;

  sequence_type sequence{};
  quality_type qualities{};

  [[nodiscard]]
  std::size_t size() const noexcept {
    return sequence.size();
  }

  [[nodiscard]]
  bool empty() const noexcept {
    return sequence.empty();
  }

  void validate() const {
    if (sequence.size() != qualities.size()) {
      throw FastqDatasetError{
          "ReadFragment sequence and quality lengths differ."};
    }
  }
};

/**
 * A logical read.
 *
 * For unpaired input, this contains fragments from one FASTQ record.
 * For paired-end input, this contains fragments from both mates.
 *
 * Mate provenance is deliberately not stored.
 */
class Read {
public:
  using fragment_type = ReadFragment;
  using fragment_container = std::vector<fragment_type>;

  Read() = default;

  explicit Read(std::string id) : id_{std::move(id)} {}

  [[nodiscard]]
  std::string const &id() const noexcept {
    return id_;
  }

  [[nodiscard]]
  fragment_container const &fragments() const noexcept {
    return fragments_;
  }

  [[nodiscard]]
  fragment_container &fragments() noexcept {
    return fragments_;
  }

  [[nodiscard]]
  std::size_t fragment_count() const noexcept {
    return fragments_.size();
  }

  /**
   * Number of retained canonical bases across all fragments.
   */
  [[nodiscard]]
  std::size_t retained_base_count() const noexcept {
    std::size_t result{};

    for (ReadFragment const &fragment : fragments_)
      result += fragment.size();

    return result;
  }

  [[nodiscard]]
  bool empty() const noexcept {
    return fragments_.empty();
  }

  void reserve_fragments(std::size_t count) { fragments_.reserve(count); }

  void add_fragment(ReadFragment fragment) {
    fragment.validate();

    if (!fragment.empty())
      fragments_.push_back(std::move(fragment));
  }

private:
  std::string id_{};
  fragment_container fragments_{};
};

/**
 * Metadata describing how a FastqDataset was loaded.
 */
struct FastqMetadata {
  FastqLayout layout{FastqLayout::unpaired};

  std::filesystem::path first_path{};
  std::filesystem::path second_path{};

  std::size_t input_record_count{};
  std::size_t logical_read_count{};

  /**
   * Number of bases discarded because they were represented by the
   * non-canonical dna5 symbol N.
   */
  std::size_t discarded_base_count{};

  bool compressed_input{};
};

/**
 * An owning, random-access collection of logical reads.
 *
 * Loading is eager. This makes the resulting dataset safe to partition into
 * independent spans after construction, provided that its size is not changed
 * while worker threads are using it.
 */
class FastqDataset {
public:
  using value_type = Read;
  using container_type = std::vector<value_type>;
  using size_type = container_type::size_type;
  using iterator = container_type::iterator;
  using const_iterator = container_type::const_iterator;

  FastqDataset() = default;

  /**
   * Load one FASTQ file as unpaired reads.
   */
  [[nodiscard]]
  static FastqDataset from_unpaired(std::filesystem::path const &path) {
    FastqDataset dataset;
    dataset.metadata_.layout = FastqLayout::unpaired;
    dataset.metadata_.first_path = path;
    dataset.metadata_.compressed_input = is_compressed_path(path);

    seqan3::sequence_file_input input{path};

    for (auto &record : input) {
      Read read{normalise_read_id(record.id())};

      append_record_fragments(read, record.sequence(), record.base_qualities(),
                              dataset.metadata_.discarded_base_count);

      dataset.reads_.push_back(std::move(read));
      ++dataset.metadata_.input_record_count;
    }

    dataset.metadata_.logical_read_count = dataset.reads_.size();
    return dataset;
  }

  /**
   * Load two FASTQ files in lockstep.
   *
   * The files must contain the same number of records. If validate_ids is
   * true, normalised read identifiers must also agree.
   */
  [[nodiscard]]
  static FastqDataset
  from_paired_files(std::filesystem::path const &forward_path,
                    std::filesystem::path const &reverse_path,
                    bool validate_ids = true) {
    FastqDataset dataset;
    dataset.metadata_.layout = FastqLayout::paired_files;
    dataset.metadata_.first_path = forward_path;
    dataset.metadata_.second_path = reverse_path;

    dataset.metadata_.compressed_input =
        is_compressed_path(forward_path) || is_compressed_path(reverse_path);

    seqan3::sequence_file_input forward_input{forward_path};
    seqan3::sequence_file_input reverse_input{reverse_path};

    auto forward_it = forward_input.begin();
    auto reverse_it = reverse_input.begin();

    auto const forward_end = forward_input.end();
    auto const reverse_end = reverse_input.end();

    std::size_t pair_index{};

    while (forward_it != forward_end && reverse_it != reverse_end) {
      auto &forward_record = *forward_it;
      auto &reverse_record = *reverse_it;

      std::string forward_id = normalise_read_id(forward_record.id());

      std::string reverse_id = normalise_read_id(reverse_record.id());

      if (validate_ids && forward_id != reverse_id) {
        throw FastqDatasetError{"Paired FASTQ identifier mismatch at pair " +
                                std::to_string(pair_index) + ": '" +
                                forward_id + "' versus '" + reverse_id + "'."};
      }

      Read read{std::move(forward_id)};

      append_record_fragments(read, forward_record.sequence(),
                              forward_record.base_qualities(),
                              dataset.metadata_.discarded_base_count);

      append_record_fragments(read, reverse_record.sequence(),
                              reverse_record.base_qualities(),
                              dataset.metadata_.discarded_base_count);

      dataset.reads_.push_back(std::move(read));

      ++forward_it;
      ++reverse_it;
      ++pair_index;
      dataset.metadata_.input_record_count += 2;
    }

    if (forward_it != forward_end || reverse_it != reverse_end) {
      throw FastqDatasetError{
          "Paired FASTQ files contain different numbers of records."};
    }

    dataset.metadata_.logical_read_count = dataset.reads_.size();
    return dataset;
  }

  /**
   * Load consecutive FASTQ records as read pairs:
   *
   *     record 0 + record 1 -> Read 0
   *     record 2 + record 3 -> Read 1
   *     ...
   *
   * An odd record count is rejected.
   */
  [[nodiscard]]
  static FastqDataset from_interleaved(std::filesystem::path const &path,
                                       bool validate_ids = true) {
    FastqDataset dataset;
    dataset.metadata_.layout = FastqLayout::interleaved;
    dataset.metadata_.first_path = path;
    dataset.metadata_.compressed_input = is_compressed_path(path);

    seqan3::sequence_file_input input{path};

    auto it = input.begin();
    auto const end = input.end();

    std::size_t pair_index{};

    while (it != end) {
      auto &forward_record = *it;
      ++it;

      if (it == end) {
        throw FastqDatasetError{
            "Interleaved FASTQ contains an odd number of records."};
      }

      /*
       * A SeqAn input record is backed by the input object and is
       * overwritten when the iterator is advanced. Therefore copy the
       * fields from the first mate before processing the second mate.
       */
      std::string forward_id = normalise_read_id(forward_record.id());

      Read read{forward_id};

      append_record_fragments(read, forward_record.sequence(),
                              forward_record.base_qualities(),
                              dataset.metadata_.discarded_base_count);

      auto &reverse_record = *it;

      std::string reverse_id = normalise_read_id(reverse_record.id());

      if (validate_ids && forward_id != reverse_id) {
        throw FastqDatasetError{
            "Interleaved FASTQ identifier mismatch at pair " +
            std::to_string(pair_index) + ": '" + forward_id + "' versus '" +
            reverse_id + "'."};
      }

      append_record_fragments(read, reverse_record.sequence(),
                              reverse_record.base_qualities(),
                              dataset.metadata_.discarded_base_count);

      dataset.reads_.push_back(std::move(read));

      ++it;
      ++pair_index;
      dataset.metadata_.input_record_count += 2;
    }

    dataset.metadata_.logical_read_count = dataset.reads_.size();
    return dataset;
  }

  [[nodiscard]]
  size_type size() const noexcept {
    return reads_.size();
  }

  [[nodiscard]]
  bool empty() const noexcept {
    return reads_.empty();
  }

  [[nodiscard]]
  Read &operator[](size_type index) noexcept {
    return reads_[index];
  }

  [[nodiscard]]
  Read const &operator[](size_type index) const noexcept {
    return reads_[index];
  }

  [[nodiscard]]
  Read &at(size_type index) {
    return reads_.at(index);
  }

  [[nodiscard]]
  Read const &at(size_type index) const {
    return reads_.at(index);
  }

  [[nodiscard]]
  iterator begin() noexcept {
    return reads_.begin();
  }

  [[nodiscard]]
  iterator end() noexcept {
    return reads_.end();
  }

  [[nodiscard]]
  const_iterator begin() const noexcept {
    return reads_.begin();
  }

  [[nodiscard]]
  const_iterator end() const noexcept {
    return reads_.end();
  }

  [[nodiscard]]
  const_iterator cbegin() const noexcept {
    return reads_.cbegin();
  }

  [[nodiscard]]
  const_iterator cend() const noexcept {
    return reads_.cend();
  }

  [[nodiscard]]
  std::span<Read> reads() noexcept {
    return reads_;
  }

  [[nodiscard]]
  std::span<Read const> reads() const noexcept {
    return reads_;
  }

  [[nodiscard]]
  FastqMetadata const &metadata() const noexcept {
    return metadata_;
  }

private:
  template <typename id_range_t>
  [[nodiscard]]
  static std::string normalise_read_id(id_range_t const &id_range) {
    std::string id;

    for (auto const symbol : id_range)
      id.push_back(seqan3::to_char(symbol));

    /*
     * FASTQ identifiers can contain a description after whitespace.
     */
    if (auto const position = id.find_first_of(" \t");
        position != std::string::npos) {
      id.resize(position);
    }

    /*
     * Common paired-read forms include:
     *
     *   sample/1 and sample/2
     *   sample.1 and sample.2
     *
     * Casava-style "sample 1:N:..." identifiers have already been
     * handled by truncating at whitespace.
     */
    if (id.size() >= 2) {
      char const separator = id[id.size() - 2];
      char const mate = id[id.size() - 1];

      if ((separator == '/' || separator == '.') &&
          (mate == '1' || mate == '2')) {
        id.resize(id.size() - 2);
      }
    }

    return id;
  }

  template <typename sequence_range_t, typename quality_range_t>
  static void append_record_fragments(Read &destination,
                                      sequence_range_t const &sequence,
                                      quality_range_t const &qualities,
                                      std::size_t &discarded_base_count) {
    auto sequence_it = std::ranges::begin(sequence);
    auto quality_it = std::ranges::begin(qualities);

    auto const sequence_end = std::ranges::end(sequence);
    auto const quality_end = std::ranges::end(qualities);

    ReadFragment current;

    for (; sequence_it != sequence_end && quality_it != quality_end;
         ++sequence_it, ++quality_it) {
      char const base = normalise_base(seqan3::to_char(*sequence_it));

      if (is_canonical_base(base)) {
        current.sequence.push_back(seqan3::dna4{}.assign_char(base));

        current.qualities.push_back(*quality_it);
      } else {
        ++discarded_base_count;

        if (!current.empty()) {
          destination.add_fragment(std::move(current));
          current = ReadFragment{};
        }
      }
    }

    if (sequence_it != sequence_end || quality_it != quality_end) {
      throw FastqDatasetError{"FASTQ sequence and quality lengths differ."};
    }

    if (!current.empty())
      destination.add_fragment(std::move(current));
  }

  [[nodiscard]]
  static char normalise_base(char base) noexcept {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
  }

  [[nodiscard]]
  static bool is_canonical_base(char base) noexcept {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
  }

  [[nodiscard]]
  static bool is_compressed_path(std::filesystem::path const &path) {
    std::string extension = path.extension().string();

    std::ranges::transform(extension, extension.begin(), [](unsigned char c) {
      return static_cast<char>(std::tolower(c));
    });

    return extension == ".gz" || extension == ".bgzf" || extension == ".bz2";
  }

  container_type reads_{};
  FastqMetadata metadata_{};
};

template <typename read_t> struct FastqDatasetRange : public SequenceContainer {
  std::span<read_t> range;

  std::size_t numFragments() const {
    return std::accumulate(
      data.cbegin(), data.cend(), (std::size_t)0,
      [](std::size_t total, const ReadVector &v) -> std::size_t {
        return total + detail::numFragments(v);
      });
  }

  virtual std::size_t numTerminals(std::size_t k) const override final {
    return (std::size_t)k * numFragments();
  }

  virtual poly_input_range<SequenceFragment> terminals() const override final {
    return poly_input_range<SequenceFragment>(
      data | std::views::join |
      std::views::transform([](const Read &r) { return r.view(); }) |
      std::views::join | std::views::transform([](const Dna4Sequence &seq) {
        return SequenceFragment(seq.cbegin(), seq.cend(), 0, true);
      }));
  }

  virtual std::size_t numKmers(std::size_t k) const override final {
    return std::accumulate(
      data.cbegin(), data.cend(), (std::size_t)0,
      [k](std::size_t total, const ReadVector &v) -> std::size_t {
        return total + detail::numKmers(v, k);
      });
  }

  virtual poly_input_range<SequenceFragment>
  fragments(std::size_t k) const override final {
    return poly_input_range<SequenceFragment>(
      data | std::views::join |
      std::views::transform([](const Read &r) { return r.view(); }) |
      std::views::join | std::views::transform([](const Dna4Sequence &seq) {
        return SequenceFragment(seq.cbegin(), seq.cend(), 0, true);
      }));
  }
};

/**
 * A lightweight deterministic partition over a FastqDataset.
 *
 * No per-chunk pointers, iterators, or spans are stored. Each range is
 * calculated from:
 *
 *     floor(i * dataset_size / chunk_count)
 *
 * The remainder is distributed across chunks while ensuring that chunk sizes
 * differ by at most one.
 *
 * The dataset must not be resized while this view or any returned spans are
 * being used.
 */
template <typename dataset_t> class BasicFastqDatasetChunkView {
private:
  using raw_dataset_type = std::remove_reference_t<dataset_t>;

  static_assert(
      std::same_as<std::remove_const_t<raw_dataset_type>, FastqDataset>,
      "BasicFastqDatasetChunkView requires FastqDataset.");

public:
  using read_type =
      std::conditional_t<std::is_const_v<raw_dataset_type>, Read const, Read>;

  using range_type = FastqDatasetRange<read_type>;

  BasicFastqDatasetChunkView(dataset_t &dataset, std::size_t chunk_count)
      : dataset_{dataset}, chunk_count_{chunk_count} {
    if (chunk_count_ == 0) {
      throw std::invalid_argument{
          "The number of chunks must be greater than zero."};
    }
  }

  [[nodiscard]]
  std::size_t chunk_count() const noexcept {
    return chunk_count_;
  }

  [[nodiscard]]
  std::size_t size() const noexcept {
    return chunk_count_;
  }

  [[nodiscard]]
  bool empty() const noexcept {
    return false;
  }

  [[nodiscard]]
  range_type range(std::size_t chunk_index) const {
    if (chunk_index >= chunk_count_) {
      throw std::out_of_range{"FASTQ dataset chunk index is out of range."};
    }

    std::size_t const begin_index = boundary(chunk_index);

    std::size_t const end_index = boundary(chunk_index + 1);

    auto dataset_span = dataset_.reads();

    return dataset_span.subspan(begin_index, end_index - begin_index);
  }

  [[nodiscard]]
  range_type operator[](std::size_t chunk_index) const {
    return range(chunk_index);
  }

private:
  /**
   * Overflow-safe equivalent of:
   *
   *     floor(index * dataset_size / chunk_count)
   *
   * It also distributes the remainder among the earliest chunks.
   */
  [[nodiscard]]
  std::size_t boundary(std::size_t index) const noexcept {
    std::size_t const element_count = dataset_.size();
    std::size_t const quotient = element_count / chunk_count_;
    std::size_t const remainder = element_count % chunk_count_;

    return index * quotient + std::min(index, remainder);
  }

  dataset_t &dataset_;
  std::size_t chunk_count_;
};

using FastqDatasetChunkView = BasicFastqDatasetChunkView<FastqDataset>;

using ConstFastqDatasetChunkView =
    BasicFastqDatasetChunkView<FastqDataset const>;

} // namespace reads
