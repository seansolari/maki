#pragma once

#include "maki/core/seq/concepts.hpp"
#include "seqan3/alphabet/views/complement.hpp"
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include <algorithm>
#include <cctype>
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

  constexpr std::size_t size() const noexcept { return sequence.size(); }
  constexpr bool empty() const noexcept { return sequence.empty(); }
  constexpr sequence_type const &data() const noexcept { return sequence; }

  void validate() const;
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

  constexpr std::string const &id() const noexcept { return id_; }

  constexpr fragment_container const &fragments() const noexcept {
    return fragments_;
  }
  constexpr fragment_container &fragments() noexcept { return fragments_; }
  constexpr std::size_t numFragments() const noexcept {
    return fragments_.size();
  }
  constexpr bool empty() const noexcept { return fragments_.empty(); }

  std::size_t length() const noexcept;
  std::size_t numKmers(std::size_t k) const noexcept;

  void reserveFragments(std::size_t count) { fragments_.reserve(count); }
  void addFragment(ReadFragment &&fragment);

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
  static FastqDataset fromUnpaired(std::filesystem::path const &path);

  /**
   * Load two FASTQ files in lockstep.
   *
   * The files must contain the same number of records. If validate_ids is
   * true, normalised read identifiers must also agree.
   */
  static FastqDataset fromPairedFiles(std::filesystem::path const &forward_path,
                                      std::filesystem::path const &reverse_path,
                                      bool validate_ids = true);

  /**
   * Load consecutive FASTQ records as read pairs:
   *
   *     record 0 + record 1 -> Read 0
   *     record 2 + record 3 -> Read 1
   *     ...
   *
   * An odd record count is rejected.
   */
  static FastqDataset fromInterleaved(std::filesystem::path const &path,
                                      bool validate_ids = true);

  constexpr size_type size() const noexcept { return reads_.size(); }
  constexpr bool empty() const noexcept { return reads_.empty(); }
  constexpr Read &operator[](size_type index) noexcept { return reads_[index]; }
  constexpr Read const &operator[](size_type index) const noexcept {
    return reads_[index];
  }
  constexpr Read &at(size_type index) { return reads_.at(index); }
  constexpr Read const &at(size_type index) const { return reads_.at(index); }
  constexpr iterator begin() noexcept { return reads_.begin(); }
  constexpr iterator end() noexcept { return reads_.end(); }
  constexpr const_iterator begin() const noexcept { return reads_.begin(); }
  constexpr const_iterator end() const noexcept { return reads_.end(); }
  constexpr const_iterator cbegin() const noexcept { return reads_.cbegin(); }
  constexpr const_iterator cend() const noexcept { return reads_.cend(); }
  constexpr std::span<Read> reads() noexcept { return reads_; }
  constexpr std::span<Read const> reads() const noexcept { return reads_; }
  constexpr FastqMetadata const &metadata() const noexcept { return metadata_; }

private:
  template <typename id_range_t>
  static std::string normaliseReadId(id_range_t const &id_range) {
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
  static void appendRecordFragments(Read &destination,
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
      char const base = normaliseBase(seqan3::to_char(*sequence_it));

      if (isCanonicalBase(base)) {
        current.sequence.push_back(seqan3::dna4{}.assign_char(base));

        current.qualities.push_back(*quality_it);
      } else {
        ++discarded_base_count;

        if (!current.empty()) {
          destination.addFragment(std::move(current));
          current = ReadFragment{};
        }
      }
    }

    if (sequence_it != sequence_end || quality_it != quality_end) {
      throw FastqDatasetError{"FASTQ sequence and quality lengths differ."};
    }

    if (!current.empty())
      destination.addFragment(std::move(current));
  }

  constexpr static char normaliseBase(char base) noexcept {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
  }

  constexpr static bool isCanonicalBase(char base) noexcept {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
  }

  constexpr static bool isCompressedPath(std::filesystem::path const &path) {
    std::string extension = path.extension().string();

    std::ranges::transform(extension, extension.begin(), [](unsigned char c) {
      return static_cast<char>(std::tolower(c));
    });

    return extension == ".gz" || extension == ".bgzf" || extension == ".bz2";
  }

  container_type reads_{};
  FastqMetadata metadata_{};
};

namespace detail {

constexpr auto flatten_fragments =
    std::views::transform([](const Read &read) { return read.fragments(); }) |
    std::views::join;

constexpr auto view_forward_reads =
    std::views::transform([](const ReadFragment &fragment) {
      return SequenceFragment(std::views::all(fragment.sequence), 0, true);
    });

constexpr auto view_reverse_reads =
    std::views::transform([](const ReadFragment &fragment) {
      return SequenceFragment(std::views::all(fragment.sequence) |
                                  std::views::reverse |
                                  seqan3::views::complement,
                              0, true);
    });

template <std::ranges::random_access_range R>
using forward_reads_view_type =
    decltype(std::declval<R>() | flatten_fragments | view_forward_reads);

template <std::ranges::random_access_range R>
using reverse_reads_view_type =
    decltype(std::declval<R>() | flatten_fragments | view_reverse_reads);

} // namespace detail

struct FastqDatasetRange {
  std::span<const Read> reads;

  /**
   * Number of fragments inserted in the range. Does not include reverse
   * complement.
   */
  std::size_t numStrandFragments() const;

  std::size_t numTerminals(std::size_t k) const;
  std::size_t numKmers(std::size_t k) const;

  using forward_view_type =
      detail::forward_reads_view_type<const decltype(reads) &>;
  using reverse_view_type =
      detail::reverse_reads_view_type<const decltype(reads) &>;

  constexpr forward_view_type forwardSequences() const {
    return reads | detail::flatten_fragments | detail::view_forward_reads;
  }
  constexpr reverse_view_type reverseSequences() const {
    return reads | detail::flatten_fragments | detail::view_reverse_reads;
  }
  constexpr forward_view_type forwardTerminals() const {
    return forwardSequences();
  }
  constexpr reverse_view_type reverseTerminals() const {
    return reverseSequences();
  }
  constexpr forward_view_type
  forwardFragments([[maybe_unused]] std::size_t k) const {
    return forwardSequences();
  }
  constexpr reverse_view_type
  reverseFragments([[maybe_unused]] std::size_t k) const {
    return reverseSequences();
  }
};

static_assert(stranded_sequence_container_like<FastqDatasetRange>);

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
class FastqDatasetChunkView {
public:
  FastqDatasetChunkView(const FastqDataset &dataset, std::size_t chunk_count);

  constexpr std::size_t chunk_count() const noexcept { return chunk_count_; }
  constexpr std::size_t size() const noexcept { return chunk_count_; }
  constexpr bool empty() const noexcept { return false; }

  FastqDatasetRange range(std::size_t chunk_index) const;
  FastqDatasetRange operator[](std::size_t chunk_index) const;

private:
  /**
   * Overflow-safe equivalent of:
   *
   *     floor(index * dataset_size / chunk_count)
   *
   * It also distributes the remainder among the earliest chunks.
   */
  std::size_t boundary(std::size_t index) const noexcept;

  const FastqDataset &dataset_;
  std::size_t chunk_count_;
};

} // namespace reads
