
#include "maki/build/io/fastq.hpp"
#include <numeric>

namespace reads {

void ReadFragment::validate() const {
  if (sequence.size() != qualities.size()) {
    throw FastqDatasetError{
        "ReadFragment sequence and quality lengths differ."};
  }
}

std::size_t Read::length() const noexcept {
  return std::accumulate(fragments_.cbegin(), fragments_.cend(), std::size_t{0},
                         [](std::size_t total, const ReadFragment &fragment) {
                           return total + fragment.size();
                         });
}

std::size_t Read::numKmers(std::size_t k) const noexcept {
  return std::accumulate(fragments_.cbegin(), fragments_.cend(), std::size_t{0},
                         [k](std::size_t total, const ReadFragment &fragment) {
                           return total + fragment.size() - k + 1;
                         });
}

void Read::addFragment(ReadFragment &&fragment) {
  fragment.validate();

  if (!fragment.empty())
    fragments_.push_back(std::move(fragment));
}

FastqDataset FastqDataset::fromUnpaired(std::filesystem::path const &path) {
  FastqDataset dataset;
  dataset.metadata_.layout = FastqLayout::unpaired;
  dataset.metadata_.first_path = path;
  dataset.metadata_.compressed_input = isCompressedPath(path);

  seqan3::sequence_file_input input{path};

  for (auto &record : input) {
    Read read{normaliseReadId(record.id())};

    appendRecordFragments(read, record.sequence(), record.base_qualities(),
                          dataset.metadata_.discarded_base_count);

    dataset.reads_.push_back(std::move(read));
    ++dataset.metadata_.input_record_count;
  }

  dataset.metadata_.logical_read_count = dataset.reads_.size();
  return dataset;
}

FastqDataset
FastqDataset::fromPairedFiles(std::filesystem::path const &forward_path,
                              std::filesystem::path const &reverse_path,
                              bool validate_ids) {
  FastqDataset dataset;
  dataset.metadata_.layout = FastqLayout::paired_files;
  dataset.metadata_.first_path = forward_path;
  dataset.metadata_.second_path = reverse_path;

  dataset.metadata_.compressed_input =
      isCompressedPath(forward_path) || isCompressedPath(reverse_path);

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

    std::string forward_id = normaliseReadId(forward_record.id());

    std::string reverse_id = normaliseReadId(reverse_record.id());

    if (validate_ids && forward_id != reverse_id) {
      throw FastqDatasetError{"Paired FASTQ identifier mismatch at pair " +
                              std::to_string(pair_index) + ": '" + forward_id +
                              "' versus '" + reverse_id + "'."};
    }

    Read read{std::move(forward_id)};

    appendRecordFragments(read, forward_record.sequence(),
                          forward_record.base_qualities(),
                          dataset.metadata_.discarded_base_count);

    appendRecordFragments(read, reverse_record.sequence(),
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

FastqDataset FastqDataset::fromInterleaved(std::filesystem::path const &path,
                                           bool validate_ids) {
  FastqDataset dataset;
  dataset.metadata_.layout = FastqLayout::interleaved;
  dataset.metadata_.first_path = path;
  dataset.metadata_.compressed_input = isCompressedPath(path);

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
    std::string forward_id = normaliseReadId(forward_record.id());

    Read read{forward_id};

    appendRecordFragments(read, forward_record.sequence(),
                          forward_record.base_qualities(),
                          dataset.metadata_.discarded_base_count);

    auto &reverse_record = *it;

    std::string reverse_id = normaliseReadId(reverse_record.id());

    if (validate_ids && forward_id != reverse_id) {
      throw FastqDatasetError{"Interleaved FASTQ identifier mismatch at pair " +
                              std::to_string(pair_index) + ": '" + forward_id +
                              "' versus '" + reverse_id + "'."};
    }

    appendRecordFragments(read, reverse_record.sequence(),
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

std::size_t FastqDatasetRange::numStrandFragments() const {
  return std::accumulate(
      reads.cbegin(), reads.cend(), (std::size_t)0,
      [](std::size_t total, const Read &read) -> std::size_t {
        return total + read.numFragments();
      });
}

std::size_t FastqDatasetRange::numTerminals(std::size_t k) const {
  return std::size_t{2} * k * numStrandFragments();
}

std::size_t FastqDatasetRange::numKmers(std::size_t k) const {
  return std::size_t{2} * std::accumulate(reads.cbegin(), reads.cend(),
                                          (std::size_t)0,
                                          [k](std::size_t total,
                                              const Read &read) -> std::size_t {
                                            return total + read.numKmers(k);
                                          });
}

FastqDatasetChunkView::FastqDatasetChunkView(const FastqDataset &dataset,
                                             std::size_t chunk_count)
    : dataset_{dataset}, chunk_count_{chunk_count} {
  if (chunk_count_ == 0) {
    throw std::invalid_argument{
        "The number of chunks must be greater than zero."};
  }
}

FastqDatasetRange FastqDatasetChunkView::range(std::size_t chunk_index) const {
  if (chunk_index >= chunk_count_) {
    throw std::out_of_range{"FASTQ dataset chunk index is out of range."};
  }

  std::size_t const begin_index = boundary(chunk_index);
  std::size_t const end_index = boundary(chunk_index + 1);

  auto dataset_span = dataset_.reads();
  return {dataset_span.subspan(begin_index, end_index - begin_index)};
}

FastqDatasetRange
FastqDatasetChunkView::operator[](std::size_t chunk_index) const {
  return range(chunk_index);
}

std::size_t FastqDatasetChunkView::boundary(std::size_t index) const noexcept {
  std::size_t const element_count = dataset_.size();
  std::size_t const quotient = element_count / chunk_count_;
  std::size_t const remainder = element_count % chunk_count_;

  return index * quotient + std::min(index, remainder);
}

} // namespace reads