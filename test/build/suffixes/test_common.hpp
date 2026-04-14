
#include "maki/build/io/fasta.hpp"
#include "maki/core/utils/tempfile.hpp"
#include <bitset>
#include <fstream>
#include <gmock/gmock.h>

using ::testing::Eq;
using ::testing::Value;

MATCHER_P(MatrixEq, expected_wrapper, "ndim::matrix") {
  auto &expected = expected_wrapper.get();
  if (arg.size() != expected.size()) {
    *result_listener << "arg.size() != expected.size() ";
    *result_listener << arg.size() << " vs " << expected.size();
    return false;
  }
  if (arg.rows() != expected.rows()) {
    *result_listener << "arg.rows() != expected.rows() ";
    *result_listener << arg.rows() << " vs " << expected.rows();
    return false;
  }
  if (arg.cols() != expected.cols()) {
    *result_listener << "arg.cols() != expected.cols() ";
    *result_listener << arg.cols() << " vs " << expected.cols();
    return false;
  }
  for (size_t r = 0; r < expected.rows(); ++r) {
    auto arg_row = arg[r];
    auto expected_row = expected[r];
    for (size_t c = 0; c < expected.cols(); ++c) {
      if (!Value(arg_row[c], Eq(expected_row[c]))) {
        *result_listener << "element[" << r << ", " << c << "] mismatch ";
        *result_listener << "0b" << std::bitset<8>(arg_row[c]) << " vs "
                         << "0b" << std::bitset<8>(expected_row[c]);
        return false;
      }
    }
  }
  return true;
}

inline std::string MakeTempPath(const char *ext) {
  return tempio::create_temporary_file(std::filesystem::temp_directory_path(),
                                       ext);
}

inline int writeToFna(const std::string &data, const std::string &fna) {
  std::ofstream output_file(fna);

  // Check if the file was successfully opened
  if (output_file.is_open()) {
    // Write the string data to the file using the insertion operator (<<)
    output_file << data;

    // Close the file
    output_file.close();
  } else {
    std::cerr << "Error: Unable to open the file." << std::endl;
    return 1; // Return an error code
  }

  return 0;
}

inline void ParseFastaToGenome(Dna4Genome &genome, const std::string &file, Colours &c, std::size_t k) {
  zstr::ifstream zis(file);
  std::string fastaData(std::istreambuf_iterator<char>(zis), {});

  // process stream data
  std::istringstream fs(fastaData);
  parseFastaStream(genome, fs, c, k);
}
