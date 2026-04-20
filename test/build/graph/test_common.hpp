#include "maki/build/io/fasta.hpp"
#include "maki/core/graph/cdbg.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <sdsl/construct.hpp>
#include <type_traits>

using ::testing::ContainerEq;
using ::testing::ElementsAreArray;

inline void ParseFastaToGenome(Dna4Genome &genome, const std::string &file,
                               Colours &c, std::size_t k) {
  zstr::ifstream zis(file);
  std::string fastaData(std::istreambuf_iterator<char>(zis), {});

  // process stream data
  std::istringstream fs(fastaData);
  parseFastaStream(genome, fs, c, k);
}

inline wavelet_matrix EdgesToWaveletMatrix(const sdsl::int_vector<4> &edges) {
  wavelet_matrix wm;
  sdsl::construct_im(wm, edges);
  return wm;
}

inline std::pair<std::vector<uint64_t>, std::vector<uint64_t>>
UnpackColours(const ColouredGraph &g) {
  std::pair<std::vector<uint64_t>, std::vector<uint64_t>> ca{};

  std::size_t C = g.carch->chunk_count();
  for (std::size_t c = 0; c < C; ++c) {
    auto view = g.carch->view(c);
    std::size_t L = view.size();
    for (std::size_t i = 0; i < L; ++i) {
      auto cols = g.cmap.colours(view.get(i));
      std::visit(
          [&](auto &&arg) {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, colour_t>) {
              // Push a single colour
              ca.first.push_back(arg);
              ca.second.push_back(1);
            } else if constexpr (std::is_same_v<T,
                                                ColourRegistry::vector_ref>) {
              // Push multiple colours
              const auto &colvec = arg.get();
              ca.first.insert(ca.first.end(), colvec.cbegin(), colvec.cend());
              ca.second.push_back(1);
              ca.second.insert(ca.second.end(), colvec.size() - 1, 0u);
            }
          },
          cols);
    }
  }

  return ca;
}