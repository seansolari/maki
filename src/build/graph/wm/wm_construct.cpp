
#include "maki/build/graph/wm/wm_construct.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/utils/tempfile.hpp"
#include <filesystem>
#include <sdsl/construct.hpp>
#include <sdsl/int_vector_mapper.hpp>

void initW(const std::string &file, const std::string &out, const std::string &tmp_base) {
  auto tmp_dir = tempio::create_temporary_directory(tmp_base);
  sdsl::int_vector_buffer<4> iv(file, std::ios::in);
  wavelet_matrix wm(iv.begin(), iv.end(), tmp_dir);
  store_to_file(wm, out);
  fs::remove_all(tmp_dir);
}
