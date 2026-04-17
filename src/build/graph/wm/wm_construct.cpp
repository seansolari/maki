
#include "maki/build/graph/wm/wm_construct.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/utils/logging.hpp"
#include "maki/core/utils/tempfile.hpp"

#include <filesystem>
#include <sdsl/construct.hpp>
#include <sdsl/int_vector_mapper.hpp>

namespace fs = std::filesystem;

void initW(const std::string &file, const std::string &out,
           const std::string &tmp_base) {

  LOG_INFO() << "Starting wavelet matrix construction";
  LOG_INFO() << "Input file: " << file;
  LOG_INFO() << "Output file: " << out;
  LOG_DEBUG() << "Temporary base directory: " << tmp_base;

  // ---------------------------------------------------------------------
  // Temporary directory
  // ---------------------------------------------------------------------
  auto tmp_dir = tempio::create_temporary_directory(tmp_base);

  LOG_DEBUG() << "Created temporary directory for WM construction: " << tmp_dir;

  // ---------------------------------------------------------------------
  // Input buffering
  // ---------------------------------------------------------------------
  sdsl::int_vector_buffer<4> iv(file, std::ios::in);

  const auto n = iv.size();
  LOG_INFO() << "Buffered integer vector with " << n << " entries";

  if (n == 0) {
    LOG_WARN() << "Input edge array is empty — resulting wavelet matrix "
               << "will also be empty";
  }

  // ---------------------------------------------------------------------
  // Wavelet matrix construction
  // ---------------------------------------------------------------------
  LOG_INFO() << "Constructing wavelet matrix";
  wavelet_matrix wm(iv.begin(), iv.end(), tmp_dir);

  LOG_DEBUG() << "Wavelet matrix constructed "
              << "(sigma = " << wm.sigma << ")";

  if (wm.sigma == 0) {
    LOG_WARN() << "Wavelet matrix sigma is zero (degenerate alphabet)";
  }

  // ---------------------------------------------------------------------
  // Persist to disk
  // ---------------------------------------------------------------------
  LOG_INFO() << "Serialising wavelet matrix to disk";
  store_to_file(wm, out);

  if (!fs::exists(out)) {
    LOG_ERROR() << "Wavelet matrix output file was not created: " << out;
  } else {
    LOG_INFO() << "Wavelet matrix written successfully";
  }

  // ---------------------------------------------------------------------
  // Cleanup
  // ---------------------------------------------------------------------
  LOG_DEBUG() << "Removing temporary directory: " << tmp_dir;
  std::error_code ec;
  fs::remove_all(tmp_dir, ec);

  if (ec) {
    LOG_WARN() << "Failed to completely remove temporary directory " << tmp_dir
               << " : " << ec.message();
  } else {
    LOG_DEBUG() << "Temporary directory removed successfully";
  }

  LOG_INFO() << "Wavelet matrix construction complete";
}