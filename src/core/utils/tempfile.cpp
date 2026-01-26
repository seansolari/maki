
#include "maki/core/utils/tempfile.hpp"
#include <cstdint>
#include <ios>
#include <random>
#include <sstream>
#include <stdexcept>

fs::path create_temporary_directory(fs::path tmp_dir,
                                    unsigned long long max_tries) {
  unsigned long long i = 0;
  std::random_device dev;
  std::mt19937 prng(dev());
  std::uniform_int_distribution<uint64_t> rand(0);
  fs::path path;
  while (true) {
    std::stringstream ss;
    ss << std::hex << rand(prng);
    path = tmp_dir / ss.str();
    // true if the directory was created.
    if (fs::create_directory(path)) {
      break;
    }
    if (i == max_tries) {
      throw std::runtime_error("could not find non-existing directory");
    }
    i++;
  }
  return path;
}

fs::path create_temporary_directory(unsigned long long max_tries) {
  return create_temporary_directory(fs::temp_directory_path(), max_tries);
}

fs::path create_temporary_file(fs::path const &base, const char *ext,
                               unsigned long long max_tries) {
  unsigned long long i = 0;
  std::random_device dev;
  std::mt19937 prng(dev());
  std::uniform_int_distribution<uint64_t> rand(0);
  fs::path path;
  while (true) {
    std::stringstream ss;
    ss << std::hex << rand(prng) << ext;
    path = base / ss.str();
    // true if the directory was created.
    if (!fs::exists(path)) {
      break;
    }
    if (i == max_tries) {
      throw std::runtime_error("could not find non-existing directory");
    }
    i++;
  }
  return path;
}
