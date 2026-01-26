
#pragma once
#include <filesystem>

namespace fs = std::filesystem;

fs::path create_temporary_directory(fs::path tmp_dir,
                                    unsigned long long max_tries = 1000);
fs::path create_temporary_directory(unsigned long long max_tries = 1000);
fs::path create_temporary_file(fs::path const &base, const char *ext,
                               unsigned long long max_tries = 1000);