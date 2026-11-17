
#pragma once
#include <filesystem>

namespace fs = std::filesystem;

namespace tempio {

fs::path create_temporary_directory(fs::path tmp_dir,
                                    unsigned long long max_tries = 1000);
fs::path create_temporary_directory(unsigned long long max_tries = 1000);
fs::path create_temporary_file(fs::path const &base, const char *ext,
                               unsigned long long max_tries = 1000);

class TempFile {
public:
  explicit TempFile(
      std::string const &contents,
      std::filesystem::path tmpdir = std::filesystem::temp_directory_path(),
      const char *ext = ".fna");

  ~TempFile();

  constexpr std::filesystem::path const &path() const { return path_; }

private:
  std::filesystem::path path_;
};

class TempFastqFile : public TempFile {
public:
  constexpr explicit TempFastqFile(
      std::string const &contents,
      std::filesystem::path tmpdir = std::filesystem::temp_directory_path())
      : TempFile(contents, tmpdir, ".fq") {}
};

class TempFastaFile : public TempFile {
public:
  constexpr explicit TempFastaFile(
      std::string const &contents,
      std::filesystem::path tmpdir = std::filesystem::temp_directory_path())
      : TempFile(contents, tmpdir, ".fna") {}
};

} // namespace tempio
