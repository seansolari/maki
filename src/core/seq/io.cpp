#include "maki/core/seq/io.hpp"
#include "maki/core/utils/logging.hpp"

#include <filesystem>

std::string toString(const Dna4Sequence &seq) {
  std::string result;
  result.reserve(seq.size());

  for (auto c : seq)
    result.push_back(seqan3::to_char(c));

  return result;
}

std::string toString(std::vector<Dna4Sequence> const &v, const char *delim) {
  size_t resultSize = 0;
  for (auto const &seq : v)
    resultSize += seq.size();
  resultSize += strlen(delim) * (v.size() - 1);

  std::string result;
  result.reserve(resultSize);

  for (auto const &seq : v) {
    if (!result.empty())
      result.append(delim);
    for (auto c : seq)
      result.push_back(seqan3::to_char(c));
  }

  return result;
}

std::string_view removeCompressedExtensions(std::string_view filePath) {
  for (auto ext : COMPRESSED_EXTENSIONS) {
    if (filePath.ends_with(ext)) {
      return filePath.substr(0, filePath.size() - ext.size());
    }
  }

  return filePath;
}

InputFileType detectFileType(std::string_view inputFile) {
  auto inputFileFmt = removeCompressedExtensions(inputFile);

  for (auto fnaExt : FNA_EXTENSIONS) {
    if (inputFileFmt.ends_with(fnaExt)) {
      return FastaFileType;
    }
  }

  for (auto gffExt : GFF_EXTENSIONS) {
    if (inputFileFmt.ends_with(gffExt)) {
      return Gff3FileType;
    }
  }

  for (auto fqExt : FQ_EXTENSIONS) {
    if (inputFileFmt.ends_with(fqExt)) {
      return FastQFileType;
    }
  }

  return UnknownFileType;
}

namespace detail {

std::size_t sepFind(const std::string &s, std::size_t t, const char sep) {
  std::size_t pos = 0;
  for (std::size_t occs = 0; pos < s.size() && occs < t; ++pos)
    if (s[pos] == sep)
      ++occs;
  if (pos == s.size())
    throw std::runtime_error(std::string{"Index out of range: "} +
                             std::to_string(t));
  return pos;
}

std::pair<std::size_t, std::size_t> sepSelect(std::string const &s,
                                              std::size_t t, const char sep) {
  std::size_t begin = sepFind(s, t, sep), end = s.find(sep, begin + 1);
  if (end == std::string::npos)
    end = s.size();
  return std::make_pair(begin, end);
}

} // namespace detail

GenomeManifest readFilePaths(const char *manifest_file, std::size_t col,
                             const char sep, InputFileType filter) {
  GenomeManifest manifest{.files = {}, .type = filter};
  std::ifstream ifh(manifest_file);

  // read lines from file
  std::string line;
  if (ifh.is_open()) {
    while (std::getline(ifh, line)) {
      auto [l, r] = detail::sepSelect(line, col, sep);
      manifest.files.emplace_back(line.substr(l, r - l));
    }
  } else {
    LOG_ERROR() << "Unable to open file" << manifest_file;
  }

  // validate all these files exist
  std::size_t missing = 0;
  for (const std::string &file : manifest.files) {
    if (!std::filesystem::exists(file)) {
      ++missing;
      LOG_ERROR() << file << " not found";
    }
    InputFileType ftype = detectFileType(file);

    if (ftype != filter) {
      ++missing;
      LOG_ERROR() << "invalid file type: " << file;
    }
  }
  if (missing) {
    LOG_ERROR() << missing << " files not found";
    throw "errors during parsing";
  } else {
    LOG_INFO() << "parsed " << manifest.files.size() << " files";
  }
  return manifest;
}

std::string_view extractSequenceName(std::string_view path) {
  size_t newStart = path.find_last_of('/', (size_t)-1);
  if (newStart != (size_t)-1) {
    path = path.substr(newStart + 1);
  }
  path = removeCompressedExtensions(path);

  for (auto fnaExt : FNA_EXTENSIONS) {
    if (path.ends_with(fnaExt)) {
      return path.substr(0, path.size() - fnaExt.size());
    }
  }

  for (auto gffExt : GFF_EXTENSIONS) {
    if (path.ends_with(gffExt)) {
      return path.substr(0, path.size() - gffExt.size());
      ;
    }
  }

  for (auto fqExt : FQ_EXTENSIONS) {
    if (path.ends_with(fqExt)) {
      return path.substr(0, path.size() - fqExt.size());
    }
  }

  return path;
}

std::ostream &operator<<(std::ostream &os, const DataFilePair &obj) {
  return os << "PairedReads(" << obj.forwardFile << ", " << obj.reverseFile
            << ")";
}

bool DataFilePair::gzCompressed() const {
  bool lCx = forwardFile.ends_with(".gz"), rCx = reverseFile.ends_with(".gz");

  if (lCx != rCx)
    throw std::runtime_error("paired files " + forwardFile + " and " +
                             reverseFile + " have different compression");
  else
    return lCx;
}
