#include "maki/core/seq/seq_io.hpp"

std::string toString(const Dna4Sequence &seq)
{
  std::string result;
  result.reserve(seq.size());

  for (auto c : seq)
    result.push_back(seqan3::to_char(c));

  return result;
}

std::string toString(std::vector<Dna4Sequence> const &v, const char *delim)
{
  size_t resultSize = 0;
  for (auto const &seq : v)
    resultSize += seq.size();
  resultSize += strlen(delim) * (v.size() - 1);

  std::string result;
  result.reserve(resultSize);

  for (auto const &seq : v)
  {
    if (!result.empty())
      result.append(delim);
    for (auto c : seq)
      result.push_back(seqan3::to_char(c));
  }

  return result;
}

std::string_view removeCompressedExtensions(std::string_view filePath)
    {
        for (auto ext : COMPRESSED_EXTENSIONS)
        {
            if (filePath.ends_with(ext))
            {
                return filePath.substr(0, filePath.size() - ext.size());
            }
        }

        return filePath;
    }

    InputFileType detectFileType(std::string_view inputFile)
    {
        auto inputFileFmt = removeCompressedExtensions(inputFile);

        for (auto fnaExt : FNA_EXTENSIONS)
        {
            if (inputFileFmt.ends_with(fnaExt))
            {
                return FastaFileType;
            }
        }

        for (auto gffExt : GFF_EXTENSIONS)
        {
            if (inputFileFmt.ends_with(gffExt))
            {
                return Gff3FileType;
            }
        }

        for (auto fqExt : FQ_EXTENSIONS)
        {
            if (inputFileFmt.ends_with(fqExt))
            {
                return FastQFileType;
            }
        }

        return UnknownFileType;
    }

    std::vector<fs::path> readFilePaths(const char *manifest_file)
  {
    std::ifstream manifest(manifest_file, InputFileType filter);
    // read lines from file
    std::vector<fs::path> files;
    std::copy(std::istream_iterator<fs::path>(manifest),
              std::istream_iterator<fs::path>(),
              std::back_inserter(files));
    // validate all these files exist
    std::size_t missing = 0;
    for (const fs::path &file : files)
    {
      if (!fs::exists(file))
      {
        ++missing;
        LOG(ERROR) << file << " not found";
      }
      InputFileType ftype = detectFileType(file);
      if (ftype != filter)
      {
        ++missing;
        LOG(ERROR) << "invalid file type: " << file;
      }
    }
    if (missing)
    {
      LOG(ERROR) << missing << " files not found";
      throw "errors during parsing";
    }
    else
    {
      LOG(INFO) << "parsed " << files.size() << " files";
    }
    return files;
  }

    std::string_view extractSequenceName(std::string_view path)
    {
        size_t newStart = path.find_last_of('/', (size_t)-1);
        if (newStart != (size_t)-1)
        {
            path = path.substr(newStart + 1);
        }
        path = removeCompressedExtensions(path);

        for (auto fnaExt : FNA_EXTENSIONS)
        {
            if (path.ends_with(fnaExt))
            {
                return path.substr(0, path.size() - fnaExt.size());
            }
        }

        for (auto gffExt : GFF_EXTENSIONS)
        {
            if (path.ends_with(gffExt))
            {
                return path.substr(0, path.size() - gffExt.size());
                ;
            }
        }

        for (auto fqExt : FQ_EXTENSIONS)
        {
            if (path.ends_with(fqExt))
            {
                return path.substr(0, path.size() - fqExt.size());
            }
        }

        return path;
    }
