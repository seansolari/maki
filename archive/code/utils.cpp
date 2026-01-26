#include <maki/utils.hpp>

#include <algorithm>
#include <zstr.hpp>

namespace tempio
{

    fs::path create_temporary_directory(unsigned long long max_tries)
    {
        auto tmp_dir = fs::temp_directory_path();
        unsigned long long i = 0;
        std::random_device dev;
        std::mt19937 prng(dev());
        std::uniform_int_distribution<uint64_t> rand(0);
        fs::path path;
        while (true)
        {
            std::stringstream ss;
            ss << std::hex << rand(prng);
            path = tmp_dir / ss.str();
            // true if the directory was created.
            if (fs::create_directory(path))
            {
                break;
            }
            if (i == max_tries)
            {
                throw std::runtime_error("could not find non-existing directory");
            }
            i++;
        }
        return path;
    }

    fs::path create_temporary_file(fs::path const &base, const char *ext, unsigned long long max_tries)
    {
        unsigned long long i = 0;
        std::random_device dev;
        std::mt19937 prng(dev());
        std::uniform_int_distribution<uint64_t> rand(0);
        fs::path path;
        while (true)
        {
            std::stringstream ss;
            ss << std::hex << rand(prng) << ext;
            path = base / ss.str();
            // true if the directory was created.
            if (!fs::exists(path))
            {
                break;
            }
            if (i == max_tries)
            {
                throw std::runtime_error("could not find non-existing directory");
            }
            i++;
        }
        return path;
    }

} // namespace tempio

size_t detail::countTsvWords(const std::string &str)
{
    return std::ranges::count(str, '\t') + 1;
}

std::string verbose::mapToString(const std::vector<std::pair<std::string, std::string>> &values)
{
    return std::transform_reduce(
        values.cbegin(),
        values.cend(),
        std::string(),
        [](const std::string &lhs, const std::string &rhs) -> std::string
        { return lhs.size() == 0 ? rhs : lhs + ";" + rhs; },
        [](const std::pair<std::string, std::string> &value) -> std::string
        { return value.first + "=" + value.second; });
}

namespace fileutils
{

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

    InputFileType detectFileType(const std::string &inputFile)
    {
        return detectFileType(std::string_view(inputFile));
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

    std::vector<fs::path> list_fasta_files(fs::path base, std::string const &ext = ".fna")
    {
        std::vector<fs::path> paths;
        for (const auto &entry : fs::directory_iterator(base))
        {
            if (fs::is_regular_file(entry) && entry.path().extension() == ext)
                paths.emplace_back(entry.path());
        }
        return paths;
    }

    std::vector<std::string> read_lines(fs::path file)
    {
        std::vector<std::string> lines;
        std::ifstream fh;
        fh.open(file.c_str());
        std::string line;
        while (std::getline(fh, line))
        {
            lines.push_back(line);
        }
        return lines;
    }

    std::vector<std::string> read_lines(fs::path file, std::function<std::string_view(std::string_view)> _formatter)
    {
        std::vector<std::string> lines;
        std::ifstream fh;
        fh.open(file.c_str());
        std::string line;
        while (std::getline(fh, line))
        {
            lines.emplace_back(_formatter(line));
        }
        return lines;
    }

} // namespace fileutils

std::ostream &stream::StreamManager::operator[](const std::string &filename) {
    if (filename.empty()) {
        return std::cout;
    } else {
        if (streams.find(filename) == streams.end()) {
            if (filename.ends_with(".gz")) {
                streams[filename] = std::make_unique<zstr::ofstream>(filename);
            } else {
                streams[filename] = std::make_unique<std::ofstream>(filename);
            }
        }
        return *streams[filename];
    }
}
