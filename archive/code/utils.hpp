#pragma once
#include <array>
#include <bit>
#include <bitset>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <random>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>
#include <glog/logging.h>
#include <oneapi/tbb/concurrent_queue.h>
#include <cereal/archives/binary.hpp>
#include <MurmurHash/MurmurHash3.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/sfstream.hpp>

namespace fs = std::filesystem;

/**
 * Memory alignment
 */

// https://stackoverflow.com/a/466242
inline constexpr uint32_t nextPowerOf2(uint32_t n) {
    --n;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    ++n;
    return n;
}

// https://stackoverflow.com/a/466242
inline constexpr uint64_t nextPowerOf2(uint64_t n) {
    --n;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n |= n >> 32;
    ++n;
    return n;
}

template <typename T>
inline constexpr uint8_t floor_log2(T __x) {
    if (__x == 0) return 0;
    return uint8_t(sizeof(T) * 8 - 1)
        - static_cast<uint8_t>(std::countl_zero(__x));
}

// number of bits required to represent a number
template <typename T>
inline constexpr uint8_t required_bits(T __x) {
    return floor_log2(__x) + 1u; 
}

template <typename T>
inline constexpr uint8_t ceil_log2(T __x) {
    if (__x == 0) return 0;
    return floor_log2(__x - 1u) + 1u;
}

// number of bytes required to represent a number
template <typename T>
inline constexpr uint8_t required_bytes(T __x) {
    return (ceil_log2(__x) + 7) / 8; 
}

template <typename T>
inline constexpr long double floor_log2_l(T __x) {
    return static_cast<long double>(floor_log2(__x));
}

template <typename T>
inline constexpr long double ceil_log2_l(T __x) {
    return static_cast<long double>(ceil_log2(__x));
}

namespace _mkimem_details
{

    // number of bytes per key
    template <uint8_t bitsPerNucl = 2>
    inline constexpr uint8_t key_size(uint8_t k) {
        return ((bitsPerNucl * k) + 7) / 8;
    }

    // number of bytes per value
    inline constexpr uint8_t value_size(uint64_t colourBits) {
        return (colourBits + /* 3 bits per edge */ 3u + 7u) / 8u;
    }

    // number of bytes per record
    inline constexpr uint32_t record_size(uint8_t keyBytes, uint8_t valueBytes) {
        return static_cast<uint32_t>(keyBytes + valueBytes);
        // return 8u * static_cast<uint32_t>(
        //     std::ceil(static_cast<float>(keyBytes + valueBytes) / 8.0f));
    }

    // estimate number of repeated k-mers for a given sequence length under uniform assumptions
    inline constexpr size_t prob_repeat(size_t seq_len_, uint8_t k_) {
        auto seq_len_dbl__ = static_cast<double>(seq_len_),
            rho_dbl__ = seq_len_dbl__ / std::pow(4.0, static_cast<double>(k_));
        return static_cast<size_t>(std::ceil(
            (seq_len_dbl__ * rho_dbl__) / (1.0 + rho_dbl__)
        ));
    }

} // namespace _mkimem_details

// copyable atomic wrapper --------------------------------------------------------

template <typename T>
struct copyable_atomic
{
    constexpr copyable_atomic() = default;

    constexpr copyable_atomic(T desired) noexcept : value(desired) {}

    constexpr copyable_atomic(const copyable_atomic &rhs) = default;

    constexpr copyable_atomic(copyable_atomic &&rhs) = default;

    inline constexpr copyable_atomic &operator=(const copyable_atomic &) = default;

    inline constexpr copyable_atomic &operator=(copyable_atomic &&) = default;

    explicit operator T() const noexcept
    { return value; }

    explicit operator bool() const noexcept
    { return value; }

    // atomic increment
    constexpr inline auto operator++() noexcept
    { return ++std::atomic_ref<T>(value); }

    inline T operator++(int) noexcept
    { return std::atomic_ref<T>(value)++; }

    inline auto operator+=(T dv) noexcept
    { return std::atomic_ref<T>(value) += dv; }

    inline void store(T newValue)
    { std::atomic_ref<T>(value).store(newValue); }

    inline T load() const
    { return std::atomic_ref<const T>(value).load(); }

    inline std::atomic_ref<T> ref()
    { return std::atomic_ref<T>(value); }

    inline std::atomic_ref<const T> ref() const
    { return std::atomic_ref<const T>(value); }

    const inline auto operator|=(T x_) noexcept
    { auto _v = std::atomic_ref<T>(value); return _v |= x_; }

protected:
    T value;
};

typedef copyable_atomic<uint64_t> copyable_atomic_uint64_t;
typedef copyable_atomic<uint32_t> copyable_atomic_uint32_t;

// copyable pair ---------------------------------------------------------------------

template <typename T>
struct copyable_pair
{
    copyable_pair() : first(), second() {}
    copyable_pair(const copyable_pair&) =default;
    copyable_pair(copyable_pair&&) =default;
    copyable_pair(const T &lhs, const T &rhs) : first(lhs), second(rhs) {}
    copyable_pair(T&& lhs, T&& rhs) : first(std::move(lhs)), second(std::move(rhs)) {}

    inline copyable_pair& operator=(const copyable_pair&) =default;
    inline copyable_pair& operator=(copyable_pair&&) =default;

    inline bool operator==(copyable_pair const &rhs) const =default;
    inline bool operator!=(copyable_pair const &rhs) const =default;

    T first, second;
};

// 32-bit hash
template <typename T>
struct copyable_pair_hash_32
{
    inline constexpr uint32_t operator()(T lhs, T rhs) const
    {
        constexpr size_t _W = sizeof(T);
        uint8_t _data[2 * _W];
        std::memcpy(_data, &lhs, _W);
        std::memcpy(_data + _W, &rhs, _W);

        uint32_t _out;
        MurmurHash3_x86_32(_data, 2 * _W, 1u, &_out);
        return _out;
    }

    // 32-bit hash
    inline constexpr uint32_t operator()(copyable_pair<T> const &v) const
    { return operator()(v.first, v.second); }
};

// 64-bit hash
template <typename T>
struct copyable_pair_hash_64
{
    constexpr uint64_t hash(copyable_pair<T> const &v) const
    {
        constexpr size_t _W = sizeof(T);
        uint8_t _data[2 * _W], _out[16];
        std::memcpy(_data, &v.first, _W);
        std::memcpy(_data + _W, &v.second, _W);
        MurmurHash3_x86_128(_data, 2 * _W, 1u, _out);

        uint64_t res;
        std::memcpy(&res, _out, 8);

        return res;
    }

    inline constexpr uint64_t operator()(copyable_pair<T> const &v) const
    { return hash(v); }

    inline constexpr bool equal(copyable_pair<T> const &lhs, copyable_pair<T> const &rhs) const
    { return lhs == rhs; }
};

template <typename Container>
inline auto __valarray_as_vec(Container _arr)
{
    return std::vector<size_t>{_arr[0], _arr[1], _arr[2], _arr[3], _arr[4]};
}

/**
 * Filesystem
 */

namespace tempio
{

    // taken from https://stackoverflow.com/a/58454949
    fs::path create_temporary_directory(unsigned long long max_tries = 1000);

    fs::path create_temporary_file(fs::path const &base, const char *ext, unsigned long long max_tries = 1000);

} // namespace tempio

namespace verbose
{

    template <std::input_iterator It, std::sentinel_for<It> S>
    std::string to_byte_string(It begin_, S end_)
    {
        std::string result;

        if (begin_ == end_)
        {
            return result;
        }

        result.append("0b" + std::bitset<8>(*begin_).to_string());
        ++begin_;

        while (begin_ != end_)
        {
            result.append(", 0b" + std::bitset<8>(*begin_).to_string());
            ++begin_;
        }

        return result;
    }

    std::string mapToString(const std::vector<std::pair<std::string, std::string>> &);

} // namespace verbose

/**
 * Utilities
 */

namespace utility
{

    template <bool B>
    struct in_place_bool_t
    {
        explicit in_place_bool_t() = default;
    };

    template <bool B>
    inline constexpr in_place_bool_t<B> in_place_bool{};

} // namespace utility

namespace detail
{

    size_t countTsvWords(const std::string &str);

} // namespace detail

namespace fileutils
{

    // Formatting ---------------------------

    enum InputFileType
    {
        FastaFileType,
        Gff3FileType,
        FastQFileType,
        UnknownFileType
    };

    constexpr std::string_view COMPRESSED_EXTENSIONS[2] = {".tar.gz", ".gz"};
    constexpr std::string_view FNA_EXTENSIONS[3] = {".fna", ".fa", ".fasta"};
    constexpr std::string_view GFF_EXTENSIONS[2] = {".gff3", ".gff"};
    constexpr std::string_view FQ_EXTENSIONS[2] = {".fq", ".fastq"};

    std::string_view removeCompressedExtensions(std::string_view filePath);

    InputFileType detectFileType(const std::string &inputFile);

    InputFileType detectFileType(std::string_view inputFile);

    std::string_view extractSequenceName(std::string_view path);

    // IO ----

    inline fs::path as_path(const char *output_path) { return fs::path(output_path); }

    std::vector<fs::path> list_fasta_files(fs::path base, std::string const &ext);

    template <typename InputIt>
    inline void write_lines(fs::path file, InputIt it, InputIt end)
    {
        std::ofstream fh(file.c_str());
        for (; it != end; ++it)
        {
            fh << it->string() << std::endl;
        }
        fh.close();
    }

    std::vector<std::string> read_lines(fs::path file);

    std::vector<std::string> read_lines(fs::path file, std::function<std::string_view(std::string_view)> _formatter);

    // Serialisation --------------------------------

    template <uint8_t t_width>
    bool store_to_file(sdsl::int_vector<t_width> const &v, std::string const &file)
    {
        sdsl::osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);

        if (!out)
        {
            LOG(ERROR) << "fileutils::store_to_file:: Could not open file `" << file << "`";
            return false;
        }

        cereal::BinaryOutputArchive archive(out);
        v.save(archive);

        out.close();
        return true;
    }

    template <class T, size_t N>
    bool store_to_file(std::array<T, N> const &arr, std::string const &file)
    {
        sdsl::osfstream out(file, std::ios::binary | std::ios::trunc | std::ios::out);

        if (!out)
        {
            LOG(ERROR) << "fileutils::store_to_file:: Could not open file `" << file << "`";
            return false;
        }

        cereal::BinaryOutputArchive archive(out);
        archive(arr);

        out.close();
        return true;
    }

    // De-serialisation --------------------------------

    template <uint8_t t_width>
    bool load_from_file(sdsl::int_vector<t_width> &v, std::string const &file)
    {
        sdsl::isfstream in(file, std::ios::binary | std::ios::in);

        if (!in)
        {
            LOG(ERROR) << "Could not load file `" << file << "`";
            return false;
        }

        cereal::BinaryInputArchive iarchive(in);
        v.load(iarchive);

        in.close();
        return true;
    }

    template <class T, size_t N>
    bool load_from_file(std::array<T, N> &arr, std::string const &file)
    {
        sdsl::isfstream in(file, std::ios::binary | std::ios::in);

        if (!in)
        {
            LOG(ERROR) << "Could not load file `" << file << "`";
            return false;
        }

        cereal::BinaryInputArchive iarchive(in);
        iarchive(arr);

        in.close();
        return true;
    }

} // namespace fileutils

namespace stream {

    template <typename T>
    class OutputHandler {
    protected:
        std::ostream& os;

    public:
        OutputHandler(std::ostream &os_) : os(os_) {}
        virtual ~OutputHandler() = default;
        virtual void write(T& data) =0;
    };

    template <typename T>
    struct DirectOutput : public OutputHandler<T> {
    };

    class StreamManager {
        std::map<std::string, std::unique_ptr<std::ostream>> streams;
    public:
        std::ostream& operator[](const std::string& filename);
    };

} // namespace stream
