#pragma once
#include <cstdint>
#include <fstream>
#include <set>
#include <string>
#include <utility>
#include <glog/logging.h>
#include <oneapi/tbb.h>
#include <maki/utils.hpp>

namespace colour {
    namespace encoding {
        template <typename T, typename U>
        struct FeatureIdAssigner {
            FeatureIdAssigner(uint64_t bitWidth, uint64_t baseId, T &nameBuffer_, std::vector<U> &featureBuffer_);
        protected:
            uint64_t _w, _base, _cid;
            T &nameBuffer;
            std::vector<U> &featureBuffer;
        public:
            uint64_t baseId() const { return _base; }
            uint64_t nextId() { return _base | (_cid++ << _w); }
            void assignName(const T &inputName) { nameBuffer = inputName; }
            void assignFeatures(std::vector<U> &&buffer) { featureBuffer = std::move(buffer); }
        };

        template <typename T, typename U>
        struct FeatureIdGenerator {
            FeatureIdGenerator(size_t numElements);
        protected:
            uint64_t _w;
            std::vector<T> featureNames;
            std::vector<std::vector<U>> features;
        public:
            FeatureIdAssigner<T, U> specify(uint64_t baseId) { return FeatureIdAssigner<T, U>(_w, baseId, featureNames[baseId], features[baseId]); }
            FeatureIdAssigner<T, U> specify(uint64_t baseId, size_t featureIndex) { return FeatureIdAssigner<T, U>(_w, baseId, featureNames[featureIndex], features[featureIndex]); }
            uint64_t getBaseFeatureWidth() const noexcept { return _w; }
            void writeFeatureNamesTo(const std::string &file) const;
            void clearFeatureNames() { std::vector<T>().swap(featureNames); }
            void writeFeaturesTo(const std::string &file) const;
            void clearFeatures() { std::vector<std::vector<U>>().swap(features); }
        };

    } // namespace encoding

} // namespace colour

class BufferValue {
    static constexpr uint64_t _msk = 0xFFu;
public:
    BufferValue() : _data(0ULL) {}
    BufferValue(uint64_t _edge) : _data(_edge) {}
    BufferValue(uint64_t _colour, uint64_t _edge) : _data((_colour << 3) | _edge) {}
    // read bytes from least- to most-significant
    explicit BufferValue(const uint8_t *src_, size_t n_);
public:
    bool operator==(const BufferValue &other) const noexcept { return _data == other._data; }
    bool operator!=(const BufferValue &other) const noexcept { return _data != other._data; }
    operator uint64_t() const { return _data; }
    uint64_t colour() const { return _data >> 3; }
    uint64_t edge() const { return _data & 0b111ULL; }
    uint64_t *data() { return &_data; }
public:
    // serialise bytes from least- to most-significant
    void flush(uint8_t *dest_, size_t n_) const;
private:
    uint64_t _data;
};

// Template definitions

template<class T, class U>
colour::encoding::FeatureIdAssigner<T, U>::FeatureIdAssigner(uint64_t bitWidth, uint64_t baseId, T &nameBuffer_, std::vector<U> &featureBuffer_)
    : _w(bitWidth)
    , _base(baseId)
    , _cid(1)
    , nameBuffer(nameBuffer_)
    , featureBuffer(featureBuffer_)
{}

template<class T, class U>
colour::encoding::FeatureIdGenerator<T, U>::FeatureIdGenerator(size_t numElements)
    : _w(ceil_log2(numElements + 1 /* null item */))
    , featureNames(numElements)
    , features(numElements)
{}

template<class T, class U>
void colour::encoding::FeatureIdGenerator<T, U>::writeFeatureNamesTo(const std::string &file) const {
    stream::StreamManager mgr;
    std::ostream &os = mgr[file];
    T::writeHeader(os);
    for (size_t i = 0; i < featureNames.size(); ++i)
        os << i + 1 << '\t'<< featureNames[i]<< '\n';
}

template<class T, class U>
void colour::encoding::FeatureIdGenerator<T, U>::writeFeaturesTo(const std::string &file) const {
    stream::StreamManager mgr;
    std::ostream &os = mgr[file];
    U::writeHeader(os);
    for (const std::vector<U> &vec : features)
        for (const U &val : vec)
            os << val << '\n';
}
