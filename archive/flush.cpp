#include <maki/flush.hpp>

size_t flush::ranges::detail::calculateLeftChunks(size_t width, size_t chunksize, oneapi::tbb::proportional_split p) {
    size_t sizeEst = width * p.left() / (p.left() + p.right());
    return (sizeEst + chunksize - 1) / chunksize;
}

