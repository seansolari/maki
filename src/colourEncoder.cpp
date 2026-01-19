#include <maki/colourEncoder.hpp>

BufferValue::BufferValue(const uint8_t *src_, size_t n_) : _data(0) {
    for (size_t i = 0; i < n_; ++i) {
        uint64_t val = *(src_ + i);
        _data |= val << (8 * i);
    }
}

void BufferValue::flush(uint8_t *dest_, size_t n_) const {
    for (size_t i = 0; i < n_; ++i) {
        uint8_t v = static_cast<uint8_t>((_data >> (8 * i)) & _msk);
        *(dest_ + i) = v;
    }
}
