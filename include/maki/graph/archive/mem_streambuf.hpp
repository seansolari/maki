
#pragma once
#include <streambuf>
#include <cstdint>

class MemStreamBuf : public std::streambuf {
public:
    MemStreamBuf(const uint8_t* base, size_t size) {
        char* p = const_cast<char*>(reinterpret_cast<const char*>(base));
        setg(p, p, p + size);
    }
};
