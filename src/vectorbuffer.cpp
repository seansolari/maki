#include <maki/vectorbuffer.hpp>

ZstdCompressor::ZstdCompressor()
    : bytes(),
      cctx(ZSTD_createCCtx(), ZstdDeleter()),
      rsl()
{
    if (cctx.get() == NULL)
    {
        LOG(ERROR) << "ZSTD_createCCtx() failed!";
        throw std::bad_alloc();
    }
}

void checkZstd(const size_t err, const char *fail_msg)
{
    if (ZSTD_isError(err))
    {
        throw std::runtime_error(std::string("Zstd check failed: ") + std::string(ZSTD_getErrorName(err)) + std::string(fail_msg));
    }
}
