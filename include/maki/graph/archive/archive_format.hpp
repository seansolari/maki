
#pragma once
#include <cstdint>

enum class Encoding : uint8_t
{
  RAW_PACKED = 1
};

#pragma pack(push, 1)
struct ChunkMeta
{
  uint64_t start_index; // prefix sum
  uint64_t elem_count;  // number of elements in chunk
  uint64_t file_offset; // start byte in file
  uint64_t byte_len;    // packed region length
  uint8_t bit_width;    // 1..64
  uint8_t encoding;     // RAW_PACKED
  uint16_t alignment;   // chunk alignment (0 = none)
  uint32_t crc32;       // per-chunk checksum over packed bytes
};

struct Footer
{
  uint64_t toc_offset;
  uint64_t chunk_count;
  uint32_t version; // = 1
  uint32_t magic;   // "CHAR"
};
#pragma pack(pop)

static constexpr uint32_t MAGIC_CHAR = 0x52414843u; // 'C','H','A','R'
static constexpr uint32_t VERSION = 1;
