

## Writing

```cpp
#include "archive_writer.hpp"
#include "raw_packer.hpp"
#include <vector>
#include <iostream>

int main() {
    using namespace charc;

    ArchiveWriter::Options opt;
    opt.preallocate_bytes = 512ull * 1024 * 1024; // optional
    opt.chunk_alignment = 4096;                   // align chunks to 4KiB

    ArchiveWriter writer("dataset.char", opt);

    // Example: write RAW_PACKED chunk
    std::vector<uint64_t> values(10'000'000);
    for (uint64_t i = 0; i < values.size(); ++i) values[i] = i & 0xFFFFu;
    std::vector<uint8_t> packed;
    pack_raw_values(values.data(), values.size(), /*bit_width=*/16, packed);
    writer.write_raw_packed(packed.data(), packed.size(), values.size(), /*bit_width=*/16);

    // Example: write SDSL chunk (assuming you have serialized bytes)
    // std::vector<uint8_t> sdsl_bytes = ...;
    // writer.write_sdsl_serialized(sdsl_bytes, elem_count, bit_width);

    writer.finalize();
    std::cout << "Wrote dataset.char\n";
}
```

## Reading

```cpp
#include "archive_reader.hpp"
#include "array_factory.hpp"
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) { std::cerr << "Usage: " << argv[0] << " dataset.char\n"; return 1; }

    try {
        charc::ArchiveReader reader(argv[1]);
        DefaultArrayFactory factory;

        std::cout << "Chunks: " << reader.chunk_count()
                  << ", total_elems: " << reader.total_elems() << "\n";

        // Random point lookup by global index
        if (reader.total_elems() > 0) {
            uint64_t idx = reader.total_elems() / 2;
            auto k = reader.find_chunk_by_global_index(idx);
            if (k >= 0) {
                const auto& meta = reader.meta(static_cast<size_t>(k));
                auto [ptr, len] = reader.chunk_region(static_cast<size_t>(k));

                // Optimize for scanning within this chunk:
                reader.advise_chunk_willneed(static_cast<size_t>(k)); // or advise_chunk_sequential(k)

                auto arr = factory.open_from_region(ptr, len, meta);
                const uint64_t local = idx - meta.start_index;
                std::cout << "idx=" << idx
                          << " in chunk " << k
                          << " local=" << local
                          << " value=" << arr->get(local) << "\n";
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
```

