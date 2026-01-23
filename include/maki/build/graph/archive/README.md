# Archiving design

## **Payload–Sink Pattern: Design & Implementation Guide**

### 1) Aim

Many applications need to **produce chunks in parallel** and **write them to multiple outputs**—e.g., a giant byte-array file, a chunked archive with CRCs, and an SDSL-compatible vector—**all in the same global order**. The challenge is to do this:

*   **Safely** (ordered writes, fault containment),
*   **Fast** (sequential I/O, mmappable files, no per-call virtual dispatch),
*   **Memory-bounded** (predictable footprint via pooling),
*   **Composable** (easy to add more sinks later without refactoring the pipeline).

We achieve this with the **Payload–Sink** pattern:

*   **Payload** = the minimal, type-safe data structure a sink needs for one chunk.
*   **Sink** = a small concrete type with a uniform interface: `write(const Payload&)`, `finalize()`.
*   **Bundle** = a single object carrying the **tuple of per-sink payloads** for a given chunk id.
*   **MultiSink** = a compile-time composite that fans a bundle’s payloads out to all sinks.
*   **Ordered Writer (Collector)** = a single-threaded in-order consumer enforcing chunk-id order and calling `MultiSink::write_bundle(...)`.
*   **Bundle Pool** = a bounded pool of reusable `Bundle`s enabling backpressure.

> Think of a restaurant kitchen: cooks (producers) assemble a **single tray** (bundle) that has the dish for each table (sink). One expediter (writer) sends trays to the dining room **in ticket order**. The tray design ensures every table gets the right plate, and no single cook needs to know where tables sit.

***

### 2) Core Abstractions

#### 2.1 Sinks as concrete types

Each sink is a concrete class with a uniform API—no base classes, no RTTI:

```cpp
struct MySink {
  using Payload = /* per-sink payload type */;
  void write(const Payload& p);  // called in strictly increasing chunk id
  void finalize();               // called once at end
};
```

#### 2.2 Map sink → payload at compile time

Define a trait to pair each sink type with the payload type it expects:

```cpp
template<class Sink> struct sink_payload { using type = /* default */; };
// Specialize for sinks:
template<> struct sink_payload<MySink> { using type = MyPayload; };
```

#### 2.3 Bundle: one submission per chunk

At **compile time**, we know the sink set. We form a bundle that carries **all** payloads for a given chunk:

```cpp
template<class... Sinks>
struct ChunkBundleT {
  uint64_t id;
  std::tuple<typename sink_payload<Sinks>::type...> payloads;
};
```

Producers **fill a single bundle and submit it once**. The writer/collector fans it out.

#### 2.4 MultiSink: composite fan-out

A composite that writes a bundle to each sink in order:

```cpp
template<class... Sinks>
class MultiSink {
public:
  explicit MultiSink(std::tuple<Sinks...> sinks);
  void write_bundle(const ChunkBundleT<Sinks...>& b); // for each i: sinks[i].write(b.payloads[i])
  void finalize();                                    // finalize each sink
};
```

#### 2.5 Ordered Writer (Collector)

Single-threaded actor that:

*   Pulls bundles from a blocking queue,
*   Buffers out-of-order arrivals in a map,
*   Emits in-order `[next_id, next_id+1, ...]`,
*   Calls `MultiSink::write_bundle`,
*   Returns bundles to the pool,
*   Calls `finalize()` at the end.

Why single-threaded? It maximizes sequential I/O and avoids inter-sink coordination complexity—heavy CPU work stays in producers.

#### 2.6 Bundle Pool (backpressure)

A bounded pool of pre-allocated bundles:

*   Controls memory usage,
*   Enables producer backpressure when sinks slow down,
*   Preserves capacity of internal vectors to minimize reallocation.

***

### 3) Policies (cross-cutting behaviors)

We use simple policy classes to control OS-level and file-level behaviors **without runtime branching**:

*   `AlignTo<N>` / `NoAlignment` — pad chunk starts (appropriate for **structured/chunked** formats; **not** for pure concatenations).
*   `Preallocate<Bytes>` / `NoPreallocate` — reduce late allocation stalls.
*   `FadviseSequential` / `NoFadvise` — give kernel hints for readahead/writeback.
*   `CRC32` / `NoChecksum` (for archives) — integrity at chunk granularity.

Policies are **template parameters** to sinks (or writers), e.g.:

```cpp
using ArchiveSinkA = ArchiveFileSink<AlignTo<4096>, Preallocate<2_GiB>, FadviseSequential, CRC32>;
```

***

### 4) Producer Interface (how to provide data to sinks)

**Producers never talk to sinks directly.** They:

1.  `auto bundle = pool.acquire();`
2.  `bundle->id = chunk_id;`
3.  Fill per-sink payloads:
    ```cpp
    auto& p0 = std::get<0>(bundle->payloads); // payload for sink0
    auto& p1 = std::get<1>(bundle->payloads); // payload for sink1
    // ...fill p0, p1...
    ```
4.  `queue.submit(std::move(bundle));`

**Payload design**: Keep payloads POD-ish and minimal (e.g., for a byte sink `std::vector<uint8_t>`, for an SDSL sink `sdsl::int_vector<>`). Include any metadata needed by the sink to write the chunk (counts, widths).

***

### 5) How to Implement a New Payload–Sink

1.  **Define the payload type** (minimal data needed per chunk):
    ```cpp
    struct MyPayload { std::vector<uint8_t> bytes; uint64_t count; };
    ```
2.  **Implement the sink class**:
    ```cpp
    template<class Align=NoAlignment, class Pre=NoPreallocate, class Adv=FadviseSequential>
    struct MySink {
      using Payload = MyPayload;
      MySink(std::string path);
      void write(const Payload& p);
      void finalize();
    };
    ```
3.  **Specialize `sink_payload<>`** so the pipeline knows what to put in the bundle:
    ```cpp
    template<class A, class P, class F>
    struct sink_payload<MySink<A,P,F>> { using type = MyPayload; };
    ```
4.  **Add the sink to the tuple** that defines the pipeline:
    ```cpp
    using Sinks = std::tuple<MySink<>, ExistingArchiveSink<>, HugeByteArraySink<> >;
    ```
5.  **Reserve per-sink payload capacity** in the pool for performance.
6.  **Test** standalone (unit-test your sink) and integrated (pipeline ordering, backpressure).

> Analogy: You’re adding a new outlet to a power strip. Define the plug shape (payload), write how the outlet handles current (sink), then plug it into the strip (tuple). The rest of the building wiring (pool, queue, writer) doesn’t change.

***

### 6) Wrapping & Runtime Control

*   **Compile-time composition** via `std::tuple<Sinks...>` gives zero-virtual overhead and inlining.
*   If you need **runtime** selection in the future, wrap sinks in `std::variant` or a simple type-erased envelope that preserves the `write/finalize` interface. The outer pipeline (pool, queue, writer, bundle) stays the same.

***

### 7) Testing Strategy

1.  **Unit tests per sink**
    *   Validate `write()` on known payloads produces exact on-disk/in-memory bytes.
    *   Validate `finalize()` side effects (TOC, CRC, sync).
    *   Error paths (short write, permission errors) propagate as exceptions.

2.  **Pipeline invariants**
    *   Out-of-order arrivals still result in ordered writes.
    *   Backpressure: small pool + many producers blocks acquire and never deadlocks.
    *   Empty payload semantics (optional sinks per chunk) are respected.

3.  **Property tests**
    *   Random payload generation → write → read back (mmapped) → verify equality.

4.  **Performance tests**
    *   Throughput under various chunk sizes (tens of MB to GB).
    *   Effect of `AlignTo<N>` and `Preallocate<T>`.
    *   Reader latencies (random point lookups, scans with `madvise`).

5.  **Robustness / Corruption**
    *   CRC validation (for archives): flip bytes, ensure detection.
    *   Bounds checks on TOC and chunk regions.

6.  **Concurrency**
    *   Many producers + single writer; repeated runs to catch heisenbugs.
    *   Parallel reader validation across threads (mmap is read-only, so lock-free).

7.  **Portability concerns**
    *   Endianness: document the on-disk endianness (we use LE).
    *   Large files: verify beyond 2^31 and 2^32 boundaries on 64-bit systems.

***

## **Chunk Archive (CHAR) Format: Design, Implementation & Usage**

### 1) Motivation & Design Choices

**Use case**: Persist huge logical arrays built from **out-of-core chunks** that arrive **roughly in order**, with:

*   **Sequential write** throughput (append-only stream),
*   **Random read** performance (mmapped),
*   **Scalable metadata** (millions of chunks feasible),
*   **Integrity**: per-chunk **CRC32**,
*   **Flexible layout**: raw-packed elements of fixed bit-width per chunk,
*   **Indexing**: O(log N) global index → chunk lookup.

**Guiding choices**:

*   **Append-only data region** → `[chunk_0][chunk_1]...[chunk_N]`.
*   A compact **TOC** at the end → `ChunkMeta[N]` contiguous array.
*   A small **footer** → `toc_offset`, `chunk_count`, `version`, `magic`.
*   **Chunk alignment** (optional) to page/hugepage boundaries for future scan speed.
*   **Preallocation** (optional) to minimize fragmentation.
*   **POSIX fadvise** for writeback hints.

This balances **write simplicity** with **read efficiency**, and separates concerns cleanly: writers enforce order, readers mmap and binary-search the TOC.

***

### 2) On-Disk Layout (Little-endian)

    [ DATA region ... variable-length chunk bytes ]
    [ TOC: ChunkMeta[chunk_count]                ]
    [ FOOTER { toc_offset, chunk_count, version, magic } ]

*   **ChunkMeta** (packed):
    *   `start_index` (global start position = prefix sum of prior chunk element counts)
    *   `elem_count`
    *   `file_offset` (byte position of this chunk in DATA region)
    *   `byte_len` (chunk size in bytes)
    *   `bit_width` (per-chunk)
    *   `encoding` (we use RAW\_PACKED = 1)
    *   `alignment` (chunk start alignment used; 0 if none)
    *   `crc32` (IEEE 802.3 polynomial over the chunk bytes)

*   **Footer**:
    *   `toc_offset`, `chunk_count`, `version` (= 1), `magic` (= "CHAR")

**Raw-packed elements**: each element uses `bit_width` bits; packed contiguously, little-endian in bytes—no padding between elements.

***

### 3) Writer (Policy-based)

**Header**: `archive_writer.hpp` (policy-based version)

*   Template policies:
    *   `AlignTo<N>` / `NoAlignment`
    *   `Preallocate<Bytes>` / `NoPreallocate`
    *   `FadviseSequential` / `NoFadvise`
    *   `CRC32` / `NoChecksum`
*   API:
    *   `append_raw_packed(data, len, elem_count, bit_width)` → returns chunk id
    *   `reserve_chunks(N)` → prevent TOC reallocs
    *   `finalize()` → writes TOC, footer, `fdatasync()`

**Why policies**: zero runtime branching, inlining, and reuse across different sinks.

**Sample**:

```cpp
#include "archive_writer.hpp"
using Writer = charc::ArchiveWriter<
  charc::AlignTo<2*1024*1024>,   // 2 MiB alignment for chunk starts
  charc::Preallocate<2ull*1024*1024*1024>, // 2 GiB prealloc
  charc::FadviseSequential,
  charc::CRC32
>;

Writer w("dataset.char");
w.reserve_chunks(500'000); // scale metadata
w.append_raw_packed(packed0.data(), packed0.size(), elems0, width0);
w.append_raw_packed(packed1.data(), packed1.size(), elems1, width1);
// ...
w.finalize();
```

***

### 4) Reader (mmap, O(log N) lookup, CRC validation)

**Header**: `archive_reader.hpp`

*   `ArchiveReader(path)` → mmaps file read-only; validates footer/TOC bounds.
*   `chunk_count()`, `total_elems()`, `meta(i)`.
*   `chunk_region(i)` → `(ptr, len)` into the mapped file (no copy).
*   `find_chunk_by_global_index(idx)` → O(log N).
*   **CRC methods**:
    *   `validate_chunk(i)` → recompute CRC32 for that chunk and compare.
    *   `validate_range(first, last)` (single-thread).
    *   `validate_all()` (single-thread).
    *   `validate_all_parallel(threads)` (parallel verification).

**Usage**:

```cpp
#include "archive_reader.hpp"
charc::ArchiveReader r("dataset.char");

auto all = r.validate_all_parallel(); // CRCs
if (!all.mismatches.empty()) { /* handle corruption */ }

// Random lookup
uint64_t idx = 123456789;
long k = r.find_chunk_by_global_index(idx);
if (k >= 0) {
  auto& m = r.meta((size_t)k);
  auto [ptr, len] = r.chunk_region((size_t)k);
  uint64_t local = idx - m.start_index;
  // Create a zero-copy reader for RAW_PACKED chunk (extract bits at local index)
  // (you can reuse your RawPackedArray or equivalent)
}
```

***

### 5) Integration with the Multi-Sink Pipeline

*   Use a **producer → bundle (per sink payloads) → ordered writer → sinks** topology.
*   For the archive, the **payload** is simply the packed byte region + element count + bit width. In our compile-time pipeline we typically pass a `PackedPayload` struct for this sink, but if you already prepack in producers, you can store the bytes right in the bundle.

**Example sink tuple**:

```cpp
#include "multisink_pipeline.hpp"
#include "archive_writer.hpp"
#include "msink_huge_bytearray_sink.hpp"  // pure concatenation (no alignment)

using ArchiveSink = /* a small sink that wraps ArchiveWriter policies */;
using BytesSink   = msink::HugeByteArraySink<msink::NoAlignment, msink::NoPreallocate, msink::FadviseSequential>;

using Sinks = std::tuple<ArchiveSink, BytesSink>;
```

> **Important**: The **HugeByteArraySink** is a literal concatenation of input buffers. It must **not** pad or align between payloads, otherwise the resulting file can’t be treated as a single continuous byte array. Keep alignment for **structured** outputs like the archive.

***

### 6) Practical Tips

*   **Chunk sizing**: Larger chunks (tens/hundreds of MB or even GB) generally improve throughput and reduce TOC overhead.
*   **Alignment**: Use `AlignTo<4096>` for general workloads; `AlignTo<2*1024*1024>` can help scan-heavy jobs by reducing TLB misses.
*   **Preallocation**: If you can estimate final file size, `Preallocate<bytes>` reduces fragmentation and write stalls.
*   **TOC scale**: `ChunkMeta` is compact; millions of entries are feasible. Call `reserve_chunks(N)` to avoid TOC reallocations.
*   **Parallel read**: The reader is thread-safe (shared read-only mmap). Use `validate_all_parallel()` for integrity sweeps.

***

### 7) End-to-End Sample (Writer + Reader)

```cpp
// Writer:
#include "archive_writer.hpp"
#include <vector>

using Writer = charc::ArchiveWriter<
  charc::AlignTo<4096>,
  charc::Preallocate<0>,
  charc::FadviseSequential,
  charc::CRC32
>;

int main() {
  Writer w("dataset.char");
  std::vector<uint8_t> packedA = /* bit-packed data ... */;
  std::vector<uint8_t> packedB = /* bit-packed data ... */;

  w.append_raw_packed(packedA.data(), packedA.size(), /*elemsA*/ 10'000'000, /*widthA*/ 16);
  w.append_raw_packed(packedB.data(), packedB.size(), /*elemsB*/ 12'000'000, /*widthB*/ 20);
  w.finalize();
}
```

```cpp
// Reader:
#include "archive_reader.hpp"
#include <iostream>

int main() {
  charc::ArchiveReader r("dataset.char");
  auto chk = r.validate_all_parallel();
  if (!chk.mismatches.empty()) {
    std::cerr << "CRC mismatches: " << chk.mismatches.size() << "\n";
  }

  // Locate a global index and interpret the chunk
  uint64_t idx = 123'456'789ull;
  long k = r.find_chunk_by_global_index(idx);
  if (k >= 0) {
    const auto& m = r.meta((size_t)k);
    auto [ptr, len] = r.chunk_region((size_t)k);
    uint64_t local = idx - m.start_index;
    // e.g., wrap in your RawPackedArray and read value at 'local'
  }
}
```

***

### 8) Why This Archive Works Well

*   **Write path**: sequential appends, OS hints, optional preallocation, minimal metadata.
*   **Read path**: one mmap; O(log N) TOC lookup; direct pointer to chunk bytes.
*   **Integrity**: per-chunk CRC32 to isolate/diagnose corruption.
*   **Extensibility**: add compression or encryption as policies; add side indexes if you need sub-chunk point queries.

***

### 9) Related Headers (cross-reference)

*   **Policies**: `sink_policies.hpp`
*   **Archive Writer**: `archive_writer.hpp` (policy-based)
*   **Archive Reader**: `archive_reader.hpp` (CRC validation)
*   **Huge Byte Array Sink**: `byte_writer.hpp` (pure concatenation)
*   **SDSL Int Vector Sinks**: `sdsl_writer.hpp` (on-disk via `int_vector_buffer<W>`)
*   **Sink Manager**: `sink_manager.hpp` (bundle, pool, queue, MultiSink, OrderedWriter)

***
