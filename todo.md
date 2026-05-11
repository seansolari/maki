# To do:

*Last updated*: Mon 11 May, 2026

## High priority

  - Python bindings for graph creation and querying.

## Low priority

  - Compress metacolour representation.
  - Use rrna genes hits to estimate genome size. In cxfy: (abundance of rrna * estimated genome size from rrna) to estimate how many read you expect, compared to how many reads have been classified. tells you whether unclassified reads belong to a bacteria you have a signal for, vs something you dont have a signal for.
  - Add Newick support if necessary.

## Planned updates

### Mon 11 May, 2026 - Python bindings

  - Coloured graph bindings
  - Weighted graph bindings

### Tue 28 Apr. 2026 - Migrate counting graph creation from FastQ

  - **DONE** - Migrate source code
    - **DONE** - Migrate read parsing code
    - **DONE** - Implement competitive terminal extraction from read data and counting graph construction in build module
      - **DONE** - Implement CountsBuffer
      - **DONE** - Implement pushNode for terminals (interleave_buffers.[hc]pp)
      - **DONE** - Implement TerminalsGate (construct_terminals.[hc]pp)
      - **DONE** - Implement SuffixwiseTerminals (construct_wdbg.[hc]pp, compare to construct_cdbg.[hc]pp)
        - **DONE** - Implement TerminalwiseSuffix using SuffixWise (construct_wdbg.[hc]pp, compare to construct_cdbg.[hc]pp)
        - **DONE** - Double check suffix plan for suffix-wise terminal buffer (constructSuffixPlan needs to start at s, not k-s?)
        - **DONE** - Refactor graph finalisation into dbg (construct_common.hpp) and [cw]dbg (construct_[cw]dbg.[hc]pp)
  - **DONE** - Migrate fastq parsing unit tests
  - **DONE** - Implement counting graph tests
    - **DONE** - wap_vector
    - **DONE** - CountsBuffer tests
    - **DONE** - TerminalsGate tests
      - **DONE** - LongSuffixGate tests
      - **DONE** - extractTerminalsDense tests
    - **DONE** - Terminal buffer suffix-specific k-mer extraction
    - **DONE** - PushNode tests
    - **DONE** - Graph k-mer counts test

### Wed 22 Apr. 2026 - Migrate merge algorithm

  - **DONE** - Migrate small merge
  - **DONE** - Implement large merge
    - **DONE** - InterleaveRange
    - **DONE** - NextIteration
    - **DONE** - ELMMergeLarge(...)
    - **DONE** - FetchRange
  - **DONE** - Migrate small merge tests
  - **DONE** - Implement large merge tests
    - **DONE** - Create equivalent tests to small merge tests
  - **DONE** - Implement generic merge test suite

### Tue 27 Jan. 2026 - Unit Testing

  - **DONE** - Colour Archive
    - **DONE** - Bit-packed array.
    - **DONE** - Writing to the archive, with various policies.
    - **DONE** - Reading from an archive.
    - **DONE** - Multiple threads reading.
  - **DONE** - Byte writer
  - **DONE** - SDSL writer
    - **DONE** - In memory.
    - **DONE** - On disk.
  - **DONE** - Vector writer
  - **DONE** - Sink manager
    - **DONE** - Factory
    - **DONE** - Pool
      - **DONE** - Preallocate
      - **DONE** - Get
      - **DONE** - Put
    - **DONE** - Multi-threaded pipeline
  - **DONE** - Tuple Map
  - **DONE** - Colour Mapping
  - **DONE** - Buffer Value
  - **DONE** - Meta Colours
  - **DONE** - FASTA parsing
  - **DONE** - Suffix tests
  - **DONE** - k-mer and terminal buffers
    - **DONE** - migrate test_terminals and test_kmers from archive
  - **DONE** - flushing tests
  - **DONE** - graph IO tests
