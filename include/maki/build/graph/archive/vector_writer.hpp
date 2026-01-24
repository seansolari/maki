
#pragma once
#include <cstdint>
#include <vector>

// -----------------------------------------------------------------------------
// In-memory sink: append payloads into a single std::vector<>
// -----------------------------------------------------------------------------
template <typename Int> class VectorInMemorySink {
public:
  using Payload = std::vector<Int>;

  std::uint64_t write(const Payload &vec) {
    const std::uint64_t n = vec.size();
    if (n == 0)
      return 0;
    const std::uint64_t old_n = acc_.size();
    acc_.reserve(old_n + n);
    acc_.insert(std::end(acc_), std::begin(vec), std::end(vec));
    return size();
  }

  void finalize() {
    // Nothing to do; vector is already ready for use.
  }

  const std::vector<Int> &data() const noexcept { return acc_; }
  std::vector<Int> &data() noexcept { return acc_; }

  std::uint64_t size() const noexcept { return acc_.size(); }

private:
  std::vector<Int> acc_; // accumulated data
};