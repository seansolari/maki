
#include "maki/build/graph/archive/counts_writer.hpp"

CountBufferSink::CountBufferSink(CountBuffer *out) : data(out) {}

std::size_t CountBufferSink::write(const CountBuffer &pld) {
  data->append(pld);
  return pld.node_count();
}
