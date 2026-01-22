
#pragma once
#include <cstdint>
#include <vector>

#include "seq_io.hpp"

class SequenceContainer
{
  virtual std::size_t numTerminals(std::size_t k_) const = 0;
};
