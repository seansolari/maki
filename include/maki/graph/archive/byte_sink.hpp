
#pragma once

class SequentialFileSink
{
public:
  explicit SequentialFileSink(const std::string &path)
  {
    f_ = std::fopen(path.c_str(), "wb");
    if (!f_)
      throw std::runtime_error("SequentialFileSink: fopen failed");
    // Optional: enlarge stdio buffer
    // setvbuf(f_, nullptr, _IOFBF, 1<<20);
  }

  ~SequentialFileSink() override
  {
    if (f_)
      std::fclose(f_);
  }

  void write_chunk(const PackedChunk &ch) override
  {
    if (ch.bytes.empty())
      return;
    size_t w = std::fwrite(ch.bytes.data(), 1, ch.bytes.size(), f_);
    if (w != ch.bytes.size())
      throw std::runtime_error("SequentialFileSink: short write");
  }
  
  void finalize() override
  {
    if (f_)
      std::fflush(f_);
  }

private:
  std::FILE *f_ = nullptr;
};