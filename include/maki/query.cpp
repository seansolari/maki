
#include <cstdint>
#include <string>

struct BuildParameters
{
  std::string databaseFolder;
  std::string queryFile;
  std::string outputFolder;
  std::size_t minReadLength = 50;
  uint8_t s = 0;
};

int main(int argc, char const *argv[])
{
  /* code */
  return 0;
}
