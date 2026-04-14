
#include <fstream>
#include <iostream>

inline int writeToFna(const std::string &data, const std::string &fna) {
  std::ofstream output_file(fna);

  // Check if the file was successfully opened
  if (output_file.is_open()) {
    // Write the string data to the file using the insertion operator (<<)
    output_file << data;

    // Close the file
    output_file.close();
  } else {
    std::cerr << "Error: Unable to open the file." << std::endl;
    return 1; // Return an error code
  }

  return 0;
}