#ifndef DGD_BENCHMARK_FILE_IO_H_
#define DGD_BENCHMARK_FILE_IO_H_

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "dgd_benchmark/benchmark_result.h"

// Checks if the given path is a valid directory.
inline bool IsValidDirectory(const std::string& path) {
  if (!std::filesystem::exists(path) || !std::filesystem::is_directory(path)) {
    std::cerr << "Error: The specified path is not a valid directory"
              << std::endl;
    return false;
  }
  return true;
}

// Gets the .obj file names from the given folder path.
inline void GetObjFileNames(const std::string& folder_path,
                            std::vector<std::string>& obj_filenames) {
  for (const auto& entry : std::filesystem::directory_iterator(folder_path)) {
    if (entry.is_regular_file() && entry.path().extension() == ".obj") {
      obj_filenames.push_back(folder_path + entry.path().filename().string());
    }
  }
}

// Writes benchmark output to a Feather file.
bool SaveBenchmarkOutputToFile(const BenchmarkResultArray& resa,
                               const std::string& filename);

#endif  // DGD_BENCHMARK_FILE_IO_H_
