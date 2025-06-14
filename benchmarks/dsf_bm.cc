#include <string>

#include "dgd/geometry/convex_set.h"
#include "helpers/benchmark_interface.h"
#include "helpers/benchmark_result.h"
#include "internal_helpers/utils.h"

// Constants.
const double position_lim = 5.0;
const int npair = 100, npose = 100;
const int ncold = 1;
// DSF settings.
const int ie_iter = 8;
const int max_iter = 30;

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <asset_path> <log_path>"
              << std::endl;
    return EXIT_FAILURE;
  }

  const std::string asset_path = argv[1], log_path = argv[2];
  if (!dgd::internal::IsValidDirectory(asset_path) ||
      !dgd::internal::IsValidDirectory(log_path)) {
    std::cerr << "Invalid asset or log directory" << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<std::string> filenames;
  dgd::internal::GetObjFileNames(asset_path, filenames);

  internal::BenchmarkInterface interface(ncold, 0);
  interface.LoadMeshesFromObjFiles(filenames);
  if (interface.nmeshes() < 2) {
    std::cerr << "At least two .obj files are required" << std::endl;
    return EXIT_FAILURE;
  }
  dsf::Settings settings{};
  settings.ie_iter = ie_iter;
  settings.max_iter = max_iter;
  interface.SetDsfSettings(settings);

  internal::BenchmarkResultArray res_dsf(npair * npose);
  internal::BenchmarkResultArray res_ie(npair * npose);
  internal::BenchmarkResultArray res_dgd(npair * npose);
  dgd::ConvexSet<3>*set1, *set2;
  int set1_idx, set2_idx;
  dgd::Transform3r tf1, tf2;
  for (int i = 0; i < npair; ++i) {
    set1_idx = interface.RandomMeshIndex();
    do {
      set2_idx = interface.RandomMeshIndex();
    } while (set2_idx == set1_idx);
    set1 = interface.vdsfs()[set1_idx].get();
    set2 = interface.vdsfs()[set2_idx].get();
    for (int j = 0; j < npose; ++j) {
      dgd::internal::SetRandomTransform<3>(tf1, tf2, -position_lim,
                                           position_lim);
      // DSF benchmark.
      interface.DsfColdStart(set1_idx, tf1, set2_idx, tf2, res_dsf);
      // IE benchmark.
      interface.IeColdStart(set1, tf1, set2, tf2, res_ie);
      // DGD benchmark.
      interface.DgdColdStart(set1, tf1, set2, tf2, res_dgd);
    }
  }

  std::cout << "Differentiable support function:" << std::endl;
  res_dsf.SaveToFile(log_path + "dsf_bm__cold_dsf.feather");
  res_dsf.PrintStatistics();

  std::cout << "Internal expanding:" << std::endl;
  res_ie.SaveToFile(log_path + "dsf_bm__cold_ie.feather");
  res_ie.PrintStatistics();

  std::cout << "DGD:" << std::endl;
  res_dgd.SaveToFile(log_path + "dsf_bm__cold_dgd.feather");
  res_dgd.PrintStatistics();
}
