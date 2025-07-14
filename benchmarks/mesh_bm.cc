#include <string>

#include "dgd/geometry/3d/mesh.h"
#include "dgd/geometry/convex_set.h"
#include "helpers/benchmark_interface.h"
#include "helpers/benchmark_result.h"
#include "internal_helpers/utils.h"

// Constants.
const double position_lim = 5.0;
const double dx_max = 0.1, ang_max = dgd::kPi / 18.0;
const int npair = 1000;
const int npose_c = 100, npose_w = 100;
const int ncold = 100, nwarm = 100;

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

  internal::BenchmarkInterface interface(ncold, nwarm);
  interface.LoadMeshesFromObjFiles(filenames);
  if (interface.nmeshes() == 0) {
    std::cerr << "No meshes were loaded" << std::endl;
    return EXIT_FAILURE;
  }

  dgd::ConvexSet<3>*set1, *set2;
  int set1_idx, set2_idx;
  dgd::Transform3r tf1, tf2;

  /**
   * Cold-start benchmarks.
   */
  internal::BenchmarkResultArray res_inc_c(npair * npose_c);
  internal::BenchmarkResultArray res_ie_c(npair * npose_c);
  internal::BenchmarkResultArray res_dgd_c(npair * npose_c);
  for (int i = 0; i < npair; ++i) {
    set1_idx = interface.RandomMeshIndex();
    set2_idx = interface.RandomMeshIndex();
    set1 = interface.meshes()[set1_idx].get();
    set2 = interface.meshes()[set2_idx].get();
    for (int j = 0; j < npose_c; ++j) {
      dgd::internal::SetRandomTransform<3>(tf1, tf2, -position_lim,
                                           position_lim);
      // INC benchmark.
      interface.IncColdStart(set1_idx, tf1, set2_idx, tf2, res_inc_c);
      // IE benchmark.
      interface.IeColdStart(set1, tf1, set2, tf2, res_ie_c);
      // DGD benchmark.
      interface.DgdColdStart(set1, tf1, set2, tf2, res_dgd_c);
    }
  }

  std::cout << "Incremental (cold-start):" << std::endl;
  res_inc_c.SaveToFile(log_path + "mesh_bm__cold_inc.feather");
  res_inc_c.PrintStatistics();

  std::cout << "Internal expanding (cold-start):" << std::endl;
  res_ie_c.SaveToFile(log_path + "mesh_bm__cold_ie.feather");
  res_ie_c.PrintStatistics();

  std::cout << "DGD (cold-start):" << std::endl;
  res_dgd_c.SaveToFile(log_path + "mesh_bm__cold_dgd.feather");
  res_dgd_c.PrintStatistics();

  /**
   * Warm-start benchmarks.
   */
  dgd::Vec3r dx;
  dgd::Rotation3r drot;
  internal::BenchmarkResultArray res_inc_w(npair * npose_w * nwarm);
  internal::BenchmarkResultArray res_dgd_w(npair * npose_w * nwarm);
  for (int i = 0; i < npair; ++i) {
    set1_idx = interface.RandomMeshIndex();
    set2_idx = interface.RandomMeshIndex();
    set1 = interface.meshes()[set1_idx].get();
    set2 = interface.meshes()[set2_idx].get();
    for (int j = 0; j < npose_w; ++j) {
      dgd::internal::SetRandomTransform<3>(tf1, tf2, -position_lim,
                                           position_lim);
      dgd::internal::SetRandomScrew(dx, drot, dx_max, ang_max);
      // INC benchmark.
      interface.IncWarmStart(set1_idx, tf1, set2_idx, tf2, dx, drot, res_inc_w);
      // DGD benchmark.
      interface.DgdWarmStart(set1, tf1, set2, tf2, dx, drot, res_dgd_w);
    }
  }

  std::cout << "Incremental (warm-start):" << std::endl;
  res_inc_w.SaveToFile(log_path + "mesh_bm__warm_inc.feather");
  res_inc_w.PrintStatistics();

  std::cout << "DGD (warm-start):" << std::endl;
  res_dgd_w.SaveToFile(log_path + "mesh_bm__warm_dgd.feather");
  res_dgd_w.PrintStatistics();
}
