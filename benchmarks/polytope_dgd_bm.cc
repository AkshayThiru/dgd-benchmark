#include <cmath>
#include <string>
#include <vector>

#include "dgd/geometry/3d/mesh.h"
#include "dgd/geometry/3d/polytope.h"
#include "dgd/geometry/convex_set.h"
#include "helpers/benchmark_interface.h"
#include "helpers/benchmark_result.h"
#include "internal_helpers/filesystem_utils.h"
#include "internal_helpers/math_utils.h"

// Constants.
const double position_lim = 5.0;
const int npair = 20;
const int npose = 100;
const int ncold = 100;

const double skew = 1e-1;

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <log_path>" << std::endl;
    return EXIT_FAILURE;
  }

  const std::string log_path = argv[1];
  if (!dgd::bench::IsValidDirectory(log_path)) {
    std::cerr << "Invalid log directory" << std::endl;
    return EXIT_FAILURE;
  }

  bench::BenchmarkInterface interface(ncold, 0);
  interface.SetRngSeed();

  using bench::DgdBcSolverType;
  using bench::DgdSolverType;

  std::shared_ptr<dgd::ConvexSet<3>> set1, set2;
  dgd::Transform3r tf1, tf2;

  const double log_vert0 = std::log(10.0);
  const double log_vert1 = std::log(250.0);
  const double log_vert2 = std::log(10000.0);
  const int npts0 = 10, npts1 = 20;
  const int nnverts = npts0 + npts1;

  std::vector<int> nverts;
  double log_gap = (log_vert1 - log_vert0) / (npts0 - 1);
  for (int i = 0; i < npts0; ++i) {
    nverts.push_back(
        static_cast<int>(std::round(std::exp(log_vert0 + i * log_gap))));
  }
  log_gap = (log_vert2 - log_vert1) / npts1;
  for (int i = 0; i < npts1; ++i) {
    nverts.push_back(
        static_cast<int>(std::round(std::exp(log_vert1 + (i + 1) * log_gap))));
  }

  /**
   * Polytope benchmarks.
   */
  bench::BenchmarkResultArray res_dgd_p(nnverts * npair * npose, true);
  for (const int nvert : nverts) {
    for (int i = 0; i < npair; ++i) {
      try {
        set1 = interface.GetRandomEllipsoidalPolytope(nvert, skew);
        set2 = interface.GetRandomEllipsoidalPolytope(nvert, skew);
      } catch (const std::runtime_error& e) {
        std::cerr << "Warning: Skipping polytope benchmark for nvert = "
                  << nvert << " due to error: " << e.what() << std::endl;
        continue;
      }
      for (int j = 0; j < npose; ++j) {
        dgd::bench::SetRandomTransforms(interface.rng(), tf1, tf2,
                                        -position_lim, position_lim);
        interface.DgdColdStart<DgdSolverType::CuttingPlane,
                               DgdBcSolverType::Cramer>(
            set1.get(), tf1, set2.get(), tf2, res_dgd_p, nvert);
      }
    }
  }

  std::cout << "Polytope (cold-start, cutting plane, Cramer's rule):"
            << std::endl;
  res_dgd_p.SaveToFile(log_path + "polytope_bm__cold_dgd_polytope.feather");
  res_dgd_p.PrintStatistics();

  /**
   * Mesh benchmarks.
   */
  bench::BenchmarkResultArray res_dgd_m(nnverts * npair * npose, true);
  for (const int nvert : nverts) {
    for (int i = 0; i < npair; ++i) {
      try {
        set1 = interface.GetRandomEllipsoidalMesh(nvert, skew);
        set2 = interface.GetRandomEllipsoidalMesh(nvert, skew);
      } catch (const std::runtime_error& e) {
        std::cerr << "Warning: Skipping polytope benchmark for nvert = "
                  << nvert << " due to error: " << e.what() << std::endl;
        continue;
      }
      for (int j = 0; j < npose; ++j) {
        dgd::bench::SetRandomTransforms(interface.rng(), tf1, tf2,
                                        -position_lim, position_lim);
        interface.DgdColdStart<DgdSolverType::CuttingPlane,
                               DgdBcSolverType::Cramer>(
            set1.get(), tf1, set2.get(), tf2, res_dgd_m, nvert);
      }
    }
  }

  std::cout << "Mesh (cold-start, cutting plane, Cramer's rule):" << std::endl;
  res_dgd_m.SaveToFile(log_path + "polytope_bm__cold_dgd_mesh.feather");
  res_dgd_m.PrintStatistics();
}
