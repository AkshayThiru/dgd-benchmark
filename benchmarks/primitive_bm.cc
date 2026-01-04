#include <string>

#include "dgd/geometry/convex_set.h"
#include "helpers/benchmark_interface.h"
#include "helpers/benchmark_result.h"
#include "internal_helpers/filesystem_utils.h"
#include "internal_helpers/math_utils.h"

// Constants.
const double position_lim = 5.0;
const double dx_max = 0.1, ang_max = dgd::kPi / 18.0;
const int npair = 1000;
const int npose_c = 100, npose_w = 100;
const int ncold = 100, nwarm = 100;

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

  bench::BenchmarkInterface interface(ncold, nwarm);
  interface.SetRngSeed();

  using bench::DgdBcSolverType;
  using bench::DgdSolverType;

  dgd::Transform3r tf1, tf2;

  /**
   * Cold-start benchmarks.
   */
  bench::BenchmarkResultArray res_ie_c(npair * npose_c);
  bench::BenchmarkResultArray res_dgd_c_cp_cr(npair * npose_c);
  bench::BenchmarkResultArray res_dgd_c_cp_lu(npair * npose_c);
  bench::BenchmarkResultArray res_dgd_c_trn(npair * npose_c);
  for (int i = 0; i < npair; ++i) {
    const auto set1 = interface.RandomCurvedPrimitiveSet();
    const auto set2 = interface.RandomCurvedPrimitiveSet();
    for (int j = 0; j < npose_c; ++j) {
      dgd::bench::SetRandomTransforms(interface.rng(), tf1, tf2, -position_lim,
                                      position_lim);
      // IE benchmark.
      interface.IeColdStart(set1.get(), tf1, set2.get(), tf2, res_ie_c);
      // DGD benchmark (cutting plane, Cramer's rule).
      interface.DgdColdStart<DgdSolverType::CuttingPlane,
                             DgdBcSolverType::Cramer>(
          set1.get(), tf1, set2.get(), tf2, res_dgd_c_cp_cr);
      // DGD benchmark (cutting plane, LU).
      interface.DgdColdStart<DgdSolverType::CuttingPlane, DgdBcSolverType::LU>(
          set1.get(), tf1, set2.get(), tf2, res_dgd_c_cp_lu);
      // DGD benchmark (trust region Newton).
      interface.DgdColdStart<DgdSolverType::TrustRegionNewton,
                             DgdBcSolverType::LU>(set1.get(), tf1, set2.get(),
                                                  tf2, res_dgd_c_trn);
    }
  }

  std::cout << "Internal expanding (cold-start):" << std::endl;
  res_ie_c.SaveToFile(log_path + "primitive_bm__cold_ie.feather");
  res_ie_c.PrintStatistics();

  std::cout << "DGD (cold-start, cutting plane, Cramer's rule):" << std::endl;
  res_dgd_c_cp_cr.SaveToFile(log_path +
                             "primitive_bm__cold_dgd_cp_cramer.feather");
  res_dgd_c_cp_cr.PrintStatistics();

  std::cout << "DGD (cold-start, cutting plane, LU):" << std::endl;
  res_dgd_c_cp_lu.SaveToFile(log_path + "primitive_bm__cold_dgd_cp_lu.feather");
  res_dgd_c_cp_lu.PrintStatistics();

  std::cout << "DGD (cold-start, trust region Newton):" << std::endl;
  res_dgd_c_trn.SaveToFile(log_path + "primitive_bm__cold_dgd_trn.feather");
  res_dgd_c_trn.PrintStatistics();

  /**
   * Warm-start benchmarks.
   */
  using dgd::WarmStartType;

  dgd::Vec3r dx;
  dgd::Rotation3r drot;
  bench::BenchmarkResultArray res_dgd_w_p_cr(npair * npose_w * nwarm);
  bench::BenchmarkResultArray res_dgd_w_p_lu(npair * npose_w * nwarm);
  bench::BenchmarkResultArray res_dgd_w_d_cr(npair * npose_w * nwarm);
  bench::BenchmarkResultArray res_dgd_w_d_lu(npair * npose_w * nwarm);
  for (int i = 0; i < npair; ++i) {
    const auto set1 = interface.RandomCurvedPrimitiveSet();
    const auto set2 = interface.RandomCurvedPrimitiveSet();
    for (int j = 0; j < npose_w; ++j) {
      dgd::bench::SetRandomTransforms(interface.rng(), tf1, tf2, -position_lim,
                                      position_lim);
      dgd::bench::SetRandomDisplacement(interface.rng(), dx, drot, dx_max,
                                        ang_max);
      // DGD benchmark (primal warm start, cutting plane, Cramer's rule).
      interface.DgdWarmStart<DgdSolverType::CuttingPlane,
                             DgdBcSolverType::Cramer>(
          set1.get(), tf1, set2.get(), tf2, dx, drot, res_dgd_w_p_cr,
          WarmStartType::Primal);
      // DGD benchmark (primal warm start, cutting plane, LU).
      interface.DgdWarmStart<DgdSolverType::CuttingPlane, DgdBcSolverType::LU>(
          set1.get(), tf1, set2.get(), tf2, dx, drot, res_dgd_w_p_lu,
          WarmStartType::Primal);
      // DGD benchmark (dual warm start, cutting plane, Cramer's rule).
      interface.DgdWarmStart<DgdSolverType::CuttingPlane,
                             DgdBcSolverType::Cramer>(
          set1.get(), tf1, set2.get(), tf2, dx, drot, res_dgd_w_d_cr,
          WarmStartType::Dual);
      // DGD benchmark (dual warm start, cutting plane, LU).
      interface.DgdWarmStart<DgdSolverType::CuttingPlane, DgdBcSolverType::LU>(
          set1.get(), tf1, set2.get(), tf2, dx, drot, res_dgd_w_d_lu,
          WarmStartType::Dual);
    }
  }

  std::cout << "DGD (primal warm start, cutting plane, Cramer's rule):"
            << std::endl;
  res_dgd_w_p_cr.SaveToFile(log_path +
                            "primitive_bm__warm_dgd_primal_cp_cramer.feather");
  res_dgd_w_p_cr.PrintStatistics();

  std::cout << "DGD (primal warm start, cutting plane, LU):" << std::endl;
  res_dgd_w_p_lu.SaveToFile(log_path +
                            "primitive_bm__warm_dgd_primal_cp_lu.feather");
  res_dgd_w_p_lu.PrintStatistics();

  std::cout << "DGD (dual warm start, cutting plane, Cramer's rule):"
            << std::endl;
  res_dgd_w_d_cr.SaveToFile(log_path +
                            "primitive_bm__warm_dgd_dual_cp_cramer.feather");
  res_dgd_w_d_cr.PrintStatistics();

  std::cout << "DGD (dual warm start, cutting plane, LU):" << std::endl;
  res_dgd_w_d_lu.SaveToFile(log_path +
                            "primitive_bm__warm_dgd_dual_cp_lu.feather");
  res_dgd_w_d_lu.PrintStatistics();
}
