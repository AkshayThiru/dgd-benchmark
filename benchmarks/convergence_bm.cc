#include <memory>

#include "dgd/geometry/geometry_3d.h"
#include "dgd/output.h"
#include "dgd/settings.h"
#include "dgd/solvers/bundle_scheme_impl.h"
#include "helpers/convergence_result.h"
#include "internal_helpers/filesystem_utils.h"
#include "internal_helpers/set_generator.h"

// Constants.
const double position_lim = 5.0;

const int npair_p = 10000, npair_m = 1000;
const int npose = 100;
const int nvert_m = 1000;

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

  dgd::bench::ConvexSetFeatureRange fr{};
  dgd::bench::ConvexSetGenerator gen(fr);
  dgd::Rng rng;
  gen.SetRngSeed();
  rng.SetSeed();

  using C = dgd::ConvexSet<3>;
  std::shared_ptr<C> set1, set2;
  dgd::Transform3r tf1, tf2;

  dgd::Settings settings;
  dgd::Output<3> out;

  auto growth_distance =
      dgd::GrowthDistanceCpTpl<3, C, C, dgd::BcSolverType::kLU>;

  /*
   * Primitive types.
   */

  {
    bench::IterationStats stats(settings.max_iter);
    bench::ConvergenceResultArray res_prim(settings.max_iter);

    for (int i = 0; i < npair_p; ++i) {
      set1 = gen.GetRandomPrimitive3DSet();
      set2 = gen.GetRandomPrimitive3DSet();
      for (int j = 0; j < npose; ++j) {
        rng.RandomTransform(-position_lim, position_lim, tf1);
        rng.RandomTransform(-position_lim, position_lim, tf2);
        growth_distance(set1.get(), tf1, set2.get(), tf2, settings, out, false);
        if (out.status == dgd::SolutionStatus::IllConditionedInputs) {
          set1->PrintInfo();
          set2->PrintInfo();
          throw std::runtime_error("Ill-conditioned input convex sets.");
        }
        if (out.status == dgd::SolutionStatus::CoincidentCenters) {
          for (auto& ub : out.ubs) ub = 1.0;
          for (auto& lb : out.lbs) lb = 1.0;
        }
        stats.AddResults(out.iter, out.ubs, out.lbs);
      }
    }

    res_prim.ComputeStatistics(stats);
    res_prim.SaveToFile(log_path + "convergence_bm__primitive.feather");
  }

  /*
   * Mesh types.
   */

  auto convergence_bm_mesh = [&](dgd::Real skew,
                                 bench::ConvergenceResultArray& res,
                                 const std::string& filename) -> void {
    bench::IterationStats stats(settings.max_iter);
    for (int i = 0; i < npair_m; ++i) {
      try {
        set1 = gen.GetRandomEllipsoidalPolytope<dgd::Mesh>(nvert_m, skew);
        set2 = gen.GetRandomEllipsoidalPolytope<dgd::Mesh>(nvert_m, skew);
      } catch (const std::exception& e) {
        continue;
      }
      for (int j = 0; j < npose; ++j) {
        rng.RandomTransform(-position_lim, position_lim, tf1);
        rng.RandomTransform(-position_lim, position_lim, tf2);
        growth_distance(set1.get(), tf1, set2.get(), tf2, settings, out, false);
        if (out.status == dgd::SolutionStatus::IllConditionedInputs) {
          set1->PrintInfo();
          set2->PrintInfo();
          throw std::runtime_error("Ill-conditioned input convex sets.");
        }
        if (out.status == dgd::SolutionStatus::MaxIterReached) {
          throw std::runtime_error("Max iterations reached");
        }
        if (out.status == dgd::SolutionStatus::CoincidentCenters) {
          for (auto& ub : out.ubs) ub = 1.0;
          for (auto& lb : out.lbs) lb = 1.0;
        }
        stats.AddResults(out.iter, out.ubs, out.lbs);
      }
    }

    res.ComputeStatistics(stats);
    res.SaveToFile(log_path + filename);
  };

  bench::ConvergenceResultArray res_mesh(settings.max_iter);

  convergence_bm_mesh(1.0, res_mesh, "convergence_bm__mesh_1.feather");
  convergence_bm_mesh(1e-1, res_mesh, "convergence_bm__mesh_1e-1.feather");
  convergence_bm_mesh(1e-2, res_mesh, "convergence_bm__mesh_1e-2.feather");
}
