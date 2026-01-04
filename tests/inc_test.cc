#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "dgd/geometry/3d/mesh.h"
#include "dgd/mesh_loader.h"
#include "dgd/utils/random.h"
#include "inc/data_types.h"
#include "inc/incremental.h"
#include "inc/polyhedron.h"
#include "inc/solution_error.h"

namespace {

void GetConvexSets(const std::string& file, inc::Polyhedron*& seth,
                   dgd::Mesh*& setv) {
  dgd::MeshLoader ml{};
  ml.LoadObj(file);

  // Construct halfspace set.
  std::vector<inc::Vec3> vec;
  std::vector<double> offset;
  std::vector<int> graph;
  inc::Vec3 center;
  ml.MakeFacetGraph(vec, offset, graph, center);
  for (std::size_t i = 0; i < vec.size(); ++i) {
    offset[i] += vec[i].dot(center);
  }
  seth = new inc::Polyhedron(vec, offset, graph);

  // Construct vertex set.
  const double inradius = ml.ComputeInradius(vec, offset, inc::Vec3::Zero());
  vec.clear();
  graph.clear();
  ml.MakeVertexGraph(vec, graph);
  for (auto& v : vec) v -= center;
  setv = new dgd::Mesh(vec, graph, inradius, 0.0);
}

}  // namespace

int main() {
  // Mesh sets.
  inc::Polyhedron *set1, *set2;
  dgd::Mesh *mesh1, *mesh2;
  GetConvexSets("../assets/006_mustard_bottle.obj", set1, mesh1);
  GetConvexSets("../assets/006_mustard_bottle.obj", set2, mesh2);

  mesh1->PrintInfo();

  dgd::Rng rng;
  rng.SetRandomSeed();

  // Rigid body transforms.
  inc::Transform3 tf1, tf2;
  inc::Vec3 v, euler;
  inc::Rotation3 dR;

  // Growth distance.
  const int ncold = 100;
  const int nwarm = 999;
  const double dt = 0.1;
  int niter_cold = 0, niter_warm = 0;

  inc::GrowthDistanceSolver solver{};
  inc::Output out{};
  inc::SolutionError err;
  double max_gap = 0.0;
  for (int i = 0; i < ncold; ++i) {
    rng.RandomTransform(-5.0, 5.0, tf1);
    rng.RandomTransform(-5.0, 5.0, tf2);
    v = dgd::Vec3r::Random();
    euler = dgd::Vec3r::Random() * dgd::kPi / 2.0 * dt;
    dR = dgd::EulerToRotation(euler);
    for (int j = 0; j < nwarm + 1; ++j) {
      solver.GrowthDistance(set1, tf1, set2, tf2, out, (j > 0));
      if (j > 0) {
        niter_warm += out.iter;
      } else {
        niter_cold += out.iter;
      }
      err = inc::ComputeSolutionError(mesh1, tf1, mesh2, tf2, out);
      max_gap = std::max(max_gap, err.prim_dual_gap);
      inc::Affine(tf1) += v * dt;
      inc::Linear(tf1) *= dR;
    }
  }

  const double avg_iter_cold = (1.0 * niter_cold) / ncold;
  const double avg_iter_warm = (1.0 * niter_warm) / (ncold * nwarm);
  std::cout << "Maximum primal-dual gap: " << max_gap << std::endl
            << "Avg. iterations (cold) : " << avg_iter_cold << std::endl
            << "Avg. iterations (warm) : " << avg_iter_warm << std::endl;

  delete set1;
  delete set2;
  delete mesh1;
  delete mesh2;
  return EXIT_SUCCESS;
}
