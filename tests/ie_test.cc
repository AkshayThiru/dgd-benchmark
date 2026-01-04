#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "dgd/data_types.h"
#include "dgd/geometry/3d/cone.h"
#include "dgd/geometry/3d/mesh.h"
#include "dgd/mesh_loader.h"
#include "dgd/utils/random.h"
#include "ie/internal_expanding.h"

int main() {
  // Cone set.
  const double ha = dgd::kPi / 6.0, radius = 1.0;
  const double height = radius / std::tan(ha);
  dgd::ConvexSet<3>* set1 = new dgd::Cone(radius, height, 0.1);

  // Mesh set.
  dgd::MeshLoader ml{};
  ml.LoadObj("../assets/006_mustard_bottle.obj");
  std::vector<dgd::Vec3r> vert;
  std::vector<int> graph;
  ml.MakeVertexGraph(vert, graph);
  dgd::Vec3r center;
  const double inradius = ml.ComputeInradius(center);
  for (auto& v : vert) v -= center;
  dgd::ConvexSet<3>* set2 = new dgd::Mesh(vert, graph, inradius, 0.1);

  dgd::Rng rng;
  rng.SetRandomSeed();

  // Rigid body transform.
  dgd::Transform3r tf1, tf2;
  rng.RandomTransform(-5.0, 5.0, tf1);
  rng.RandomTransform(-5.0, 5.0, tf2);

  // Growth distance.
  const ie::Settings settings{};
  ie::Output out{};
  const double gd{ie::GrowthDistance(set1, tf1, set2, tf2, settings, out)};
  auto err = ie::ComputeSolutionError(set1, tf1, set2, tf2, out);

  std::cout << "Growth distance               : " << gd << std::endl
            << "Primal infeasibility error (m): " << err.prim_infeas_err
            << std::endl
            << "Dual infeasibility error      : " << err.dual_infeas_err
            << std::endl
            << "Primal dual gap               : " << err.prim_dual_gap
            << std::endl;

  delete set1;
  delete set2;
  return EXIT_SUCCESS;
}
