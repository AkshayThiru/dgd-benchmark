#include <dgd/data_types.h>
#include <dgd/geometry/3d/cone.h>
#include <dgd/geometry/3d/mesh.h>
#include <dgd/mesh_loader.h>
#include <dgd/utils.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "ie/internal_expanding.h"

int main() {
  // Cone set.
  const double ha{dgd::kPi / 6.0}, radius{1.0};
  const double height{radius / std::tan(ha)};
  dgd::ConvexSet<3>* set1 = new dgd::Cone(radius, height, 0.1);

  // Mesh set.
  dgd::MeshLoader ml{};
  ml.LoadOBJ("../assets/rock_lowpoly.obj");
  std::vector<dgd::Vec3f> vert;
  std::vector<int> graph;
  ml.MakeVertexGraph(vert, graph);
  dgd::Vec3f center;
  const double inradius{ml.ComputeInradius(center)};
  for (auto& v : vert) v -= center;
  dgd::ConvexSet<3>* set2 = new dgd::Mesh(vert, graph, 0.1, inradius);

  // Rigid body transform.
  dgd::Transform3f tf1, tf2;
  dgd::RandomRigidBodyTransform<3>(-5.0, 5.0, tf1);
  dgd::RandomRigidBodyTransform<3>(-5.0, 5.0, tf2);

  // Growth distance.
  const ie::Settings settings{};
  ie::Output out{};
  const double gd{ie::GrowthDistance(set1, tf1, set2, tf2, settings, out)};
  auto err = ie::ComputeSolutionError(set1, tf1, set2, tf2, out);

  std::cout << "Growth distance         : " << gd << std::endl
            << "Primal feasibility error: " << err.prim_feas_err << " m"
            << std::endl
            << "Dual feasibility error  : " << err.dual_feas_err << std::endl
            << "Primal dual gap         : " << err.prim_dual_gap << std::endl;

  delete set1;
  delete set2;
  return EXIT_SUCCESS;
}
