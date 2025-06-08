#include <iostream>
#include <string>
#include <vector>

#include "dgd/mesh_loader.h"
#include "dgd/utils.h"
#include "dsf/dsf_collision.h"
#include "dsf/dsf_interface.h"

double LoadObj(std::string filename, std::vector<dsf::Vec3>& vert) {
  dgd::MeshLoader ml{};
  ml.LoadObj(filename);

  vert.clear();
  std::vector<int> graph;
  const bool valid = ml.MakeVertexGraph(vert, graph);
  if (!valid) throw std::runtime_error("Qhull error");
  dsf::Vec3 center;
  const double inradius = ml.ComputeInradius(center);
  for (auto& v : vert) v -= center;
  return inradius;
}

int main() {
  // Mesh sets.
  std::vector<dsf::Vec3> vert;
  double inradius, margin = 0.0;
  constexpr unsigned int exp = 16;

  inradius = LoadObj("../assets/rock_lowpoly.obj", vert);

  auto set1 = new dsf::VDSFInterface<exp>(vert, inradius, margin);
  auto set2 = new dsf::VDSFInterface<exp>(vert, inradius, margin);

  set1->VDSFPtr()->PrintInfo();

  // Rigid bosy transforms.
  dsf::Transform3 tf1, tf2;
  dgd::RandomRigidBodyTransform<3>(-5.0, 5.0, tf1);
  dgd::RandomRigidBodyTransform<3>(-5.0, 5.0, tf2);

  // Growth distance.
  const dsf::Settings settings{};
  dsf::Output out{};
  double gd = dsf::GrowthDistance(set1->VDSFPtr(), tf1, set2->VDSFPtr(), tf2,
                                  settings, out);
  auto err = dsf::ComputeSolutionError(set1->VDSFPtr(), tf1, set2->VDSFPtr(),
                                       tf2, out);

  std::cout << "Growth distance         : " << gd << std::endl
            << "Primal feasibility error: " << err.prim_feas_err << " m"
            << std::endl
            << "Dual feasibility error  : " << err.dual_feas_err << std::endl
            << "Primal dual gap         : " << err.prim_dual_gap << std::endl;

  delete set1;
  delete set2;
  return EXIT_SUCCESS;
}
