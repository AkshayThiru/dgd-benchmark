#include <iostream>
#include <string>
#include <vector>

#include "dcf/dcf_collision.h"
#include "dcf/dsf_interface.h"
#include "dgd/mesh_loader.h"
#include "dgd/utils.h"

double LoadObj(std::string filename, std::vector<dcf::Vec3>& vert) {
  dgd::MeshLoader ml{};
  ml.LoadObj(filename);

  vert.clear();
  std::vector<int> graph;
  const bool valid = ml.MakeVertexGraph(vert, graph);
  if (!valid) throw std::runtime_error("Qhull error");
  dcf::Vec3 center;
  const double inradius = ml.ComputeInradius(center);
  for (auto& v : vert) v -= center;
  return inradius;
}

int main() {
  // Mesh sets.
  std::vector<dcf::Vec3> vert;
  double inradius, margin = 0.0;
  constexpr unsigned int exp = 16;

  inradius = LoadObj("../assets/006_mustard_bottle.obj", vert);
  auto set1 = new dcf::VDSFInterface<exp>(vert, inradius, margin);
  inradius = LoadObj("../assets/006_mustard_bottle.obj", vert);
  auto set2 = new dcf::VDSFInterface<exp>(vert, inradius, margin);

  set1->VDSFPtr()->PrintInfo();

  // Rigid body transforms.
  dcf::Transform3 tf1, tf2;
  dgd::RandomRigidBodyTransform<3>(-5.0, 5.0, tf1);
  dgd::RandomRigidBodyTransform<3>(-5.0, 5.0, tf2);

  // Growth distance.
  dcf::Settings settings{};
  dcf::Output out{};
  double gd = dcf::GrowthDistance(set1->VDSFPtr(), tf1, set2->VDSFPtr(), tf2,
                                  settings, out);
  auto err = dcf::ComputeSolutionError(set1->VDSFPtr(), tf1, set2->VDSFPtr(),
                                       tf2, out);

  std::cout << "Growth distance         : " << gd << std::endl
            << "Primal feasibility error: " << err.prim_feas_err << " m"
            << std::endl
            << "Dual feasibility error  : " << err.dual_feas_err << std::endl
            << "Primal dual gap         : " << err.prim_dual_gap << std::endl
            << "Iterations              : " << out.iter << std::endl;

  delete set1;
  delete set2;
  return EXIT_SUCCESS;
}
