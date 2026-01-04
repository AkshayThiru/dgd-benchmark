#include <iostream>
#include <string>
#include <vector>

#include "dcf/dcf_collision.h"
#include "dcf/dsf_interface.h"
#include "dgd/mesh_loader.h"
#include "dgd/utils/random.h"

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

  set1->PrintInfo();

  dgd::Rng rng;
  rng.SetRandomSeed();

  // Rigid body transforms.
  dcf::Transform3 tf1, tf2;
  rng.RandomTransform(-5.0, 5.0, tf1);
  rng.RandomTransform(-5.0, 5.0, tf2);

  // Growth distance.
  dcf::Settings settings{};
  dcf::Output out{};
  double gd = dcf::GrowthDistance(set1->VDSFPtr(), tf1, set2->VDSFPtr(), tf2,
                                  settings, out);
  auto err = dcf::ComputeSolutionError(set1->VDSFPtr(), tf1, set2->VDSFPtr(),
                                       tf2, out);

  std::cout << "Growth distance               : " << gd << std::endl
            << "Primal infeasibility error (m): " << err.prim_infeas_err
            << std::endl
            << "Dual infeasibility error      : " << err.dual_infeas_err
            << std::endl
            << "Primal dual gap               : " << err.prim_dual_gap
            << std::endl
            << "Iterations                    : " << out.iter << std::endl;

  delete set1;
  delete set2;
  return EXIT_SUCCESS;
}
