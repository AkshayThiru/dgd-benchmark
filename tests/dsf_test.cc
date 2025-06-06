#include <iostream>
#include <string>
#include <vector>

#include "dgd/utils.h"
#include "dgd_benchmark/dsf/dsf_collision.h"
#include "dgd_benchmark/dsf/dsf_interface.h"

int main() {
  // Mesh sets.
  std::vector<dsf::Vec3> vert;
  double inradius, margin = 0.0;
  constexpr unsigned int exp = 16;

  inradius = dsf::LoadObj("../assets/rock_lowpoly.obj", vert);

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
