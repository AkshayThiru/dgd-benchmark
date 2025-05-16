#include <dgd/data_types.h>
#include <dgd/utils.h>

#include <string>
#include <vector>

#include "dsf/dsf_collision.h"
#include "dsf/dsf_interface.h"

int main() {
  std::vector<dgd::Vec3f> vert;
  double inradius, margin{0.0};
  constexpr unsigned int exp{16};

  inradius = dsf::LoadOBJ("../assets/rock_lowpoly.obj", vert);

  auto set1 = new dsf::VDSFInterface<exp>(vert, inradius, margin);
  auto set2 = new dsf::VDSFInterface<exp>(vert, inradius, margin);

  set1->VDSFPtr()->PrintInfo();

  dgd::Transform3f tf1, tf2;
  dgd::RandomRigidBodyTransform<3>(-5.0, 5.0, tf1);
  dgd::RandomRigidBodyTransform<3>(-5.0, 5.0, tf2);

  const dsf::Settings settings{};
  dsf::Output out{};
  double gd{dsf::GrowthDistance(set1->VDSFPtr(), tf1, set2->VDSFPtr(), tf2, out,
                                settings)};
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
