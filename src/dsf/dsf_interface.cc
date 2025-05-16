#include "dsf/dsf_interface.h"

#include <dgd/mesh_loader.h>

#include <stdexcept>

namespace dsf {

double LoadOBJ(const std::string& file, std::vector<dgd::Vec3f>& vert) {
  dgd::MeshLoader ml{};
  ml.LoadOBJ(file);

  vert.clear();
  std::vector<int> graph;
  bool valid{ml.MakeVertexGraph(vert, graph)};
  if (!valid) throw std::runtime_error("Qhull error");
  dgd::Vec3f center;
  const double inradius{ml.ComputeInradius(center)};
  for (auto& v : vert) v -= center;
  return inradius;
}

}  // namespace dsf
