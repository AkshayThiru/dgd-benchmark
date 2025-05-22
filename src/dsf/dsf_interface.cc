#include "dsf/dsf_interface.h"

#include <dgd/mesh_loader.h>

#include <stdexcept>

namespace dsf {

double LoadOBJ(const std::string& file, std::vector<Vec3>& vert) {
  dgd::MeshLoader ml{};
  ml.LoadOBJ(file);

  vert.clear();
  std::vector<int> graph;
  bool valid{ml.MakeVertexGraph(vert, graph)};
  if (!valid) throw std::runtime_error("Qhull error");
  Vec3 center;
  const double inradius{ml.ComputeInradius(center)};
  for (auto& v : vert) v -= center;
  return inradius;
}

}  // namespace dsf
