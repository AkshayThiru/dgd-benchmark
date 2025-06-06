#include "dgd_benchmark/dsf/dsf_interface.h"

#include <stdexcept>
#include <string>
#include <vector>

#include "dgd/mesh_loader.h"

namespace dsf {

double LoadObj(const std::string& file, std::vector<Vec3>& vert) {
  dgd::MeshLoader ml{};
  ml.LoadObj(file);

  vert.clear();
  std::vector<int> graph;
  bool valid = ml.MakeVertexGraph(vert, graph);
  if (!valid) throw std::runtime_error("Qhull error");
  Vec3 center;
  const double inradius = ml.ComputeInradius(center);
  for (auto& v : vert) v -= center;
  return inradius;
}

}  // namespace dsf
