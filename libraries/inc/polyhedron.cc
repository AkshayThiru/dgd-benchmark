#include "inc/polyhedron.h"

#include <cassert>
#include <vector>

#include "inc/data_types.h"

namespace inc {

void Polyhedron::SetPolyhedron(const std::vector<Vec3>& normal,
                               const std::vector<double>& offset,
                               const std::vector<int>& graph) {
  assert(graph[0] >= 0);
  nfacet_ = graph[0];

  gamma_.resize(nfacet_);
  for (int i = 0; i < nfacet_; ++i) {
    // Hyperplane constraint: normal[i].dot(z) + offset[i] <= 0.
    gamma_[i] = -normal[i] / offset[i];
  }

  degree_.resize(nfacet_);
  for (int i = 0; i < nfacet_ - 1; ++i) {
    degree_[i] = graph[2 + i + 1] - graph[2 + i] - 1;
  }
  const int nridge = graph[1];
  degree_[nfacet_ - 1] = nfacet_ + 2 * nridge - graph[2 + nfacet_ - 1] - 1;

  graph_ = graph;
  fnid_off_ = graph_.begin() + 2;
  fnid_adr_ = fnid_off_ + nfacet_;
}

}  // namespace inc
