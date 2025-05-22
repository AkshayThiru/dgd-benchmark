#ifndef INC_POLYHEDRON_H_
#define INC_POLYHEDRON_H_

#include <cassert>
#include <cmath>
#include <vector>

#include "inc/data_types.h"

namespace inc {

class Polyhedron {
 public:
  explicit Polyhedron() {};

  explicit Polyhedron(const std::vector<Vec3>& normal,
                      const std::vector<double>& offset,
                      const std::vector<int>& graph);

  ~Polyhedron() = default;

  void SetPolyhedron(const std::vector<Vec3>& normal,
                     const std::vector<double>& offset,
                     const std::vector<int>& graph);

  // Checks if nfid is a neighbour of fid. If so, returns its neighbour index;
  // otherwise returns -1.
  int NeighbourIndex(int fid, int nfid) const;

  // Ray cast from the origin to the polyhedron boundary; computes a facet id
  // (fid) intersecting with the ray, and returns gamma_[fid].dot(ray).
  // The current value of fid is used as warm start.
  double RayCast(const Vec3& ray, int& fid) const;

  // Computes the largest step size s (along with facet nfid) such that:
  // gamma_[id].dot(x + s * dx) <= sigma - s,
  // for all neighbours id of fid. If skip_neighbour = true, the neighbour
  // facet nfid_ is skipped.
  template <bool skip_neighbour = false>
  double LocalFaceStepSize(const Vec3& x, const Vec3& dx, double sigma, int fid,
                           int& nfid, [[maybe_unused]] int nfid_ = -2) const;

  // Same as LocalFaceStepSize, except that only neighbours of the facet with
  // the smaller degree among fid1 and fid2 are checked.
  double LocalEdgeStepSize(const Vec3& x, const Vec3& dx, double sigma,
                           int fid1, int fid2, int& nfid) const;

  // Same as LocalFaceStepSize, except that step size feasibility is checked for
  // all facets of the polyhedron.
  double GlobalStepSize(const Vec3& x, const Vec3& dx, double sigma, int fid,
                        int& gfid) const;

  // Checks if (x, sigma) is feasible, by checking all neighbours of fid.
  // If skip_neighbour = true, the neighbour facet nfid_ is skipped.
  template <bool skip_neighbour = false>
  bool LocalFaceFeasible(const Vec3& x, double sigma, int fid,
                         [[maybe_unused]] int nfid_ = -2) const;

  // Same as LocalFaceFeasible, except that only neighbours of the facet
  // with the smaller degree among fid1 and fid2 are checked.
  bool LocalEdgeFeasible(const Vec3& x, double sigma, int fid1, int fid2) const;

  const std::vector<Vec3>& Gamma() const;

  int nfacet() const;

 private:
  // Facet normal vectors (gamma_[i].dot(z) <= 1).
  std::vector<Vec3> gamma_;
  std::vector<int> graph_, degree_;
  std::vector<int>::const_iterator fnid_off_, fnid_adr_;
  int nfacet_;
};

inline Polyhedron::Polyhedron(const std::vector<Vec3>& normal,
                              const std::vector<double>& offset,
                              const std::vector<int>& graph) {
  SetPolyhedron(normal, offset, graph);
}

inline int Polyhedron::NeighbourIndex(int fid, int nfid) const {
  assert((fid >= 0) && (fid < nfacet_));
  const auto adr = fnid_adr_ + *(fnid_off_ + fid);
  for (int id, i{0}; (id = *(adr + i)) != -1; ++i) {
    if (id == nfid) return i;
  }
  return -1;
}

inline double Polyhedron::RayCast(const Vec3& ray, int& fid) const {
  double val, max{gamma_[fid].dot(ray)};
  auto adr = fnid_adr_ + *(fnid_off_ + fid);
  int id{*adr};
  do {
    val = gamma_[id].dot(ray);
    if (val > max) {
      max = val;
      fid = id;
      adr = fnid_adr_ + *(fnid_off_ + fid);
    } else {
      ++adr;
    }
  } while ((id = *adr) != -1);
  assert(max > 0.0);
  return max;
}

template <bool skip_neighbour>
inline double Polyhedron::LocalFaceStepSize(const Vec3& x, const Vec3& dx,
                                            double sigma, int fid, int& nfid,
                                            int nfid_) const {
  assert((fid >= 0) && (fid < nfacet_));
  auto adr = fnid_adr_ + *(fnid_off_ + fid);
  Vec3 gamma;
  double den, num, val, min{kInf};
  int id{*adr};
  do {
    if constexpr (skip_neighbour) {
      if (id == nfid_) continue;
    }
    gamma = gamma_[id];
    den = 1.0 + gamma.dot(dx);
    if (den > 0.0) {
      num = std::max(sigma - gamma.dot(x), 0.0);  // For robustness.
      // if (num > 0.0) { // Not necessary.
      val = num / den;
      if (val < min) {
        min = val;
        nfid = id;
      }
      // }
    }
  } while ((id = *(++adr)) != -1);
  assert(min >= 0.0);
  return min;
}

inline double Polyhedron::LocalEdgeStepSize(const Vec3& x, const Vec3& dx,
                                            double sigma, int fid1, int fid2,
                                            int& nfid) const {
  if (degree_[fid1] < degree_[fid2]) {
    return LocalFaceStepSize<true>(x, dx, sigma, fid1, nfid, fid2);
  } else {
    return LocalFaceStepSize<true>(x, dx, sigma, fid2, nfid, fid1);
  }
}

inline double Polyhedron::GlobalStepSize(const Vec3& x, const Vec3& dx,
                                         double sigma, int fid,
                                         int& gfid) const {
  assert((fid >= 0) && (fid < nfacet_));
  auto adr = fnid_adr_ + *(fnid_off_ + fid);
  Vec3 gamma;
  double den, num, val, min{kInf};
  int id{*adr};
  do {
    // if (id == fid) {  // Not required.
    //   ++adr;
    //   continue;
    // }
    gamma = gamma_[id];
    den = 1.0 + gamma.dot(dx);
    if (den > 0.0) {
      num = std::max(sigma - gamma.dot(x), 0.0);  // For robustness.
      // if (num > 0.0) { // Not necessary.
      val = num / den;
      if (val < min) {
        min = val;
        gfid = id;
        adr = fnid_adr_ + *(fnid_off_ + gfid);
        continue;
      }
      // }
    }
    ++adr;
  } while ((id = *adr) != -1);
  assert(min >= 0.0);
  return min;
}

template <bool skip_neighbour = false>
inline bool Polyhedron::LocalFaceFeasible(const Vec3& x, double sigma, int fid,
                                          int nfid_) const {
  assert((fid >= 0) && (fid < nfacet_));
  auto adr = fnid_adr_ + *(fnid_off_ + fid);
  int id{*adr};
  do {
    if constexpr (skip_neighbour) {
      if (id == nfid_) continue;
    }
    if (gamma_[id].dot(x) > sigma) {
      return false;
    }
  } while ((id = *(++adr)) != -1);
  return true;
}

inline bool Polyhedron::LocalEdgeFeasible(const Vec3& x, double sigma, int fid1,
                                          int fid2) const {
  if (degree_[fid1] < degree_[fid2]) {
    return LocalFaceFeasible<true>(x, sigma, fid1, fid2);
  } else {
    return LocalFaceFeasible<true>(x, sigma, fid2, fid1);
  }
}

inline const std::vector<Vec3>& Polyhedron::Gamma() const { return gamma_; }

inline int Polyhedron::nfacet() const { return nfacet_; }

}  // namespace inc

#endif  // INC_POLYHEDRON_H_
