#ifndef DGD_BENCHMARK_DSF_DSF_INTERFACE_H_
#define DGD_BENCHMARK_DSF_DSF_INTERFACE_H_

#include <memory>
#include <string>
#include <vector>

#include "dgd/geometry/convex_set.h"
#include "dsf/dsf.h"

namespace dsf {

// DSF interface class.
template <int exp>
class VDSFInterface : public dgd::ConvexSet<3> {
 public:
  explicit VDSFInterface(const std::vector<Vec3>& vert, double inradius,
                         double margin);

  ~VDSFInterface() {};

  double SupportFunction(
      const Vec3& n, Vec3& sp,
      dgd::SupportFunctionHint<3>* /*hint*/ = nullptr) const final override;

  bool RequireUnitNormal() const final override;

  VDSF<exp>* VDSFPtr() const;

 private:
  std::unique_ptr<VDSF<exp>> vdsf_;
  const double margin_;
  const double inradius_;
};

template <int exp>
inline VDSFInterface<exp>::VDSFInterface(const std::vector<Vec3>& vert,
                                         double inradius, double margin)
    : margin_(margin), inradius_(inradius) {
  vdsf_ = std::make_unique<VDSF<exp>>(vert);
}

template <int exp>
inline double VDSFInterface<exp>::SupportFunction(
    const Vec3& n, Vec3& sp, dgd::SupportFunctionHint<3>* /*hint*/) const {
  vdsf_->SupportFunction(n, Vec3::Zero(), Rotation3::Identity(), sp);
  sp += margin_ * n;
  return n.dot(sp);
}

template <int exp>
inline bool VDSFInterface<exp>::RequireUnitNormal() const {
  return (margin_ > 0.0);
}

template <int exp>
inline VDSF<exp>* VDSFInterface<exp>::VDSFPtr() const {
  return vdsf_.get();
}

}  // namespace dsf

#endif  // DGD_BENCHMARK_DSF_DSF_INTERFACE_H_
