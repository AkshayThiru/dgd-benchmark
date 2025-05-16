#ifndef DSF_DSF_INTERFACE_H_
#define DSF_DSF_INTERFACE_H_

#include <dgd/data_types.h>
#include <dgd/geometry/convex_set.h>

#include <memory>
#include <string>
#include <vector>

#include "dsf/dsf.h"

namespace dsf {

// Load mesh object from file and return vertices and inradius.
double LoadOBJ(const std::string& file, std::vector<dgd::Vec3f>& vert);

// DSF interface class.
template <int exp>
class VDSFInterface : public dgd::ConvexSet<3> {
 public:
  explicit VDSFInterface(const std::vector<dgd::Vec3f>& vert, double inradius,
                         double margin);

  ~VDSFInterface() {};

  double SupportFunction(
      const dgd::Vec3f& n, dgd::Vec3f& sp,
      dgd::SupportFunctionHint<3>* /*hint*/ = nullptr) const final override;

  bool RequireUnitNormal() const final override;

  VDSF<exp>* VDSFPtr() const;

 private:
  std::unique_ptr<VDSF<exp>> vdsf_;
  const double margin_;
  const double inradius_;
};

template <int exp>
inline VDSFInterface<exp>::VDSFInterface(const std::vector<dgd::Vec3f>& vert,
                                         double inradius, double margin)
    : margin_(margin), inradius_(inradius) {
  vdsf_ = std::make_unique<VDSF<exp>>(vert);
}

template <int exp>
inline double VDSFInterface<exp>::SupportFunction(
    const dgd::Vec3f& n, dgd::Vec3f& sp,
    dgd::SupportFunctionHint<3>* /*hint*/) const {
  vdsf_->SupportFunc(n, dgd::Vec3f::Zero(), dgd::Rot3f::Identity(), sp);
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

#endif  // DSF_DSF_INTERFACE_H_
