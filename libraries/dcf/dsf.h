#ifndef DGD_BENCHMARK_DCF_DSF_H_
#define DGD_BENCHMARK_DCF_DSF_H_

#include <cmath>
#include <cstddef>
#include <vector>

#include "dcf/precompiled.h"
#include "dcf/utils.h"

namespace dcf {

// Differentiable support function base class.
class DSF {
 public:
  DSF() {}

  virtual ~DSF() {}

  virtual void SupportFunction(const Vec3& /*x*/, const Vec3& /*pos*/,
                               const Rotation3& /*R*/, Vec3& /*s*/) {}

  virtual void SupportFunction(const Vec3& /*x*/, const Vec3& /*pos*/,
                               const Rotation3& /*R*/, Vec3& /*s*/,
                               Mat3& /*dsdx*/) {}

  virtual void SupportFunction(const Vec3& /*x*/, const Vec3& /*pos*/,
                               const Rotation3& /*R*/, Vec3& /*s*/,
                               Mat3& /*dsdx*/, Mat<3, 6>& /*dsdq*/) {}

  virtual void PrintInfo() {}
};

// Vertex DSF class.
template <int exp>
class VDSF : public DSF {
 public:
  VDSF(const std::vector<Vec3>& vert)
      : DSF(), vert_(vert), nvert_(vert_.size()) {
    static_assert(exp > 1, "Exponent must be greater than 1.");

    vvt_.resize(nvert_);
    for (std::size_t i = 0; i < nvert_; ++i) {
      vvt_[i] = vert_[i] * vert_[i].transpose();
    }
    vx_.resize(nvert_);
  }

  ~VDSF() {}

  void SupportFunction(const Vec3& x, const Vec3& pos, const Rotation3& R,
                       Vec3& s) final override {
    const Vec3 xt = R.transpose() * x;
    double max_vx = 0.0;
    for (std::size_t i = 0; i < nvert_; ++i) {
      vx_[i] = std::max(vert_[i].dot(xt), 0.0);
      max_vx = std::max(max_vx, vx_[i]);
    }

    s.setZero();
    double rvx, powa;
    const double max_vxi = 1.0 / max_vx;
    double powa_sum = 0.0;
    for (size_t i = 0; i < nvert_; ++i) {
      rvx = vx_[i] * max_vxi;
      if (rvx > kThresh) {
        powa = Power<exp - 1>(rvx);
        s += powa * vert_[i];
        powa_sum += rvx * powa;
      }
    }
    double h = std::pow(powa_sum, 1.0 / exp);
    double h2 = h / powa_sum;
    s *= h2;
    s = pos + R * s;
  }

  void SupportFunction(const Vec3& x, const Vec3& pos, const Rotation3& R,
                       Vec3& s, Mat3& dsdx) final override {
    const Vec3 xt = R.transpose() * x;
    double max_vx = 0.0;
    for (std::size_t i = 0; i < nvert_; ++i) {
      vx_[i] = std::max(vert_[i].dot(xt), 0.0);
      max_vx = std::max(max_vx, vx_[i]);
    }

    s.setZero();
    dsdx.setZero();
    double rvx, powa, powa2;
    const double max_vxi = 1.0 / max_vx;
    double powa_sum = 0.0;
    for (std::size_t i = 0; i < nvert_; ++i) {
      rvx = vx_[i] * max_vxi;
      if (rvx > kThresh) {
        powa2 = Power<exp - 2>(rvx);
        powa = powa2 * rvx;
        s += powa * vert_[i];
        dsdx += powa2 * vvt_[i];
        powa_sum += rvx * powa;
      }
    }
    double h = std::pow(powa_sum, 1.0 / exp);
    double h2 = h / powa_sum;
    s *= h2;
    dsdx = ((exp - 1) * h2 * max_vxi) * dsdx;
    dsdx -= ((exp - 1) / (h * max_vx)) * s * s.transpose();
    s = pos + R * s;
    dsdx = R * dsdx * R.transpose();
  }

  void SupportFunction(const Vec3& x, const Vec3& pos, const Rotation3& R,
                       Vec3& s, Mat3& dsdx, Mat<3, 6>& dsdq) final override {
    SupportFunction(x, pos, R, s, dsdx);
    dsdq.block<3, 3>(0, 0) = Mat3::Identity();
    Mat3 Ws_pos, Wx;
    Skew(s - pos, Ws_pos);
    Skew(x, Wx);
    dsdq.block<3, 3>(0, 3) = (-Ws_pos + dsdx * Wx) * R;
  }

  void PrintInfo() final override {
    std::cout << "===== VDSF info =====" << std::endl;
    std::cout << "nvert: " << nvert_ << std::endl;
    std::cout << "exp  : " << exp << std::endl;
    std::cout << "=====================" << std::endl;
  }

 private:
  std::vector<Mat3> vvt_;
  const std::vector<Vec3> vert_;
  std::vector<double> vx_;
  const std::size_t nvert_;

  static constexpr double kThresh = 0.1;
};

}  // namespace dcf

#endif  // DGD_BENCHMARK_DCF_DSF_H_
