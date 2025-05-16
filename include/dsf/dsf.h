#ifndef DSF_DSF_H_
#define DSF_DSF_H_

#include <cmath>
#include <vector>

#include "dsf/precompiled.h"
#include "dsf/utils.h"

namespace dsf {

// Differentiable support function base class.
class DSF {
 public:
  DSF() {}

  virtual ~DSF() {}

  virtual void SupportFunc(const Vector3d& /*x*/, const Vector3d& /*pos*/,
                           const Matrix3d& /*R*/, Vector3d& /*s*/) {}

  virtual void SupportFunc(const Vector3d& /*x*/, const Vector3d& /*pos*/,
                           const Matrix3d& /*R*/, Vector3d& /*s*/,
                           Matrix3d& /*dsdx*/) {}

  virtual void SupportFunc(const Vector3d& /*x*/, const Vector3d& /*pos*/,
                           const Matrix3d& /*R*/, Vector3d& /*s*/,
                           Matrix3d& /*dsdx*/, Matrix<double, 3, 6>& /*dsdq*/) {
  }

  virtual void PrintInfo() {}
};

// Vertex DSF class.
template <int exp>
class VDSF : public DSF {
 private:
  std::vector<Matrix3d> vvt_;
  const std::vector<Vector3d> vert_;
  std::vector<double> vx_;
  const int N_v_;

  static constexpr double thresh_{0.1};

 public:
  VDSF(const std::vector<Vector3d>& vert)
      : DSF(), vert_(vert), N_v_(static_cast<int>(vert_.size())) {
    static_assert(exp > 1, "Exponent must be greater than 1.");

    vvt_.resize(N_v_);
    for (int i = 0; i < N_v_; ++i) {
      vvt_[i] = vert_[i] * vert_[i].transpose();
    }
    vx_.resize(N_v_);
  }

  ~VDSF() {}

  void SupportFunc(const Vector3d& x, const Vector3d& pos, const Matrix3d& R,
                   Vector3d& s) final override {
    const Vector3d xt{R.transpose() * x};
    double max_vx{0.0};
    for (int i = 0; i < N_v_; ++i) {
      vx_[i] = std::max(vert_[i].dot(xt), 0.0);
      max_vx = std::max(max_vx, vx_[i]);
    }

    s.setZero();
    double rvx, powa;
    const double max_vxi{1.0 / max_vx};
    double powa_sum{0.0};
    for (int i = 0; i < N_v_; ++i) {
      rvx = vx_[i] * max_vxi;
      if (rvx > thresh_) {
        powa = Power<exp - 1>(rvx);
        s += powa * vert_[i];
        powa_sum += rvx * powa;
      }
    }
    double h{std::pow(powa_sum, 1.0 / exp)};
    double h2{h / powa_sum};
    s *= h2;
    s = pos + R * s;
  }

  void SupportFunc(const Vector3d& x, const Vector3d& pos, const Matrix3d& R,
                   Vector3d& s, Matrix3d& dsdx) final override {
    const Vector3d xt{R.transpose() * x};
    double max_vx{0.0};
    for (int i = 0; i < N_v_; ++i) {
      vx_[i] = std::max(vert_[i].dot(xt), 0.0);
      max_vx = std::max(max_vx, vx_[i]);
    }

    s.setZero();
    dsdx.setZero();
    double rvx, powa, powa2;
    const double max_vxi{1.0 / max_vx};
    double powa_sum{0.0};
    for (int i = 0; i < N_v_; ++i) {
      rvx = vx_[i] * max_vxi;
      if (rvx > thresh_) {
        powa2 = Power<exp - 2>(rvx);
        powa = powa2 * rvx;
        s += powa * vert_[i];
        dsdx += powa2 * vvt_[i];
        powa_sum += rvx * powa;
      }
    }
    double h{std::pow(powa_sum, 1.0 / exp)};
    double h2{h / powa_sum};
    s *= h2;
    dsdx = ((exp - 1) * h2 * max_vxi) * dsdx;
    dsdx -= ((exp - 1) / (h * max_vx)) * s * s.transpose();
    s = pos + R * s;
    dsdx = R * dsdx * R.transpose();
  }

  void SupportFunc(const Vector3d& x, const Vector3d& pos, const Matrix3d& R,
                   Vector3d& s, Matrix3d& dsdx,
                   Matrix<double, 3, 6>& dsdq) final override {
    SupportFunc(x, pos, R, s, dsdx);
    dsdq.block<3, 3>(0, 0) = Matrix3d::Identity();
    Matrix3d Ws_pos, Wx;
    Skew(s - pos, Ws_pos);
    Skew(x, Wx);
    dsdq.block<3, 3>(0, 3) = (-Ws_pos + dsdx * Wx) * R;
  }

  void PrintInfo() final override {
    std::cout << "===== VDSF info =====" << std::endl;
    std::cout << "N_v: " << N_v_ << std::endl;
    std::cout << "exp: " << exp << std::endl;
    std::cout << "=====================" << std::endl;
  }
};

}  // namespace dsf

#endif  // DSF_DSF_H_
