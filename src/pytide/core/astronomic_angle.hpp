// Copyright (c) 2020 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
#pragma once
#include <limits>
#include "math.hpp"

/// Astronomical angles.
class AstronomicAngle {
 protected:
  /// Hour angle of mean sun.
  double t_;
  /// Longitude of moon's node
  double n_;
  /// Mean longitude of the sun.
  double h_;
  /// Mean longitude of the moon.
  double s_;
  /// Mean longitude of solar perigee.
  double p1_;
  /// Mean longitude of lunar perigee.
  double p_;
  /// Obliquity of lunar orbit with respect to earth's equator
  double i_;
  /// Longitude in moon's orbit of lunar intersection
  double xi_;
  /// Right ascension of lunar intersection
  double nu_;
  /// Factor in amplitude of constituent L₂
  double x1ra_;
  /// Term in argument of constituent L₂
  double r_;
  /// Term in argument of lunisolar constituent K₁
  double nuprim_;
  /// Term in argument of lunisolar constituent K₂
  double nusec_;

 public:
  /// Default constructor.
  AstronomicAngle() {
    t_ = n_ = h_ = s_ = p1_ = p_ = i_ = xi_ = nu_ = x1ra_ = r_ = nuprim_ =
        nusec_ = std::numeric_limits<double>::quiet_NaN();
  }

  /// Initialize some astronomic data useful for nodal corrections.
  ///
  /// @param date Desired UTC time
  explicit AstronomicAngle(const double epoch);

  /// Returns the hour angle of mean sun.
  double t() const noexcept { return t_; }

  /// Returns the longitude of moon's node (radians)
  double n() const noexcept { return n_; }

  /// Mean longitude of the sun (radians)
  double h() const noexcept { return h_; }

  /// Mean longitude of the moon (radians)
  double s() const noexcept { return s_; }

  /// Mean longitude of solar perigee (radians)
  double p1() const noexcept { return p1_; }

  /// Mean longitude of lunar perigee (radians)
  double p() const noexcept { return p_; }

  /// Obliquity of lunar orbit with respect to earth's equator (radians)
  double i() const noexcept { return i_; }

  /// Longitude in moon's orbit of lunar intersection (radians)
  double xi() const noexcept { return xi_; }

  /// Right ascension of lunar intersection (radians)
  double nu() const noexcept { return nu_; }

  /// Factor in amplitude of constituent L₂
  double x1ra() const noexcept { return x1ra_; }

  /// Term in argument of constituent L₂ (radians)
  double r() const noexcept { return r_; }

  /// Term in argument of lunisolar constituent K₁ (radians)
  double nuprim() const noexcept { return nuprim_; }

  /// Term in argument of lunisolar constituent K₂ (radians)
  double nusec() const noexcept { return nusec_; }

  /// f_o1  = sin I cos² ½I / 0.3800
  double f_o1() const {
    return std::sin(i_) * sqr(std::cos(i_ * 0.5)) / 0.3800;
  }

  /// f_oo1  = sin I sin² ½I / 0.0164
  double f_oo1() const {
    return std::sin(i_) * sqr(std::sin(i_ * 0.5)) / 0.01640;
  }

  /// f_1  = 1
  constexpr double f_1() const noexcept { return 1; }

  /// f_j1  = sin 2I / 0.7214
  double f_j1() const { return std::sin(2.0 * i_) / 0.7214; }

  /// f_m13 = (1 -10 sin² ½I + 15 sin⁴ ½I) cos² ½I / 0.5873
  double f_m13() const {
    return f_o1() * std::sqrt(2.310 + 1.435 * std::cos(2.0 * (p_ - xi_)));
  }

  /// f_m2 = cos⁴ ½I / 0.9154
  double f_m2() const { return pow4(std::cos(i_ * 0.5)) / 0.9154; }

  /// f_m3 = cos⁶ ½I / 0.8758
  double f_m3() const { return std::pow(std::cos(i_ * 0.5), 6.0) / 0.8758; }

  /// f_mf = sin² I / 0.1578
  double f_mf() const { return sqr(std::sin(i_)) / 0.1578; }

  /// f_mm = (2/3 - sin² I) / 0.5021
  double f_mm() const { return (2.0 / 3.0 - sqr(std::sin(i_))) / 0.5021; }

  /// f_m22 = f²(M₂)
  double f_m22() const { return sqr(f_m2()); }

  /// f_m23 = f(M₂)³
  double f_m23() const { return pow3(f_m2()); }

  /// f_m24 = f(M₂)⁴
  double f_m24() const { return pow4(f_m2()); }

  /// f_k1 = √(0.8965 sin² 2I+0.6001 sin 2I cos ν + 0.1006)
  double f_k1() const {
    return std::sqrt(0.8965 * sqr(std::sin(2.0 * i_)) +
                     0.6001 * sin(2.0 * i_) * std::cos(nu_) + 0.1006);
  }

  /// f_k2 = √(19.0444 sin⁴ I + 2.7702 sin² I cos 2ν + 0.0981)
  double f_k2() const {
    return sqrt(19.0444 * pow4(std::sin(i_)) +
                2.7702 * sqr(std::sin(i_)) * std::cos(2.0 * nu_) + 0.0981);
  }

  /// f_kj2 = sin² I / 0.1565 (formula #79)
  double f_kj2() const { return sqr(std::sin(i_)) / 0.1565; }

  /// f_l2 = _f_m2() * 1 / Ra
  double f_l2() const { return f_m2() * x1ra_; }

  /// f_m2_k2 = f(M₂) * f(K₂)
  double f_m2_k2() const { return f_m2() * f_k2(); }

  /// f_m2_k1 = f(M₂) * f(K₁)
  double f_m2_k1() const { return f_m2() * f_k1(); }

  /// f_m2_o1 = f(M₂) * f(O₁)
  double f_m2_o1() const { return f_m2() * f_o1(); }

  /// f_m22_k1 = f(M₂)² * f(K₁)
  double f_m22_k1() const { return f_m22() * f_k1(); }

  /// f_m23_k2 = f(M₂)³ * f(K₂)
  double f_m23_k2() const { return f_m23() * f_k2(); }
};
