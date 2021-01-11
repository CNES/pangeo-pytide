// Copyright (c) 2020 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
#pragma once
#include <array>
#include <complex>
#include <memory>
#include <string>
#include <vector>
#include <Eigen/Core>
#include "angle.hpp"
#include "astronomic_angle.hpp"
#include "math.hpp"

/// @brief Wave definition
class Wave : public std::enable_shared_from_this<Wave> {
 public:
  /// Typename to a function pointer for calculate the nodal factor
  using NodalFactor = double (AstronomicAngle::*)() const;

  /// @brief Possible type of tidal wave.
  enum TidalType {
    kLongPeriod = 0,  //!< Long period tidal waves
    kShortPeriod      //!< Short period tidal waves
  };

  /// @brief Index to access the wave in the internal table
  enum Ident : size_t {
    kMm = 0,        //!< %Mm
    kMf = 1,        //!< %Mf
    kMtm = 2,       //!< %Mtm
    kMsqm = 3,      //!< %Msqm
    k2Q1 = 4,       //!< 2Q₁
    kSigma1 = 5,    //!< σ₁
    kQ1 = 6,        //!< Q₁
    kRho1 = 7,      //!< ρ₁
    kO1 = 8,        //!< O₁
    kMP1 = 9,       //!< MP₁
    kM11 = 10,      //!< M₁₁
    kM12 = 11,      //!< M₁₂
    kM13 = 12,      //!< M₁₃
    kChi1 = 13,     //!< χ₁
    kPi1 = 14,      //!< π₁
    kP1 = 15,       //!< P₁
    kS1 = 16,       //!< S₁
    kK1 = 17,       //!< K₁
    kPsi1 = 18,     //!< ψ₁
    kPhi1 = 19,     //!< φ₁
    kTheta1 = 20,   //!< θ₁
    kJ1 = 21,       //!< J₁
    kOO1 = 22,      //!< OO₁
    kMNS2 = 23,     //!< MNS₂
    kEps2 = 24,     //!< ε₂
    k2N2 = 25,      //!< 2N₂
    kMu2 = 26,      //!< µ₂
    k2MS2 = 27,     //!< 2MS₂
    kN2 = 28,       //!< N₂
    kNu2 = 29,      //!< ν₂
    kM2 = 30,       //!< M₂
    kMKS2 = 31,     //!< MKS₂
    kLambda2 = 32,  //!< λ₂
    kL2 = 33,       //!< L₂
    k2MN2 = 34,     //!< 2MN₂
    kT2 = 35,       //!< T₂
    kS2 = 36,       //!< S₂
    kR2 = 37,       //!< R₂
    kK2 = 38,       //!< K₂
    kMSN2 = 39,     //!< MSN₂
    kEta2 = 40,     //!< η₂
    k2SM2 = 41,     //!< 2SM₂
    kMO3 = 42,      //!< MO₃
    k2MK3 = 43,     //!< 2MK₃
    kM3 = 44,       //!< M₃
    kMK3 = 45,      //!< MK₃
    kN4 = 46,       //!< N₄
    kMN4 = 47,      //!< MN₄
    kM4 = 48,       //!< M₄
    kSN4 = 49,      //!< SN₄
    kMS4 = 50,      //!< MS₄
    kMK4 = 51,      //!< MK₄
    kS4 = 52,       //!< S₄
    kSK4 = 53,      //!< SK₄
    kR4 = 54,       //!< R₄
    k2MN6 = 55,     //!< 2MN₆
    kM6 = 56,       //!< M₆
    kMSN6 = 57,     //!< MSN₆
    k2MS6 = 58,     //!< 2MS₆
    k2MK6 = 59,     //!< 2MK₆
    k2SM6 = 60,     //!< 2SM₆
    kMSK6 = 61,     //!< MSK₆
    kS6 = 62,       //!< S₆
    kM8 = 63,       //!< %M8
    kMSf = 64,      //!< %MSf
    kSsa = 65,      //!< %Ssa
    kSa = 66,       //!< %Sa
  };

 protected:
  /// nodal correction for phase
  double u_{std::numeric_limits<double>::quiet_NaN()};

 private:
  /// Wave ident
  Ident ident_;

  /// Type of tide.
  TidalType type_;

  /// Function to call for computing the node factor
  NodalFactor calculate_node_factor_;

  /// Wave frequency.
  double freq_;

  /// greenwich argument
  double v_{std::numeric_limits<double>::quiet_NaN()};

  /// Nodal correction for amplitude.
  double f_{std::numeric_limits<double>::quiet_NaN()};

  /// Harmonic constituents (T, s, h, p, N′, p₁, shift, ξ, ν, ν′, ν″)
  std::array<int16_t, 11> argument_;

  /// Computes the wave frequency from the doodson arguments
  ///
  /// @param t Mean solar angle relative to Greenwich
  /// @param s moon's mean longitude
  /// @param h sun's mean longitude
  /// @param p longitude of moon's perigee
  /// @param n longitude of moon's ascending node
  /// @param p1 longitude of sun's perigee
  static constexpr double frequency(const int16_t t, const int16_t s,
                                    const int16_t h, const int16_t p,
                                    const int16_t n, const int16_t p1) {
    return ((frequency::tau() + frequency::s() - frequency::h()) * t +
            frequency::s() * s + frequency::h() * h + frequency::p() * p +
            frequency::n() * n + frequency::p1() * p1) *
           360;
  }

 public:
  /// Initializes the properties of the wave (frequency, doodson's coefficients,
  /// etc.).
  ///
  /// @param ident Index of the wave in the internal table
  /// @param t Mean solar angle relative to Greenwich
  /// @param s moon's mean longitude
  /// @param h sun's mean longitude
  /// @param p longitude of moon's perigee
  /// @param n longitude of moon's ascending node
  /// @param p1 longitude of sun's perigee
  /// @param shift TODO
  /// @param eps Coefficient for the longitude in moon's orbit of lunar
  ///   intersection
  /// @param nu Coefficient for the right ascension of lunar intersection
  /// @param nuprim Coefficient for the term in argument of lunisolar
  ///   constituent K₁
  /// @param nusec Coefficient for the term in argument of lunisolar constituent
  ///   K₂
  /// @param type Type of tidal wave
  /// @param calculate_node_factor Function used to calculate the nodal factor
  Wave(const Ident ident, const int16_t t, const int16_t s, const int16_t h,
       const int16_t p, const int16_t n, const int16_t p1, const int16_t shift,
       const int16_t eps, const int16_t nu, const int16_t nuprim,
       const int16_t nusec, TidalType type, NodalFactor calculate_node_factor)
      : ident_(ident),
        type_(type),
        calculate_node_factor_(calculate_node_factor),
        freq_(radians(frequency(t, s, h, p, n, p1)) / 3600.0) {
    argument_[0] = t;
    argument_[1] = s;
    argument_[2] = h;
    argument_[3] = p;
    argument_[4] = n;
    argument_[5] = p1;
    argument_[6] = shift;
    argument_[7] = eps;
    argument_[8] = nu;
    argument_[9] = nuprim;
    argument_[10] = nusec;
  }

  /// Default destructor
  virtual ~Wave() = default;

  /// Default copy constructor
  Wave(const Wave&) = default;

  /// Default copy assignment operator
  Wave& operator=(const Wave&) = default;

  /// Move constructor
  Wave(Wave&&) noexcept = default;

  /// Move assignment operator
  Wave& operator=(Wave&&) noexcept = default;

  /// Compute nodal corrections from SCHUREMAN (1958).
  ///
  /// @param a Astronomic angle
  void nodal_a(const AstronomicAngle& a) { f_ = (a.*calculate_node_factor_)(); }

  /// Compute nodal corrections from SCHUREMAN (1958).
  ///
  /// @param a Astronomic angle
  virtual void nodal_g(const AstronomicAngle& a);

  /// Gets the wave ident
  constexpr Ident ident() const noexcept { return ident_; }

  /// Gets the wave frequency (radians per seconds)
  constexpr double freq() const noexcept { return freq_; }

  /// Gets the wave type
  constexpr TidalType type() const noexcept { return type_; }

  /// Gets v (greenwich argument) + u (nodal correction for phase)
  double vu() const noexcept { return std::fmod(v_ + u_, two_pi<double>()); }

  /// Gets v0 (greenwich argument)
  double v() const noexcept { return v_; }

  /// Gets the nodal correction for amplitude
  constexpr double f() const noexcept { return f_; }

  /// Gets the nodal correction for phase
  constexpr double u() const noexcept { return u_; }

  /// Gets the wave name
  std::string name() const;
};

/// Mm
///
/// V = s - p;
/// u = 0;
/// f = f(Mm)
class Mm : public Wave {
 public:
  Mm()
      : Wave(kMm, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, kLongPeriod,
             &AstronomicAngle::f_mm) {}
};

/// Mf
///
/// V = 2s
/// u = -2ξ
/// f = f(Mf)
class Mf : public Wave {
 public:
  Mf()
      : Wave(kMf, 0, 2, 0, 0, 0, 0, 0, -2, 0, 0, 0, kLongPeriod,
             &AstronomicAngle::f_mf) {}
};

/// Mtm
///
/// V = 3s - p
/// u = -2ξ
/// f = f(Mf)
class Mtm : public Wave {
 public:
  Mtm()
      : Wave(kMtm, 0, 3, 0, -1, 0, 0, 0, -2, 0, 0, 0, kLongPeriod,
             &AstronomicAngle::f_mf) {}
};

/// Msqm
///
/// V = 4s - 2h
/// u = -2ξ
/// f = f(Mf)
class Msqm : public Wave {
 public:
  Msqm()
      : Wave(kMsqm, 0, 4, -2, 0, 0, 0, 0, -2, 0, 0, 0, kLongPeriod,
             &AstronomicAngle::f_mf) {}
};

/// Ssa
///
/// V = 2h
/// u = 0
/// f = 1
class Ssa : public Wave {
 public:
  Ssa()
      : Wave(kSsa, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, kLongPeriod,
             &AstronomicAngle::f_1) {}
};

/// Sa
///
/// V = h
/// u = 0
/// f = 1
class Sa : public Wave {
 public:
  Sa()
      : Wave(kSa, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, kLongPeriod,
             &AstronomicAngle::f_1) {}
};

/// 2Q₁
///
/// V = T - 4s + h + 2p + 90°
/// u = +2ξ - ν
/// f = f(O₁)
class _2Q1 : public Wave {
 public:
  _2Q1()
      : Wave(k2Q1, 1, -4, 1, 2, 0, 0, 1, 2, -1, 0, 0, kShortPeriod,
             &AstronomicAngle::f_o1) {}
};

/// σ₁
///
/// V = T - 4s + 3h + 90°
/// u = +2ξ - ν
/// f = f(O₁)
class Sigma1 : public Wave {
 public:
  Sigma1()
      : Wave(kSigma1, 1, -4, 3, 0, 0, 0, 1, 2, -1, 0, 0, kShortPeriod,
             &AstronomicAngle::f_o1) {}
};

/// Q₁
///
/// V = T - 3s + h + p + 90°
/// u = +2ξ - ν
/// f = f(O₁)
class Q1 : public Wave {
 public:
  Q1()
      : Wave(kQ1, 1, -3, 1, 1, 0, 0, 1, 2, -1, 0, 0, kShortPeriod,
             &AstronomicAngle::f_o1) {}
};

/// ρ₁
///
/// V = T - 3s + 3h - p + 90°
/// u = +2ξ - ν
/// f = f(O₁)
class Rho1 : public Wave {
 public:
  Rho1()
      : Wave(kRho1, 1, -3, 3, -1, 0, 0, 1, 2, -1, 0, 0, kShortPeriod,
             &AstronomicAngle::f_o1) {}
};

// O₁
///
/// V = T - 2s + h + 90°
/// u = +2ξ - ν
/// f = f(O₁)
///
class O1 : public Wave {
 public:
  O1()
      : Wave(kO1, 1, -2, 1, 0, 0, 0, 1, 2, -1, 0, 0, kShortPeriod,
             &AstronomicAngle::f_o1) {}
};

/// MP₁
///
/// V = T - 2s + 3h - 90°
/// u = -ν
/// f = f(J₁)
class MP1 : public Wave {
 public:
  MP1()
      : Wave(kMP1, 1, -2, 3, 0, 0, 0, -1, 0, -1, 0, 0, kShortPeriod,
             &AstronomicAngle::f_j1) {}
};

/// M₁₂ (Formula A16)
///
/// V = T - s + h - p - 90°
/// u = +2ξ - ν
/// f = f(O₁)
class M12 : public Wave {
 public:
  M12()
      : Wave(kM12, 1, -1, 1, -1, 0, 0, -1, 2, -1, 0, 0, kShortPeriod,
             &AstronomicAngle::f_o1) {}
};

/// M₁₃ (= M11 + M12)
///
/// V = T - s + h + p - 90
/// u = -ν
/// f = f(M₁₃)
class M13 : public Wave {
 public:
  M13()
      : Wave(kM13, 1, -1, 1, 1, 0, 0, -1, 0, -1, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m13) {}

  /// Compute nodal corrections from SCHUREMAN (1958).
  ///
  /// @param a Astronomic angle
  void nodal_g(const AstronomicAngle& a) final {
    Wave::nodal_g(a);
    u_ -= radians(1.0 /
                  std::sqrt(2.310 + 1.435 * std::cos(2 * (a.p() - a.xi()))));
  }
};

/// M₁₁ (Formula A23)
///
/// V = T - s + h + p - 90°
/// u = -ν
/// f = f(J₁)
class M11 : public Wave {
 public:
  M11()
      : Wave(kM11, 1, -1, 1, 1, 0, 0, -1, 0, -1, 0, 0, kShortPeriod,
             &AstronomicAngle::f_j1) {}
};

/// χ₁
///
/// V = T - s + 3h - p - 90°
/// u = -ν
/// f = f(J₁)
class Chi1 : public Wave {
 public:
  Chi1()
      : Wave(kChi1, 1, -1, 3, -1, 0, 0, -1, 0, -1, 0, 0, kShortPeriod,
             &AstronomicAngle::f_j1) {}
};

/// π₁
///
/// V = T - 2h + p1 + 90°
/// u = 0
/// f = 1
class Pi1 : public Wave {
 public:
  Pi1()
      : Wave(kPi1, 1, 0, -2, 0, 0, 1, 1, 0, 0, 0, 0, kShortPeriod,
             &AstronomicAngle::f_1) {}
};

/// P₁
///
/// V = T - h + 90°
/// u = 0
/// f = 1
class P1 : public Wave {
 public:
  P1()
      : Wave(kP1, 1, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, kShortPeriod,
             &AstronomicAngle::f_1) {}
};

/// S₁
///
/// V = T
/// u = 0
/// f = 1
class S1 : public Wave {
 public:
  S1()
      : Wave(kS1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, kShortPeriod,
             &AstronomicAngle::f_1) {}
};

/// K₁
///
/// V = T + h - 90°
/// u = - ν'
/// f = f(k₁)
class K1 : public Wave {
 public:
  K1()
      : Wave(kK1, 1, 0, 1, 0, 0, 0, -1, 0, 0, -1, 0, kShortPeriod,
             &AstronomicAngle::f_k1) {}
};

/// ψ₁
///
/// V = T + 2h - p1 - 90°
/// u = 0
/// f = 1
class Psi1 : public Wave {
 public:
  Psi1()
      : Wave(kPsi1, 1, 0, 2, 0, 0, -1, -1, 0, 0, 0, 0, kShortPeriod,
             &AstronomicAngle::f_1) {}
};

/// φ₁
///
/// V = T + 3h - 90°
/// u = 0
/// f = 1
class Phi1 : public Wave {
 public:
  Phi1()
      : Wave(kPhi1, 1, 0, 3, 0, 0, 0, -1, 0, 0, 0, 0, kShortPeriod,
             &AstronomicAngle::f_1) {}
};

/// θ₁
///
/// V = T + s - h + p - 90°
/// u = -ν
/// f = f(J₁)
class Theta1 : public Wave {
 public:
  Theta1()
      : Wave(kTheta1, 1, 1, -1, 1, 0, 0, -1, 0, -1, 0, 0, kShortPeriod,
             &AstronomicAngle::f_j1) {}
};

/// J₁
///
/// V = T + s + h - p - 90°
/// u = -ν
/// f = f(J₁)
class J1 : public Wave {
 public:
  J1()
      : Wave(kJ1, 1, 1, 1, -1, 0, 0, -1, 0, -1, 0, 0, kShortPeriod,
             &AstronomicAngle::f_j1) {}
};

/// OO₁
///
/// V = T + 2s + h - 90°
/// u = -2ξ - ν
/// f = f(OO₁)
class OO1 : public Wave {
 public:
  OO1()
      : Wave(kOO1, 1, 2, 1, 0, 0, 0, -1, -2, -1, 0, 0, kShortPeriod,
             &AstronomicAngle::f_oo1) {}
};

/// MNS₂ = M₂ + N₂ + S₂
///
/// V = 2T - 5s + 4h + p
/// u = +4ξ - 4ν
/// f = f(M₂)²
class MNS2 : public Wave {
 public:
  MNS2()
      : Wave(kMNS2, 2, -5, 4, 1, 0, 0, 0, 4, -4, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m22) {}
};

/// ε₂
///
/// V = 2T - 5s + 4h + p
/// u = +2ξ - 2ν
/// f = f(M₂)
class Eps2 : public Wave {
 public:
  Eps2()
      : Wave(kEps2, 2, -5, 4, 1, 0, 0, 0, 2, -2, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m2) {}
};

/// 2N₂
///
/// V = 2T - 4s + 2h + 2p
/// u = +2ξ - 2ν
/// f = f(M₂)
class _2N2 : public Wave {
 public:
  _2N2()
      : Wave(k2N2, 2, -4, 2, 2, 0, 0, 0, 2, -2, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m2) {}
};

/// µ₂
///
/// V = 2T - 4s + 4h
/// u = +2ξ - 2ν
/// f = f(M₂)
class Mu2 : public Wave {
 public:
  Mu2()
      : Wave(kMu2, 2, -4, 4, 0, 0, 0, 0, 2, -2, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m2) {}
};

/// 2MS₂ = 2M₂ - S₂
///
/// V = 2T - 4s + 4h
/// u = +4ξ - 4ν
/// f = f(M₂)²
class _2MS2 : public Wave {
 public:
  _2MS2()
      : Wave(k2MS2, 2, -4, 4, 0, 0, 0, 0, 4, -4, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m22) {}
};

/// N₂
///
/// V = 2T - 3s + 2h + p
/// u = +2ξ - 2ν
/// f = f(M₂)
class N2 : public Wave {
 public:
  N2()
      : Wave(kN2, 2, -3, 2, 1, 0, 0, 0, 2, -2, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m2) {}
};

/// ν₂
///
/// V = 2T - 3s + 4h - p
/// u = +2ξ - 2ν
/// f = f(M₂)
class Nu2 : public Wave {
 public:
  Nu2()
      : Wave(kNu2, 2, -3, 4, -1, 0, 0, 0, 2, -2, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m2) {}
};

/// M₂
///
/// V = 2T - 2s + 2h
/// u = +2ξ - 2ν
/// f = f(M₂)
class M2 : public Wave {
 public:
  M2()
      : Wave(kM2, 2, -2, 2, 0, 0, 0, 0, 2, -2, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m2) {}
};

/// MKS₂ = M₂ + K₂ - S₂
///
/// V = 2T - 2s + 4h
/// u = +2ξ - 2ν -2ν''
/// f = f(M₂) × f(K₂)
class MKS2 : public Wave {
 public:
  MKS2()
      : Wave(kMKS2, 2, -2, 4, 0, 0, 0, 0, 2, -2, 0, -2, kShortPeriod,
             &AstronomicAngle::f_m2_k2) {}
};

/// λ₂
///
/// V = 2T - s + p + 180°
/// u = +2ξ - 2ν
/// f = f(M₂)
class Lambda2 : public Wave {
 public:
  Lambda2()
      : Wave(kLambda2, 2, -1, 0, 1, 0, 0, 2, 2, -2, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m2) {}
};

/// L₂
///
/// V = 2T - s + 2h - p + 180°
/// u = +2ξ - 2ν - R
/// f = f(L₂)
class L2 : public Wave {
 public:
  L2()
      : Wave(kL2, 2, -1, 2, -1, 0, 0, 2, 2, -2, 0, 0, kShortPeriod,
             &AstronomicAngle::f_l2) {}

  /// Compute nodal corrections from SCHUREMAN (1958).
  ///
  /// @param a Astronomic angle
  void nodal_g(const AstronomicAngle& a) final {
    Wave::nodal_g(a);
    u_ -= a.r();
  }
};

/// 2MN₂ = 2M₂ - N₂
///
/// V = 2T - s + 2h - p + 180°
/// u = +2ξ - 2ν
/// f = f(M₂)³
class _2MN2 : public Wave {
 public:
  _2MN2()
      : Wave(k2MN2, 2, -1, 2, -1, 0, 0, 2, 2, -2, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m23) {}
};

/// T₂
///
/// V = 2T - h + p₁
/// u = 0
/// f = 1
class T2 : public Wave {
 public:
  T2()
      : Wave(kT2, 2, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, kShortPeriod,
             &AstronomicAngle::f_1) {}
};

/// S₂
///
/// V = 2T
/// u = 0
/// f = 1
class S2 : public Wave {
 public:
  S2()
      : Wave(kS2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, kShortPeriod,
             &AstronomicAngle::f_1) {}
};

/// R₂
///
/// V = 2T + h - p1 + 180°
/// u = 0
/// f = 1
class R2 : public Wave {
 public:
  R2()
      : Wave(kR2, 2, 0, 1, 0, 0, -1, 2, 0, 0, 0, 0, kShortPeriod,
             &AstronomicAngle::f_1) {}
};

/// K₂
///
/// V = 2T + 2h
/// u = -2ν″
/// f = f(K₂)
class K2 : public Wave {
 public:
  K2()
      : Wave(kK2, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, -2, kShortPeriod,
             &AstronomicAngle::f_k2) {}
};

/// MSN₂ = M2 + S2 - N2
///
/// V = 2T + s -p
/// u = 0
/// f = f(M₂)²
class MSN2 : public Wave {
 public:
  MSN2()
      : Wave(kMSN2, 2, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m22) {}
};

/// η₂ = KJ₂
///
/// V = 2T + s + 2h - p
/// u = -2ν
/// f = f(KJ₂)
class Eta2 : public Wave {
 public:
  Eta2()
      : Wave(kEta2, 2, 1, 2, -1, 0, 0, 0, 0, -2, 0, 0, kShortPeriod,
             &AstronomicAngle::f_kj2) {}
};

/// 2SM₂ = 2S₂ - M₂
///
/// V = 2T + 2s - 2h
/// u = -2ξ + 2ν
/// f = f(M₂)
class _2SM2 : public Wave {
 public:
  _2SM2()
      : Wave(k2SM2, 2, 2, -2, 0, 0, 0, 0, 2, -2, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m2) {}
};

/// MO₃ = M₂ + O₁
///
/// V = 3T - 4s + 3h + 90°
/// u = 4ξ - 3ν
/// f = f(M₂) × f(O₁)
class MO3 : public Wave {
 public:
  MO3()
      : Wave(kMO3, 3, -4, 3, 0, 0, 0, 1, 4, -3, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m2_o1) {}
};

/// 2MK₃ = 2M₂ - K₁
///
/// V = 3T - 4s + 3h + 90°
/// u = 4ξ - 4ν + ν′
/// f = f(M₂)² × f(K₁)
class _2MK3 : public Wave {
 public:
  _2MK3()
      : Wave(k2MK3, 3, -4, 3, 0, 0, 0, 1, 4, -4, 1, 0, kShortPeriod,
             &AstronomicAngle::f_m22_k1) {}
};

/// M₃
///
/// V = 3T - 3s + 3h
/// u = +3ξ - 3ν
/// f = f(M₃)
class M3 : public Wave {
 public:
  M3()
      : Wave(kM3, 3, -3, 3, 0, 0, 0, 0, 3, -3, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m3) {}
};

/// MK₃ = M₂ + K₁
///
/// V = 3T - 2s + 3h - 90°
/// u = 2ξ - 2ν - ν′
/// f = f(M₂) × f(K₁)
class MK3 : public Wave {
 public:
  MK3()
      : Wave(kMK3, 3, -2, 3, 0, 0, 0, -1, 2, -2, -1, 0, kShortPeriod,
             &AstronomicAngle::f_m2_k1) {}
};

/// N4 = N₂ + N₂
///
/// V = 4T - 6s + 4h + 2p
/// u = +4ξ - 4ν
/// f = f(M₂)²
class N4 : public Wave {
 public:
  N4()
      : Wave(kN4, 4, -6, 4, 2, 0, 0, 0, 4, -4, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m22) {}
};

/// MN₄ = M₂ + N₂
///
/// V = 4T - 5s + 4h + p
/// u = +4ξ - 4ν
/// f = f(M₂)²
class MN4 : public Wave {
 public:
  MN4()
      : Wave(kMN4, 4, -5, 4, 1, 0, 0, 0, 4, -4, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m22) {}
};

/// M₄ = 2M₂
///
/// V = 4T - 4s + 4h
/// u = +4ξ - 4ν
/// f = f²(M₂)
class M4 : public Wave {
 public:
  M4()
      : Wave(kM4, 4, -4, 4, 0, 0, 0, 0, 4, -4, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m22) {}
};

/// SN₄ = S₂ + N₂
///
/// V = 4T - 3s + 2h + p
/// u = 2ξ - 2ν
/// f = f(M₂)
class SN4 : public Wave {
 public:
  SN4()
      : Wave(kSN4, 4, -3, 2, 1, 0, 0, 0, 2, -2, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m2) {}
};

/// MS₄ = M₂ + S₂
///
/// V = 4T - 2s + 2h
/// u = +2ξ - 2ν
/// f = f(M₂)
class MS4 : public Wave {
 public:
  MS4()
      : Wave(kMS4, 4, -2, 2, 0, 0, 0, 0, 2, -2, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m2) {}
};

/// MK₄ = M₂ + K₂
///
/// V = 4T - 2s + 4h
/// u = 2ξ - 2ν - 2ν''
/// f = f(MK₄)
class MK4 : public Wave {
 public:
  MK4()
      : Wave(kMK4, 4, -2, 4, 0, 0, 0, 0, 2, -2, -2, 0, kShortPeriod,
             &AstronomicAngle::f_m2_k2) {}
};

/// S₄ = S₂ + S₂
///
/// V = 4T
/// u = 0
/// f = 1
class S4 : public Wave {
 public:
  S4()
      : Wave(kS4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, kShortPeriod,
             &AstronomicAngle::f_1) {}
};

/// SK₄ = S₂ + K₂
///
/// V = 4T + 2h
/// u = -2ν''
/// f = f(K₂)
class SK4 : public Wave {
 public:
  SK4()
      : Wave(kSK4, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, -2, kShortPeriod,
             &AstronomicAngle::f_k2) {}
};

/// R₄ = R₂ + R₂
///
/// V = 4T + 2h - 2p1
/// u = 0
/// f = 1
class R4 : public Wave {
 public:
  R4()
      : Wave(kR4, 4, 0, 2, 0, 0, -2, 0, 0, 0, 0, 0, kShortPeriod,
             &AstronomicAngle::f_1) {}
};

/// 2MN₆ = 2M₂ + N₂
///
/// V = 6T - 7s + 6h + p
/// u = 6ξ - 6ν
/// f = f(M₂)³
class _2MN6 : public Wave {
 public:
  _2MN6()
      : Wave(k2MN6, 6, -7, 6, 1, 0, 0, 0, 6, -6, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m23) {}
};

/// M₆ = 3M₂
///
/// V = 6T - 6s + 6h
/// u = +6ξ - 6ν
/// f = f(M₂)³
class M6 : public Wave {
 public:
  M6()
      : Wave(kM6, 6, -6, 6, 0, 0, 0, 0, 6, -6, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m23) {}
};

/// MSN₆ = M₂ + S₂ + N₂
///
/// V = 6T - 5s + 4h + p
/// u = 4ξ - 4ν
/// f = f(M₂)²
class MSN6 : public Wave {
 public:
  MSN6()
      : Wave(kMSN6, 6, -5, 4, 1, 0, 0, 0, 4, -4, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m22) {}
};

/// 2MS₆ = 2M₂ + S₂
///
/// V = 6T - 4s + 4h
/// u = 4ξ - 4ν
/// f = f(M₂)²
class _2MS6 : public Wave {
 public:
  _2MS6()
      : Wave(k2MS6, 6, -4, 4, 0, 0, 0, 0, 4, -4, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m22) {}
};

/// 2MK₆ = 2M₂ + K₂
///
/// V = 6T - 4s + 6h
/// u = 4ξ - 4ν - 2ν''
/// f = f(M₂)² × f(K₂)
class _2MK6 : public Wave {
 public:
  _2MK6()
      : Wave(k2MK6, 6, -4, 6, 0, 0, 0, 0, 4, -4, 0, -2, kShortPeriod,
             &AstronomicAngle::f_m23_k2) {}
};

/// 2SM₆ = 2S₂ + M₂
///
/// V = 6T - 2s + 2h
/// u = 2ξ - 2ν
/// f = f(M₂)
class _2SM6 : public Wave {
 public:
  _2SM6()
      : Wave(k2SM6, 6, -2, 2, 0, 0, 0, 0, 2, -2, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m2) {}
};

/// MSK₆ = M₂ + K₂ + S₂
///
/// V = 6T - 2s + 4h
/// u = 2ξ - 2ν - 2ν''
/// f = f(M₂) × f(K₂)
class MSK6 : public Wave {
 public:
  MSK6()
      : Wave(kMSK6, 6, -2, 4, 0, 0, 0, 0, 2, -2, -2, 0, kShortPeriod,
             &AstronomicAngle::f_m2_k2) {}
};

/// S₆ = 3S₂
///
/// V = 6T
/// u = 0
/// f = 1
class S6 : public Wave {
 public:
  S6()
      : Wave(kS6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, kShortPeriod,
             &AstronomicAngle::f_1) {}
};

/// M₈ = 4M₂
///
/// V = 8T - 8s + 8h
/// u = 8ξ - 8ν
/// f = f(M₂)⁴
class M8 : public Wave {
 public:
  M8()
      : Wave(kM8, 8, -8, 8, 0, 0, 0, 0, 8, -8, 0, 0, kShortPeriod,
             &AstronomicAngle::f_m24) {}
};

/// MSf = M₂ - S₂
///
/// V = 2s - 2h
/// u = 2ξ - 2ν
/// f = f(M₂) * f(S2) = f(M₂)
///
/// @warning Same frequency as MSf LP : 2s -2h
class MSf : public Wave {
 public:
  MSf()
      : Wave(kMSf, 0, 2, -2, 0, 0, 0, 0, 2, -2, 0, 0, kLongPeriod,
             &AstronomicAngle::f_m2) {}
};

/// Properties of tide waves computed
class WaveTable {
 private:
  std::vector<std::shared_ptr<Wave>> waves_{};

  /// Compute nodal corrections from SCHUREMAN (1958).
  ///
  /// Indexes used in this routine are internal to the code
  /// and corresponds to the "original" ondes.dat file.
  ///
  /// @param a Astronomic angle
  void nodal_a(const AstronomicAngle& a) {
    for (auto& item : waves_) {
      item->nodal_a(a);
    }
  }

  /// Compute nodal corrections from SCHUREMAN (1958).
  ///
  /// Indexes used in this routine are internal to the code and corresponds to
  /// the "original" ondes.dat file.
  ///
  /// @param a Astronomic angle
  void nodal_g(const AstronomicAngle& a) {
    for (auto& item : waves_) {
      item->nodal_g(a);
    }
  }

 public:
  /// Default constructor
  WaveTable(const std::vector<std::string>& waves = {});

  /// Copy constructor
  WaveTable(const WaveTable& wt) {
    auto waves = std::vector<std::string>();
    waves.reserve(wt.size());
    for (const auto& item : wt.waves_) {
      waves.push_back(item->name());
    }
    *this = WaveTable(waves);
  }

  /// Gets the tidal waves known
  static std::vector<std::string> known_constituents();

  /// Compute nodal corrections.
  ///
  /// @param epoch Desired UTC time expressed in number of seconds elapsed since
  /// 1970-01-01T00:00:00
  /// @return the astronomic angle, indicating the date on which the tide is to
  /// be calculated.
  AstronomicAngle compute_nodal_corrections(const double epoch) {
    auto angles = AstronomicAngle(epoch);

    nodal_a(angles);
    nodal_g(angles);

    return angles;
  }

  /// Gets the wave properties
  std::shared_ptr<Wave> wave(const Wave::Ident ident) const {
    auto it = find(ident);
    return it != end() ? *it : nullptr;
  }

  /// Gets the wave properties
  const std::shared_ptr<Wave>& operator[](const size_t index) const {
    return waves_.at(index);
  }

  /// Gets the wave properties
  std::shared_ptr<Wave> wave(const std::string& ident) const {
    auto it = find(ident);
    return it != end() ? *it : nullptr;
  }

  /// Gets the wave properties
  std::shared_ptr<Wave>& wave(const Wave::Ident ident) { return waves_[ident]; }

  /// Returns an iterator to the beginning of the wave table
  std::vector<std::shared_ptr<Wave>>::const_iterator begin() const {
    return waves_.begin();
  }

  /// Returns an iterator to the end of the wave table
  std::vector<std::shared_ptr<Wave>>::const_iterator end() const {
    return waves_.end();
  }

  /// Returns an iterator to the beginning of the wave table
  std::vector<std::shared_ptr<Wave>>::iterator begin() {
    return waves_.begin();
  }

  /// Returns an iterator to the end of the wave table
  std::vector<std::shared_ptr<Wave>>::iterator end() { return waves_.end(); }

  /// Searches the properties of a wave from its name.
  std::vector<std::shared_ptr<Wave>>::const_iterator find(
      const std::string& name) const {
    for (auto it = begin(); it != end(); ++it) {
      if (name == (*it)->name()) {
        return it;
      }
    }
    return end();
  }

  /// Searches the properties of a wave from its identifier.
  std::vector<std::shared_ptr<Wave>>::const_iterator find(
      const Wave::Ident& ident) const {
    for (auto it = begin(); it != end(); ++it) {
      if (ident == (*it)->ident()) {
        return it;
      }
    }
    return end();
  }

  /// Returns the size of the table
  size_t size() const { return waves_.size(); }

  static Eigen::VectorXcd harmonic_analysis(
      const Eigen::Ref<const Eigen::VectorXd>& h,
      const Eigen::Ref<const Eigen::MatrixXd>& f,
      const Eigen::Ref<const Eigen::MatrixXd>& vu);
};