// Copyright (c) 2019 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
#pragma once
#include <cmath>

/// PI
template <typename T>
constexpr T pi() {
  return std::atan2(T(0), T(-1));
}

/// PI/2
template <typename T>
constexpr T pi_2() {
  return 0.5 * pi<T>();
}

/// 2 * PI
template <typename T>
constexpr T two_pi() {
  return T(2) * pi<T>();
}

/// Square a number.
///
/// @return \f$x^2\f$
template <typename T>
constexpr T sqr(const T& x) {
  return x * x;
}

/// Power 3.
///
/// @return \f$x^3\f$
template <typename T>
constexpr T pow3(const T& x) {
  return sqr(x) * x;
}

/// Power 4.
///
/// @return \f$x^4\f$
template <typename T>
constexpr T pow4(const T& x) {
  return sqr(x) * sqr(x);
}

/// Convert angle x from radians to degrees.
template <typename T>
constexpr T radians(const T& x) {
  return x * pi<T>() / T(180);
}

/// Convert angle x from degrees to radians.
template <typename T>
constexpr T degrees(const T& x) {
  return x * T(180) / pi<T>();
}

/// Normalize an angle.
///
/// @param x The angle in degrees.
/// @param min Minimum circle value
/// @return the angle reduced to the range [min, 360 + min[
template <typename T>
constexpr T normalize_angle(const T& x, const T& min = T(-180)) {
  T result = std::remainder(x - min, T(360));
  if (result < 0) {
    result += T(360);
  }
  result += min;
  return result;
}