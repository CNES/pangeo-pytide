// Copyright (c) 2020 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
#include "astronomic_angle.hpp"
#include "math.hpp"

/// Handle the number of days elapsed since 1900-01-01T00:00:00.000 UTC
static auto const REFERENCE = -25567;

AstronomicAngle::AstronomicAngle(const double epoch) {
  auto jc = ((epoch / 86400.0) - REFERENCE) / 36525.0;
  auto jc_squared = sqr(jc);

  // T mean solar angle relative to Greenwich
  t_ =
      std::fmod(pi<double>() + two_pi<double>() * 36525 * jc, two_pi<double>());

  // SCHUREMAN FORMULAE P. 162 (oder 2 is enough)

  // Longitude of moon's node (N)
  n_ = radians(normalize_angle(
      259.1560564 - 1934.1423972 * jc + 0.0021056 * jc_squared, 0.0));

  // Mean longitude of sun (h)
  h_ = radians(normalize_angle(
      280.1895014 + 36000.768925 * jc + 0.0003025 * jc_squared, 0.0));

  // Mean longitude of moon (s)
  s_ = radians(normalize_angle(
      277.0256206 + 481267.892 * jc + 0.002525 * jc_squared, 0.0));

  // Longitude of solar perigee (p₁)
  p1_ = radians(normalize_angle(
      281.2208569 + 1.719175 * jc + 0.0004528 * jc_squared, 0.0));

  // Longitude of lunar perigee (p)
  p_ = radians(normalize_angle(
      334.3837215 + 4069.0322056 * jc - 0.0103444 * jc_squared, 0.0));

  // SCHUREMAN FORMULAE P. 156
  auto u = 0.913694997 - 0.035692561 * std::cos(n_);

  // Inclination of the moon's orbit to the celestial equator
  i_ = std::acos(u);

  auto tgn2 = std::tan(n_ * 0.5);
  auto at1 = std::atan(1.01883 * tgn2);
  auto at2 = std::atan(0.64412 * tgn2);

  // Longitude in moon's orbit of lunar intersection
  xi_ = -at1 - at2 + n_;

  if (n_ > pi<double>()) {
    xi_ -= two_pi<double>();
  }

  // Right ascension of lunar intersection
  nu_ = at1 - at2;

  // for constituents l2, k1, k2
  auto tgi2 = std::tan(i_ * 0.5);

  // SCHUREMAN P. 41 (191)
  // Mean longitude  of the lunar perigee reckoned from the lunar intersection
  auto p = p_ - xi_;

  // SCHUREMAN P. 44 (213)
  x1ra_ =
      std::sqrt(1.0 - 12.0 * sqr(tgi2) * std::cos(2.0 * p) + 36.0 * pow4(tgi2));

  // SCHUREMAN P. 41 (196)
  r_ = std::atan(std::sin(2.0 * p) /
                 (1.0 / (6.0 * sqr(tgi2)) - std::cos(2.0 * p)));

  // SCHUREMAN P. 45 (224)
  nuprim_ = std::atan(std::sin(2.0 * i_) * std::sin(nu_) /
                      (std::sin(2.0 * i_) * std::cos(nu_) + 0.3347));

  // SCHUREMAN P. 46 (232)
  nusec_ = 0.5 * std::atan((sqr(std::sin(i_)) * std::sin(2.0 * nu_)) /
                           (sqr(std::sin(i_)) * std::cos(2.0 * nu_) + 0.0727));
}
