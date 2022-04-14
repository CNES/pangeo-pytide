// Copyright (c) 2020 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
#pragma once

namespace speed {

/// @def julian_century
/// @brief Number of seconds in one Julian century
static constexpr double julian_century = 3155760000.0;

/// @brief Compute the speed in degree by hour for the moon's mean longitude (s)
constexpr double s() noexcept {
  return (((1336.0 * 360.0 + 307.892) / julian_century)) * 3600.0;
}

/// @brief Compute the speed in degree by hour for the sun's mean longitude (h)
constexpr double h() noexcept {
  return (((100.0 * 360.0 + 0.769) / julian_century)) * 3600.0;
}

/// @brief Compute the speed in degree by hour for the longitude of moon's
/// perigee (p)
constexpr double p() noexcept {
  return (((11.0 * 360.0 + 109.032) / julian_century)) * 3600.0;
}

/// @brief Compute the speed in degree by hour for the longitude of moon's
/// ascending node (N′)
constexpr double n() noexcept {
  return (((-5.0 * 360.0 - 134.142) / julian_century)) * 3600.0;
}

/// @brief Compute the speed in degree by hour for the longitude of sun's
/// perigee (p₁)
constexpr double p1() noexcept { return ((1.719 / julian_century)) * 3600; }

/// @brief Compute the speed in degree by hour for the local mean lunar time (τ)
constexpr double tau() noexcept { return 15.0 - s() + h(); }

} // namespace speed

namespace frequency {

/// @brief Compute the frequency in degree by hour for the moon's mean longitude
/// (s)
constexpr double s() noexcept { return 1.0 / ((15.0 / speed::s()) * 24.0); }

/// @brief Compute the frequency in degree by hour for the sun's mean longitude
/// (h)
constexpr double h() noexcept { return 1.0 / ((15.0 / speed::h()) * 24.0); }

/// @brief Compute the frequency in degree by hour for the longitude of moon's
/// perigee (p)
constexpr double p() noexcept { return 1.0 / ((15.0 / speed::p()) * 24.0); }

/// @brief Compute the frequency in degree by hour for the longitude of moon's
/// ascending node (N′)
constexpr double n() noexcept { return 1.0 / ((15.0 / speed::n()) * 24.0); }

/// @brief Compute the frequency in degree by hour for the longitude of sun's
/// perigee (p₁)
constexpr double p1() noexcept { return 1.0 / ((15.0 / speed::p1()) * 24.0); }

/// @brief Compute the frequency in degree by hour for the local mean lunar time
/// (τ)
constexpr double tau() noexcept { return 1.0 / ((15.0 / speed::tau()) * 24.0); }

} // namespace frequency
