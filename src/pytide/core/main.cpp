// Copyright (c) 2019 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <datetime.h>
#include "astronomic_angle.hpp"
#include "wave.hpp"

namespace py = pybind11;

/// Returns number of days since civil 1970-01-01.  Negative values indicate
/// days prior to 1970-01-01.
///
/// http://howardhinnant.github.io/date_algorithms.html
static inline int64_t days_from_civil(int y, unsigned m, unsigned d) {
  y -= static_cast<int>(m <= 2);
  const auto era = (y >= 0 ? y : y - 399) / 400;
  const auto yoe = static_cast<unsigned>(y - era * 400);            // [0, 399]
  const auto doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;  // [0, 365]
  const auto doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;  // [0, 146096]
  return era * 146097LL + static_cast<int64_t>(doe) - 719468LL;
}

/// Return POSIX timestamp as double
double timestamp(py::handle datetime) {
  if (!datetime) {
    throw std::invalid_argument(
        "a datetime.datetime is required (got type null)");
  }

  if (PyDateTime_Check(datetime.ptr())) {
    if (reinterpret_cast<_PyDateTime_BaseTZInfo*>(datetime.ptr())->hastzinfo) {
      throw std::invalid_argument(
          "only the naive datetime object can be converted to timestamp");
    }

    auto sec = days_from_civil(PyDateTime_GET_YEAR(datetime.ptr()),
                               PyDateTime_GET_MONTH(datetime.ptr()),
                               PyDateTime_GET_DAY(datetime.ptr())) *
               86400;
    sec += PyDateTime_DATE_GET_HOUR(datetime.ptr()) * 3600 +
           PyDateTime_DATE_GET_MINUTE(datetime.ptr()) * 60 +
           PyDateTime_DATE_GET_SECOND(datetime.ptr());

    return sec + PyDateTime_DATE_GET_MICROSECOND(datetime.ptr()) * 1e-6;
  }
  throw std::invalid_argument(
      "a datetime.datetime is required (got type " +
      std::string(py::str(datetime.get_type().attr("__name__"))) + ")");
}

/// Calculates the tide of a given time series.
static py::array_t<double> tide_from_time_series(
    WaveTable& self, const Eigen::Ref<const Eigen::VectorXd>& epoch,
    const Eigen::Ref<const Eigen::VectorXcd>& wave) {
  if (static_cast<size_t>(wave.rows()) != self.size()) {
    throw std::invalid_argument(
        "wave must contain as many items as tidal constituents loaded");
  }

  py::array_t<double, py::array::c_style> result(
      py::array::ShapeContainer{epoch.rows()});
  auto _result = result.mutable_unchecked<1>();
  {
    py::gil_scoped_release release;

    for (py::ssize_t ix = 0; ix < epoch.rows(); ++ix) {
      double tide = 0;
      self.compute_nodal_corrections(epoch(ix));

      for (size_t jx = 0; jx < self.size(); ++jx) {
        const auto& item = self[jx];
        double phi = item->vu();

        tide += item->f() * (wave(jx).real() * std::cos(phi) +
                             wave(jx).imag() * std::sin(phi));
      }
      _result(ix) = tide;
    }
  }
  return result;
}

/// Calculates the tide for a given date from a grid describing the properties
/// of tidal waves over an area of interest.
static py::array_t<double> tide_from_mapping(
    WaveTable& self, const double epoch,
    const Eigen::Ref<const Eigen::MatrixXcd>& wave) {
  if (static_cast<size_t>(wave.rows()) != self.size()) {
    throw std::invalid_argument(
        "the first dimension of wave must contain as many items as "
        "tidal constituents loaded");
  }
  self.compute_nodal_corrections(epoch);
  py::array_t<double, py::array::c_style> result(
      py::array::ShapeContainer{wave.cols()});
  auto _result = result.mutable_unchecked<1>();
  {
    py::gil_scoped_release release;

    for (py::ssize_t ix = 0; ix < wave.cols(); ++ix) {
      double tide = 0;
      for (size_t jx = 0; jx < self.size(); ++jx) {
        const auto& item = self[jx];
        double phi = item->vu();

        tide += item->f() * (wave(jx, ix).real() * std::cos(phi) +
                             wave(jx, ix).imag() * std::sin(phi));
      }
      _result(ix) = tide;
    }
  }
  return result;
}

PYBIND11_MODULE(core, m) {
  if (!PyDateTimeAPI) {
    PyDateTime_IMPORT;
  }

  m.doc() = R"__doc__(
Core module
-----------
)__doc__";

  m.def(
      "timestamp",
      [](py::handle datetime) -> double { return timestamp(datetime); },
      py::arg("date"), "Return POSIX timestamp as float");

  py::class_<AstronomicAngle>(m, "AstronomicAngle")
      .def(py::init<>())
      .def(py::init<const double>(), py::arg("epoch"), R"__doc__(
Initialize some astronomic data useful for nodal corrections.

Args:
  epoch (float, optional): Desired UTC time)__doc__")
      .def_property_readonly("t", &AstronomicAngle::t, "Hour angle of mean sun")
      .def_property_readonly("n", &AstronomicAngle::n,
                             "Longitude of moon's node")
      .def_property_readonly("h", &AstronomicAngle::h,
                             "Mean longitude of the sun")
      .def_property_readonly("s", &AstronomicAngle::s,
                             "Mean longitude of the moon")
      .def_property_readonly("p1", &AstronomicAngle::p1,
                             "Mean longitude of solar perigee")
      .def_property_readonly("p", &AstronomicAngle::p,
                             "Mean longitude of lunar perigee")
      .def_property_readonly(
          "i", &AstronomicAngle::i,
          "Obliquity of lunar orbit with respect to earth's equator")
      .def_property_readonly("xi", &AstronomicAngle::xi,
                             "Longitude in moon's orbit of lunar intersection")
      .def_property_readonly("nu", &AstronomicAngle::nu,
                             "Right ascension of lunar intersection")
      .def_property_readonly("x1ra", &AstronomicAngle::x1ra,
                             "Factor in amplitude of constituent :math:`L_{2}`")
      .def_property_readonly("r", &AstronomicAngle::r,
                             "Term in argument of constituent :math:`L_{2}`")
      .def_property_readonly(
          "nuprim", &AstronomicAngle::nuprim,
          "Term in argument of lunisolar constituent :math:`K_{1}`")
      .def_property_readonly(
          "nusec", &AstronomicAngle::nusec,
          "Term in argument of lunisolar constituent :math:`K_{2}`");

  py::class_<Wave, std::shared_ptr<Wave>> wave(m, "Wave", "Wave definition");
  py::enum_<Wave::Ident>(wave, "Ident")
      .value("kMm", Wave::kMm, ":math:`Mm`")
      .value("kMf", Wave::kMf, ":math:`Mf`")
      .value("kMtm", Wave::kMtm, ":math:`Mtm`")
      .value("kMsqm", Wave::kMsqm, ":math:`Msqm`")
      .value("k2Q1", Wave::k2Q1, ":math:`2Q_{1}`")
      .value("kSigma1", Wave::kSigma1, ":math:`\\sigma_{1}`")
      .value("kQ1", Wave::kQ1, ":math:`Q_{1}`")
      .value("kRho1", Wave::kRho1, ":math:`\\rho_{1}`")
      .value("kO1", Wave::kO1, ":math:`O_{1}`")
      .value("kMP1", Wave::kMP1, ":math:`MP_{1}`")
      .value("kM11", Wave::kM11, ":math:`M_{11}`")
      .value("kM12", Wave::kM12, ":math:`M_{12}`")
      .value("kM13", Wave::kM13, ":math:`M_{13}`")
      .value("kChi1", Wave::kChi1, ":math:`\\chi_{1}`")
      .value("kPi1", Wave::kPi1, ":math:`\\pi_{1}`")
      .value("kP1", Wave::kP1, ":math:`P_{1}`")
      .value("kS1", Wave::kS1, ":math:`S_{1}`")
      .value("kK1", Wave::kK1, ":math:`K_{1}`")
      .value("kPsi1", Wave::kPsi1, ":math:`\\psi_{1}`")
      .value("kPhi1", Wave::kPhi1, ":math:`\\varphi_{1}`")
      .value("kTheta1", Wave::kTheta1, ":math:`\\theta_{1}`")
      .value("kJ1", Wave::kJ1, ":math:`J_{1}`")
      .value("kOO1", Wave::kOO1, ":math:`OO_{1}`")
      .value("kMNS2", Wave::kMNS2, ":math:`MNS_{2}`")
      .value("kEps2", Wave::kEps2, ":math:`\\varepsilon_{2}`")
      .value("k2N2", Wave::k2N2, ":math:`2N_{2}`")
      .value("kMu2", Wave::kMu2, ":math:`\\upsilon_{2}`")
      .value("k2MS2", Wave::k2MS2, ":math:`2MS_{2}`")
      .value("kN2", Wave::kN2, ":math:`N_{2}`")
      .value("kNu2", Wave::kNu2, ":math:`\\nu_{2}`")
      .value("kM2", Wave::kM2, ":math:`M_{2}`")
      .value("kMKS2", Wave::kMKS2, ":math:`MKS_{2}`")
      .value("kLambda2", Wave::kLambda2, ":math:`\\lambda_{2}`")
      .value("kL2", Wave::kL2, ":math:`L_{2}`")
      .value("k2MN2", Wave::k2MN2, ":math:`2MN_{2}`")
      .value("kT2", Wave::kT2, ":math:`T_{2}`")
      .value("kS2", Wave::kS2, ":math:`S_{2}`")
      .value("kR2", Wave::kR2, ":math:`R_{2}`")
      .value("kK2", Wave::kK2, ":math:`K_{2}`")
      .value("kMSN2", Wave::kMSN2, ":math:`MSN_{2}`")
      .value("kEta2", Wave::kEta2, ":math:`\\eta_{2}`")
      .value("k2SM2", Wave::k2SM2, ":math:`2SM_{2}`")
      .value("kMO3", Wave::kMO3, ":math:`MO_{3}`")
      .value("k2MK3", Wave::k2MK3, ":math:`2MK_{3}`")
      .value("kM3", Wave::kM3, ":math:`M_{3}`")
      .value("kMK3", Wave::kMK3, ":math:`MK_{3}`")
      .value("kN4", Wave::kN4, ":math:`N_{4}`")
      .value("kMN4", Wave::kMN4, ":math:`MN_{4}`")
      .value("kM4", Wave::kM4, ":math:`M_{4}`")
      .value("kSN4", Wave::kSN4, ":math:`SN_{4}`")
      .value("kMS4", Wave::kMS4, ":math:`MS_{4}`")
      .value("kMK4", Wave::kMK4, ":math:`MK_{4}`")
      .value("kS4", Wave::kS4, ":math:`S_{4}`")
      .value("kSK4", Wave::kSK4, ":math:`SK_{4}`")
      .value("kR4", Wave::kR4, ":math:`R_{4}`")
      .value("k2MN6", Wave::k2MN6, ":math:`2MN_{6}`")
      .value("kM6", Wave::kM6, ":math:`M_{6}`")
      .value("kMSN6", Wave::kMSN6, ":math:`MSN_{6}`")
      .value("k2MS6", Wave::k2MS6, ":math:`2MS_{6}`")
      .value("k2MK6", Wave::k2MK6, ":math:`2MK_{6}`")
      .value("k2SM6", Wave::k2SM6, ":math:`2SM_{6}`")
      .value("kMSK6", Wave::kMSK6, ":math:`MSK_{6}`")
      .value("kS6", Wave::kS6, ":math:`S_{6}`")
      .value("kM8", Wave::kM8, ":math:`M_{8}`")
      .value("kMSf", Wave::kMSf, ":math:`MSf`")
      .value("kSsa", Wave::kSsa, ":math:`Ssa`")
      .value("kSa", Wave::kSa, ":math:`Sa`");

  py::enum_<Wave::TidalType>(wave, "TidalType", "Possible type of tidal wave")
      .value("kLongPeriod", Wave::kLongPeriod, "Long period tidal waves")
      .value("kShortPeriod", Wave::kShortPeriod, "Short period tidal waves");

  wave.def_property_readonly("ident", &Wave::ident, R"__doc__(
Gets the wave ident
)__doc__")
      .def_property_readonly("freq", &Wave::freq, R"__doc__(
Gets the wave frequency (radians per seconds)
)__doc__")
      .def_property_readonly("type", &Wave::type, R"__doc__(
Gets the wave type
)__doc__")
      .def_property_readonly("f", &Wave::f, R"__doc__(
Gets the nodal correction for amplitude
)__doc__")
      .def_property_readonly("u", &Wave::u, R"__doc__(
Gets the nodal correction for phase
)__doc__")
      .def("nodal_a", &Wave::nodal_a, py::arg("a"), R"__doc__(
Compute nodal corrections from SCHUREMAN (1958).

Args:
  a (pytide.core.AstronomicAngle): Astronomic angle
)__doc__")
      .def("nodal_g", &Wave::nodal_g, py::arg("a"), R"__doc__(
Compute nodal corrections from SCHUREMAN (1958).

Args:
  a (pytide.core.AstronomicAngle): Astronomic angle
)__doc__")
      .def("vu", &Wave::vu, R"__doc__(
Gets :math:`v` (greenwich argument) + :math:`u` (nodal correction for phase)
)__doc__")
      .def("v", &Wave::v, R"__doc__(
Gets :math:`v` (greenwich argument)
)__doc__")
      .def("name", &Wave::name, R"__doc__(
Gets the wave name
)__doc__");

  py::class_<WaveTable>(m, "WaveTable", "Properties of tide waves computed")
      .def(py::init<std::vector<std::string>>(),
           py::arg("waves") = std::vector<std::string>{})
      .def_static("known_constituents", &WaveTable::known_constituents,
                  "Gets the tidal waves known by this object")
      .def(
          "compute_nodal_corrections",
          [](WaveTable& self, double epoch) -> AstronomicAngle {
            return self.compute_nodal_corrections(epoch);
          },
          py::arg("epoch"), R"__doc__(
Compute nodal corrections.

Args:
  epoch (float): Desired UTC time expressed in number of seconds elapsed since
    1970-01-01T00:00:00
Returns:
  pytide.core.AstronomicAngle: The astronomic angle, indicating the date on
    which the tide is to be calculated.
)__doc__")
      .def(
          "compute_nodal_modulations",
          [](WaveTable& self, py::array_t<double>& epoch) -> py::tuple {
            if (epoch.ndim() != 1) {
              throw std::invalid_argument(
                  "epoch must be a one-dimensional array");
            }
            py::ssize_t size = self.size();
            py::array_t<double, py::array::c_style> f(
                py::array::ShapeContainer{size, epoch.size()});
            py::array_t<double, py::array::c_style> vu(
                py::array::ShapeContainer{size, epoch.size()});
            {
              py::gil_scoped_release release;

              auto _epoch = epoch.mutable_unchecked<1>();
              auto _f = f.mutable_unchecked<2>();
              auto _vu = vu.mutable_unchecked<2>();

              for (py::ssize_t ix = 0; ix < epoch.shape(0); ++ix) {
                self.compute_nodal_corrections(_epoch(ix));

                for (std::size_t jx = 0; jx < self.size(); ++jx) {
                  _f(jx, ix) = self[jx]->f();
                  _vu(jx, ix) = self[jx]->vu();
                }
              }
            }
            return py::make_tuple(f, vu);
          },
          py::arg("epoch"), R"__doc__(
Compute nodal modulations for amplitude and phase.

Args:
  epoch (numpy.ndarray): Desired UTC time expressed in number of seconds
    elapsed since 1970-01-01T00:00:00
Returns:
  tuple: the nodal correction for amplitude, v (greenwich argument) + u
  (nodal correction for phase)
)__doc__")
      .def(
          "wave",
          [](const WaveTable& self, const Wave::Ident ident)
              -> std::shared_ptr<Wave> { return self.wave(ident); },
          py::arg("ident"), "Gets the wave properties")
      .def(
          "wave",
          [](const WaveTable& self, const std::string& ident)
              -> std::shared_ptr<Wave> { return self.wave(ident); },
          py::arg("ident"), "Gets the wave properties")
      .def(
          "tide_from_tide_series",
          [](WaveTable& self, const Eigen::Ref<const Eigen::VectorXd>& epoch,
             const Eigen::Ref<const Eigen::VectorXcd>& wave)
              -> py::array_t<double> {
            return tide_from_time_series(self, epoch, wave);
          },
          py::arg("epoch"), py::arg("wave"), R"__doc__(
Calculates the tide of a given time series.

Args:
  epoch (numpy.ndarray): Time series dates.
  wave (numpy.ndarray): Tidal wave properties.

Return:
  numpy.ndarray:
    The tide calculated for the time series provided.
)__doc__")
      .def(
          "tide_from_mapping",
          [](WaveTable& self, const double epoch,
             const Eigen::Ref<const Eigen::MatrixXcd>& wave)
              -> py::array_t<double> {
            return tide_from_mapping(self, epoch, wave);
          },
          py::arg("epoch"), py::arg("wave"), R"__doc__(
Calculates the tide for a given tidal wave mapping.

Args:
  epoch (float): Mapping date
  wave (numpy.ndarray): A matrix containing the wave properties for each point
    on the map.

Return:
  numpy.ndarray:
    The tide calculated on the area of interest provided.
)__doc__")
      .def_static(
          "harmonic_analysis",
          [](const Eigen::Ref<const Eigen::VectorXd>& h,
             const Eigen::Ref<const Eigen::MatrixXd>& f,
             const Eigen::Ref<const Eigen::MatrixXd>& vu) -> Eigen::VectorXcd {
            py::gil_scoped_release release;
            return WaveTable::harmonic_analysis(h, f, vu);
          },
          py::arg("h"), py::arg("f"), py::arg("vu"), R"__doc__(
Harmonic Analysis

Args:
  h (numpy.ndarray): Sea level.
  f (numpy.ndarray): Nodal correction coefficient applied to the amplitude of
    the constituents analyzed.
  vu (numpy.ndarray): Astronomical argument at time :math:`t` + the nodal
    correction coefficient applied to the phase of the constituents
    analyzed

Returns:
  numpy.ndarray: The complex number representing the different reconstructed
  waves.
)__doc__")
      .def("__len__", [](const WaveTable& self) { return self.size(); })
      .def(
          "__getitem__",
          [](const WaveTable& self, size_t index) -> std::shared_ptr<Wave> {
            return self[index];
          },
          py::arg("index"))
      .def(
          "__getitem__",
          [](const WaveTable& self,
             py::slice slice) -> std::vector<std::shared_ptr<Wave>> {
            size_t start, stop, step, slicelength;
            if (!slice.compute(self.size(), &start, &stop, &step,
                               &slicelength)) {
              throw py::error_already_set();
            }
            auto result = std::vector<std::shared_ptr<Wave>>(slicelength);
            for (size_t ix = 0; ix < slicelength; ++ix) {
              result[ix] = self[start];
              start += step;
            }
            return result;
          },
          py::arg("slice"))
      .def(
          "__iter__",
          [](const WaveTable& self) {
            return py::make_iterator(self.begin(), self.end());
          },
          py::keep_alive<0, 1>());
}