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

PYBIND11_MODULE(core, m) {
  if (!PyDateTimeAPI) {
    PyDateTime_IMPORT;
  }

  m.def("timestamp",
        [](py::handle datetime) -> double { return timestamp(datetime); },
        py::arg("date"), "Return POSIX timestamp as float");

  py::class_<AstronomicAngle>(m, "AstronomicAngle")
      .def(py::init<>())
      .def(py::init<const double>(), py::arg("epoch"))
      .def_property_readonly("t", &AstronomicAngle::t)
      .def_property_readonly("n", &AstronomicAngle::n)
      .def_property_readonly("h", &AstronomicAngle::h)
      .def_property_readonly("s", &AstronomicAngle::s)
      .def_property_readonly("p1", &AstronomicAngle::p1)
      .def_property_readonly("p", &AstronomicAngle::p)
      .def_property_readonly("i", &AstronomicAngle::i)
      .def_property_readonly("xi", &AstronomicAngle::xi)
      .def_property_readonly("nu", &AstronomicAngle::nu)
      .def_property_readonly("x1ra", &AstronomicAngle::x1ra)
      .def_property_readonly("r", &AstronomicAngle::r)
      .def_property_readonly("nuprim", &AstronomicAngle::nuprim)
      .def_property_readonly("nusec", &AstronomicAngle::nusec);

  py::class_<Wave, std::shared_ptr<Wave>> wave(m, "Wave");
  py::enum_<Wave::Ident>(wave, "Ident")
      .value("kMm", Wave::kMm)
      .value("kMf", Wave::kMf)
      .value("kMtm", Wave::kMtm)
      .value("kMsqm", Wave::kMsqm)
      .value("k2Q1", Wave::k2Q1)
      .value("kSigma1", Wave::kSigma1)
      .value("kQ1", Wave::kQ1)
      .value("kRho1", Wave::kRho1)
      .value("kO1", Wave::kO1)
      .value("kMP1", Wave::kMP1)
      .value("kM11", Wave::kM11)
      .value("kM12", Wave::kM12)
      .value("kM13", Wave::kM13)
      .value("kChi1", Wave::kChi1)
      .value("kPi1", Wave::kPi1)
      .value("kP1", Wave::kP1)
      .value("kS1", Wave::kS1)
      .value("kK1", Wave::kK1)
      .value("kPsi1", Wave::kPsi1)
      .value("kPhi1", Wave::kPhi1)
      .value("kTheta1", Wave::kTheta1)
      .value("kJ1", Wave::kJ1)
      .value("kOO1", Wave::kOO1)
      .value("kMNS2", Wave::kMNS2)
      .value("kEps2", Wave::kEps2)
      .value("k2N2", Wave::k2N2)
      .value("kMu2", Wave::kMu2)
      .value("k2MS2", Wave::k2MS2)
      .value("kN2", Wave::kN2)
      .value("kNu2", Wave::kNu2)
      .value("kM2", Wave::kM2)
      .value("kMKS2", Wave::kMKS2)
      .value("kLambda2", Wave::kLambda2)
      .value("kL2", Wave::kL2)
      .value("k2MN2", Wave::k2MN2)
      .value("kT2", Wave::kT2)
      .value("kS2", Wave::kS2)
      .value("kR2", Wave::kR2)
      .value("kK2", Wave::kK2)
      .value("kMSN2", Wave::kMSN2)
      .value("kEta2", Wave::kEta2)
      .value("k2SM2", Wave::k2SM2)
      .value("kMO3", Wave::kMO3)
      .value("k2MK3", Wave::k2MK3)
      .value("kM3", Wave::kM3)
      .value("kMK3", Wave::kMK3)
      .value("kN4", Wave::kN4)
      .value("kMN4", Wave::kMN4)
      .value("kM4", Wave::kM4)
      .value("kSN4", Wave::kSN4)
      .value("kMS4", Wave::kMS4)
      .value("kMK4", Wave::kMK4)
      .value("kS4", Wave::kS4)
      .value("kSK4", Wave::kSK4)
      .value("kR4", Wave::kR4)
      .value("k2MN6", Wave::k2MN6)
      .value("kM6", Wave::kM6)
      .value("kMSN6", Wave::kMSN6)
      .value("k2MS6", Wave::k2MS6)
      .value("k2MK6", Wave::k2MK6)
      .value("k2SM6", Wave::k2SM6)
      .value("kMSK6", Wave::kMSK6)
      .value("kS6", Wave::kS6)
      .value("kM8", Wave::kM8)
      .value("kMSf", Wave::kMSf)
      .value("kSsa", Wave::kSsa)
      .value("kSa", Wave::kSa);

  py::enum_<Wave::TidalType>(wave, "TidalType")
      .value("kLongPeriod", Wave::kLongPeriod)
      .value("kShortPeriod", Wave::kShortPeriod);

  wave.def_property_readonly("ident", &Wave::ident)
      .def_property_readonly("freq", &Wave::freq)
      .def_property_readonly("type", &Wave::type)
      .def_property_readonly("f", &Wave::f)
      .def_property_readonly("u", &Wave::u)
      .def("nodal_a", &Wave::nodal_a, py::arg("a"))
      .def("nodal_g", &Wave::nodal_g, py::arg("a"))
      .def("v0u", &Wave::v0u)
      .def("v0", &Wave::v0)
      .def("name", &Wave::name);

  py::class_<WaveTable>(m, "WaveTable")
      .def(py::init<std::vector<std::string>>(),
           py::arg("waves") = std::vector<std::string>{})
      .def_static("known_constituents", &WaveTable::known_constituents)
      .def("compute_nodal_corrections",
           [](WaveTable& self, double epoch) -> AstronomicAngle {
             return self.compute_nodal_corrections(epoch);
           },
           py::arg("epoch"))
      .def("compute_nodal_corrections",
           [](WaveTable& self, py::array_t<double>& epoch) -> py::tuple {
             if (epoch.ndim() != 1) {
               throw std::invalid_argument(
                   "epoch must be a one-dimensional array");
             }
             py::ssize_t size = self.size();
             py::array_t<double, py::array::c_style> f(
                 py::array::ShapeContainer{epoch.size(), size});
             py::array_t<double, py::array::c_style> v0u(
                 py::array::ShapeContainer{epoch.size(), size});
             {
               py::gil_scoped_release release;

               auto _epoch = epoch.mutable_unchecked<1>();
               auto _f = f.mutable_unchecked<2>();
               auto _v0u = v0u.mutable_unchecked<2>();

               for (py::ssize_t ix = 0; ix < epoch.shape(0); ++ix) {
                 self.compute_nodal_corrections(_epoch(ix));

                 for (std::size_t jx = 0; jx < self.size(); ++jx) {
                   _f(ix, jx) = self[jx]->f();
                   _v0u(ix, jx) = self[jx]->v0u();
                 }
               }
             }
             return py::make_tuple(f, v0u);
           },
           py::arg("epoch"))
      .def("wave",
           [](const WaveTable& self, const Wave::Ident ident)
               -> std::shared_ptr<Wave> { return self.wave(ident); },
           py::arg("ident"))
      .def("wave",
           [](const WaveTable& self, const std::string& ident)
               -> std::shared_ptr<Wave> { return self.wave(ident); },
           py::arg("ident"))
      .def("tide",
           [](WaveTable& self, py::array_t<double>& epoch,
              py::array_t<std::complex<double>>& wave) -> py::array_t<double> {
             if (epoch.ndim() != 1) {
               throw std::invalid_argument(
                   "epoch must be a one-dimensional array");
             }
             if (wave.ndim() != 1) {
               throw std::invalid_argument(
                   "wave must be a one-dimensional array");
             }
             if (static_cast<size_t>(wave.shape(0)) != self.size()) {
               throw std::invalid_argument(
                   "wave must contain as many items as tidal constituents "
                   "loaded");
             }
             py::array_t<double, py::array::c_style> result(
                 py::array::ShapeContainer{epoch.size()});
             {
               py::gil_scoped_release release;

               auto _epoch = epoch.mutable_unchecked<1>();
               auto _wave = wave.mutable_unchecked<1>();
               auto _result = result.mutable_unchecked<1>();

               for (py::ssize_t ix = 0; ix < epoch.shape(0); ++ix) {
                 self.compute_nodal_corrections(_epoch(ix));
                 double tide = 0;

                 for (size_t jx = 0; jx < self.size(); ++jx) {
                   const auto& item = self[jx];
                   double phi = item->v0u();
                   tide += item->f() * (_wave(jx).real() * std::cos(phi) +
                                        _wave(jx).imag() * std::sin(phi));
                 }
                 _result(ix) = tide;
               }
             }
             return result;
           },
           py::arg("epoch"), py::arg("wave"))
      .def_static(
          "harmonic_analysis",
          [](const Eigen::Ref<const Eigen::VectorXd>& h,
             const Eigen::Ref<const Eigen::MatrixXd>& f,
             const Eigen::Ref<const Eigen::MatrixXd>& v0u) -> Eigen::VectorXcd {
            py::gil_scoped_release release;
            return WaveTable::harmonic_analysis(h, f, v0u);
          },
          py::arg("h"), py::arg("f"), py::arg("v0u"))
      .def("__len__", [](const WaveTable& self) { return self.size(); })
      .def("__getitem__",
           [](const WaveTable& self, size_t index) -> std::shared_ptr<Wave> {
             return self[index];
           },
           py::arg("index"))
      .def("__getitem__",
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
      .def("__iter__",
           [](const WaveTable& self) {
             return py::make_iterator(self.begin(), self.end());
           },
           py::keep_alive<0, 1>());
}