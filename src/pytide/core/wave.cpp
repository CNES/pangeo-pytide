// Copyright (c) 2022 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
#include "wave.hpp"
#include "astronomic_angle.hpp"
#include "math.hpp"
#include <Eigen/Dense>

std::string Wave::name() const {
  switch (ident_) {
  case Wave::kO1:
    return "O1";
  case Wave::kP1:
    return "P1";
  case Wave::kK1:
    return "K1";
  case Wave::k2N2:
    return "2N2";
  case Wave::kMu2:
    return "Mu2";
  case Wave::kN2:
    return "N2";
  case Wave::kNu2:
    return "Nu2";
  case Wave::kM2:
    return "M2";
  case Wave::kL2:
    return "L2";
  case Wave::kT2:
    return "T2";
  case Wave::kS2:
    return "S2";
  case Wave::kK2:
    return "K2";
  case Wave::kM4:
    return "M4";
  case Wave::kS1:
    return "S1";
  case Wave::kQ1:
    return "Q1";
  case Wave::kMm:
    return "Mm";
  case Wave::kMf:
    return "Mf";
  case Wave::kMtm:
    return "Mtm";
  case Wave::kMsqm:
    return "Msqm";
  case Wave::kEps2:
    return "Eps2";
  case Wave::kLambda2:
    return "Lambda2";
  case Wave::kEta2:
    return "Eta2";
  case Wave::k2Q1:
    return "2Q1";
  case Wave::kSigma1:
    return "Sigma1";
  case Wave::kRho1:
    return "Rho1";
  case Wave::kM11:
    return "M11";
  case Wave::kM12:
    return "M12";
  case Wave::kChi1:
    return "Chi1";
  case Wave::kPi1:
    return "Pi1";
  case Wave::kPhi1:
    return "Phi1";
  case Wave::kTheta1:
    return "Theta1";
  case Wave::kJ1:
    return "J1";
  case Wave::kOO1:
    return "OO1";
  case Wave::kM3:
    return "M3";
  case Wave::kM6:
    return "M6";
  case Wave::kMN4:
    return "MN4";
  case Wave::kMS4:
    return "MS4";
  case Wave::kN4:
    return "N4";
  case Wave::kR2:
    return "R2";
  case Wave::kR4:
    return "R4";
  case Wave::kS4:
    return "S4";
  case Wave::kMNS2:
    return "MNS2";
  case Wave::kM13:
    return "M13";
  case Wave::kMK4:
    return "MK4";
  case Wave::kSN4:
    return "SN4";
  case Wave::kSK4:
    return "SK4";
  case Wave::k2MN6:
    return "2MN6";
  case Wave::k2MS6:
    return "2MS6";
  case Wave::k2MK6:
    return "2MK6";
  case Wave::kMSN6:
    return "MSN6";
  case Wave::k2SM6:
    return "2SM6";
  case Wave::kMSK6:
    return "MSK6";
  case Wave::kMP1:
    return "MP1";
  case Wave::k2SM2:
    return "2SM2";
  case Wave::kPsi1:
    return "Psi1";
  case Wave::k2MS2:
    return "2MS2";
  case Wave::kMKS2:
    return "MKS2";
  case Wave::k2MN2:
    return "2MN2";
  case Wave::kMSN2:
    return "MSN2";
  case Wave::kMO3:
    return "MO3";
  case Wave::k2MK3:
    return "2MK3";
  case Wave::kMK3:
    return "MK3";
  case Wave::kS6:
    return "S6";
  case Wave::kM8:
    return "M8";
  case Wave::kMSf:
    return "MSf";
  case Wave::kSsa:
    return "Ssa";
  case Wave::kSa:
    return "Sa";
  default:
    return "unknown";
  }
}

std::vector<std::string> WaveTable::known_constituents() {
  return {"O1",   "P1",   "K1",   "2N2",  "Mu2",     "N2",   "Nu2",    "M2",
          "L2",   "T2",   "S2",   "K2",   "M4",      "S1",   "Q1",     "Mm",
          "Mf",   "Mtm",  "Msqm", "Eps2", "Lambda2", "Eta2", "2Q1",    "Sigma1",
          "Rho1", "M11",  "M12",  "Chi1", "Pi1",     "Phi1", "Theta1", "J1",
          "OO1",  "M3",   "M6",   "MN4",  "MS4",     "N4",   "R2",     "R4",
          "S4",   "MNS2", "M13",  "MK4",  "SN4",     "SK4",  "2MN6",   "2MS6",
          "2MK6", "MSN6", "2SM6", "MSK6", "MP1",     "2SM2", "Psi1",   "2MS2",
          "MKS2", "2MN2", "MSN2", "MO3",  "2MK3",    "MK3",  "S6",     "M8",
          "MSf",  "Ssa",  "Sa"};
}

static std::shared_ptr<Wave> wave_factory(const std::string &name) {
  if (name == "O1") {
    return std::shared_ptr<Wave>(new O1());
  } else if (name == "P1") {
    return std::shared_ptr<Wave>(new P1());
  } else if (name == "K1") {
    return std::shared_ptr<Wave>(new K1());
  } else if (name == "2N2") {
    return std::shared_ptr<Wave>(new _2N2());
  } else if (name == "Mu2") {
    return std::shared_ptr<Wave>(new Mu2());
  } else if (name == "N2") {
    return std::shared_ptr<Wave>(new N2());
  } else if (name == "Nu2") {
    return std::shared_ptr<Wave>(new Nu2());
  } else if (name == "M2") {
    return std::shared_ptr<Wave>(new M2());
  } else if (name == "L2") {
    return std::shared_ptr<Wave>(new L2());
  } else if (name == "T2") {
    return std::shared_ptr<Wave>(new T2());
  } else if (name == "S2") {
    return std::shared_ptr<Wave>(new S2());
  } else if (name == "K2") {
    return std::shared_ptr<Wave>(new K2());
  } else if (name == "M4") {
    return std::shared_ptr<Wave>(new M4());
  } else if (name == "S1") {
    return std::shared_ptr<Wave>(new S1());
  } else if (name == "Q1") {
    return std::shared_ptr<Wave>(new Q1());
  } else if (name == "Mm") {
    return std::shared_ptr<Wave>(new Mm());
  } else if (name == "Mf") {
    return std::shared_ptr<Wave>(new Mf());
  } else if (name == "Mtm") {
    return std::shared_ptr<Wave>(new Mtm());
  } else if (name == "Msqm") {
    return std::shared_ptr<Wave>(new Msqm());
  } else if (name == "Eps2") {
    return std::shared_ptr<Wave>(new Eps2());
  } else if (name == "Lambda2") {
    return std::shared_ptr<Wave>(new Lambda2());
  } else if (name == "Eta2") {
    return std::shared_ptr<Wave>(new Eta2());
  } else if (name == "2Q1") {
    return std::shared_ptr<Wave>(new _2Q1());
  } else if (name == "Sigma1") {
    return std::shared_ptr<Wave>(new Sigma1());
  } else if (name == "Rho1") {
    return std::shared_ptr<Wave>(new Rho1());
  } else if (name == "M11") {
    return std::shared_ptr<Wave>(new M11());
  } else if (name == "M12") {
    return std::shared_ptr<Wave>(new M12());
  } else if (name == "Chi1") {
    return std::shared_ptr<Wave>(new Chi1());
  } else if (name == "Pi1") {
    return std::shared_ptr<Wave>(new Pi1());
  } else if (name == "Phi1") {
    return std::shared_ptr<Wave>(new Phi1());
  } else if (name == "Theta1") {
    return std::shared_ptr<Wave>(new Theta1());
  } else if (name == "J1") {
    return std::shared_ptr<Wave>(new J1());
  } else if (name == "OO1") {
    return std::shared_ptr<Wave>(new OO1());
  } else if (name == "M3") {
    return std::shared_ptr<Wave>(new M3());
  } else if (name == "M6") {
    return std::shared_ptr<Wave>(new M6());
  } else if (name == "MN4") {
    return std::shared_ptr<Wave>(new MN4());
  } else if (name == "MS4") {
    return std::shared_ptr<Wave>(new MS4());
  } else if (name == "N4") {
    return std::shared_ptr<Wave>(new N4());
  } else if (name == "R2") {
    return std::shared_ptr<Wave>(new R2());
  } else if (name == "R4") {
    return std::shared_ptr<Wave>(new R4());
  } else if (name == "S4") {
    return std::shared_ptr<Wave>(new S4());
  } else if (name == "MNS2") {
    return std::shared_ptr<Wave>(new MNS2());
  } else if (name == "M13") {
    return std::shared_ptr<Wave>(new M13());
  } else if (name == "MK4") {
    return std::shared_ptr<Wave>(new MK4());
  } else if (name == "SN4") {
    return std::shared_ptr<Wave>(new SN4());
  } else if (name == "SK4") {
    return std::shared_ptr<Wave>(new SK4());
  } else if (name == "2MN6") {
    return std::shared_ptr<Wave>(new _2MN6());
  } else if (name == "2MS6") {
    return std::shared_ptr<Wave>(new _2MS6());
  } else if (name == "2MK6") {
    return std::shared_ptr<Wave>(new _2MK6());
  } else if (name == "MSN6") {
    return std::shared_ptr<Wave>(new MSN6());
  } else if (name == "2SM6") {
    return std::shared_ptr<Wave>(new _2SM6());
  } else if (name == "MSK6") {
    return std::shared_ptr<Wave>(new MSK6());
  } else if (name == "MP1") {
    return std::shared_ptr<Wave>(new MP1());
  } else if (name == "2SM2") {
    return std::shared_ptr<Wave>(new _2SM2());
  } else if (name == "Psi1") {
    return std::shared_ptr<Wave>(new Psi1());
  } else if (name == "2MS2") {
    return std::shared_ptr<Wave>(new _2MS2());
  } else if (name == "MKS2") {
    return std::shared_ptr<Wave>(new MKS2());
  } else if (name == "2MN2") {
    return std::shared_ptr<Wave>(new _2MN2());
  } else if (name == "MSN2") {
    return std::shared_ptr<Wave>(new MSN2());
  } else if (name == "MO3") {
    return std::shared_ptr<Wave>(new MO3());
  } else if (name == "2MK3") {
    return std::shared_ptr<Wave>(new _2MK3());
  } else if (name == "MK3") {
    return std::shared_ptr<Wave>(new MK3());
  } else if (name == "S6") {
    return std::shared_ptr<Wave>(new S6());
  } else if (name == "M8") {
    return std::shared_ptr<Wave>(new M8());
  } else if (name == "MSf") {
    return std::shared_ptr<Wave>(new MSf());
  } else if (name == "Ssa") {
    return std::shared_ptr<Wave>(new Ssa());
  } else if (name == "Sa") {
    return std::shared_ptr<Wave>(new Sa());
  }

  throw std::runtime_error("The tidal wave is unknown: " + name);
}

WaveTable::WaveTable(const std::vector<std::string> &waves) {
  const auto &wave_list = waves.empty() ? known_constituents() : waves;
  for (auto &item : wave_list) {
    waves_.emplace_back(wave_factory(item));
  }
}

void Wave::nodal_g(const AstronomicAngle &a) {
  v_ = argument_[0] * a.t() + argument_[1] * a.s() + argument_[2] * a.h() +
       argument_[3] * a.p() + argument_[5] * a.p1() +
       argument_[6] * pi_2<double>();
  u_ = argument_[7] * a.xi() + argument_[8] * a.nu() +
       argument_[9] * a.nuprim() + argument_[10] * a.nusec();
}

Eigen::VectorXcd
WaveTable::harmonic_analysis(const Eigen::Ref<const Eigen::VectorXd> &h,
                             const Eigen::Ref<const Eigen::MatrixXd> &f,
                             const Eigen::Ref<const Eigen::MatrixXd> &vu) {
  if (f.rows() != vu.rows() || f.cols() != vu.cols()) {
    throw std::invalid_argument(
        "f and vu could not be broadcast together with shape (" +
        std::to_string(f.rows()) + ", " + std::to_string(f.cols()) + ") (" +
        std::to_string(vu.rows()) + ", " + std::to_string(vu.cols()) + ")");
  }

  if (h.rows() != f.cols() || h.rows() != vu.cols()) {
    throw std::invalid_argument(
        "f, vu could not be broadcast with h with shape (" +
        std::to_string(f.rows()) + ", " + std::to_string(f.cols()) + ") (" +
        std::to_string(h.cols()) + ")");
  }
  auto w_size = f.rows();
  auto result = Eigen::VectorXcd(w_size);

  if (h.hasNaN()) {
    result.fill(std::complex<double>(std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::quiet_NaN()));
    return result;
  }

  auto H = Eigen::MatrixXd(w_size << 1, h.rows());

  H.topRows(w_size) = f.array() * vu.array().cos();
  H.bottomRows(w_size) = f.array() * vu.array().sin();

  auto solution = ((H * H.transpose()).inverse() * H) * h;
  result.real() = solution.topRows(w_size);
  result.imag() = solution.bottomRows(w_size);

  return result;
}
