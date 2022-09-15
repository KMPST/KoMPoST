/*
 * Copyright (c) 2018, Aleksi Kurkela, Aleksas Mazeliauskas, Jean-Francois
 * Paquet, Soeren Schlichting and Derek Teaney
 * All rights reserved.
 *
 * KoMPoST is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/KMPST/KoMPoST/
 */

#include <gsl/gsl_math.h>
#include <iostream>
#ifndef M_HBARC
#define M_HBARC 0.197326979
#endif

namespace GreensFunctions {

class LowKArguments;

//////////////////////////////////////////////////

namespace EnergyPerturbations {

namespace FreeStreaming {
namespace CoordinateSpace {
double Gs(double dX, double dT);
double Gv(double dX, double dT);
double Gd(double dX, double dT);
double Gr(double dX, double dT);
}
} // FreeStreaming

namespace KineticTheory {
namespace CoordinateSpace {
double Gs(double dX, double dT, double ScalingVarOut);
double Gv(double dX, double dT, double ScalingVarOut);
double Gd(double dX, double dT, double ScalingVarOut);
double Gr(double dX, double dT, double ScalingVarOut);
}
} // KineticTheory

namespace LowKLimit {
namespace CoordinateSpace {
double Gs(double dX, double dT, const LowKArguments &in);
double Gv(double dX, double dT, const LowKArguments &in);
double Gd(double dX, double dT, const LowKArguments &in);
double Gr(double dX, double dT, const LowKArguments &in);
}
} // LowKLimit

} // EnergyPerturbations

//////////////////////////////////////////////////

namespace MomentumPerturbations {

namespace FreeStreaming {
namespace CoordinateSpace {
double Hv(double dX, double dT);
double Hd(double dX, double dT);
double Hr(double dX, double dT);
double Htd(double dX, double dT);
double Htm(double dX, double dT);
double Htr(double dX, double dT);
}
} // FreeStreaming

namespace KineticTheory {
namespace CoordinateSpace {
double Hv(double dX, double dT, double ScalingVarOut);
double Hd(double dX, double dT, double ScalingVarOut);
double Hr(double dX, double dT, double ScalingVarOut);
double Htd(double dX, double dT, double ScalingVarOut);
double Htm(double dX, double dT, double ScalingVarOut);
double Htr(double dX, double dT, double ScalingVarOut);
}
} // Kinetic Theory

} // MomentumPerturbations

//////////////////////////////////////////////////

void Setup(double Sigma, int NumberOfPoints, int ENERGY_PERTURBATIONS,
           int MOMENTUM_PERTURBATIONS);

void Output(int ENERGY_PERTURBATIONS, int MOMENTUM_PERTURBATIONS);

class LowKArguments {

  // private:
public:
  // Number of DOF
  const int NuG = 16;

  double tauF;  // Final time in fm
  double tauIn; // Initial time in fm

  double E0;
  double EF;

  double TXX0;
  double TXXF;

  double EtaByS;
  double SigmaBG; // Sigma in fm

public:
  LowKArguments(double tf, double tin, double e0, double ef, double txx0,
                double txxf, double etabys, double sigmabg)
      : tauF(tf), tauIn(tin), E0(e0), EF(ef), TXX0(txx0), TXXF(txxf),
        EtaByS(etabys), SigmaBG(sigmabg) {}

  double S(double dX, double dT) const {
    const double &s = SigmaBG;
    double rs = dX / s;
    return 1. / (M_PI * s * s) * exp(-rs * rs);
  }

  double DS(double dX, double dT) const {
    const double &s = SigmaBG;
    double rs = dX / s;
    return dT / (M_PI * s * s * s) * exp(-rs * rs) * (-2. * rs);
  }

  double DSOverR(double dX, double dT) const {
    const double &s = SigmaBG;
    double rs = dX / s;
    return dT * dT / (M_PI * s * s * s * s) * exp(-rs * rs) * (-2.);
  }

  double DDS(double dX, double dT) const {
    const double &s = SigmaBG;
    double rs = dX / s;
    return dT * dT / (M_PI * s * s * s * s) * exp(-rs * rs) *
           (4 * rs * rs - 2.);
  }

  double Ass() const {
    double dT = tauF - tauIn;
    double s2 = SigmaBG * SigmaBG / (dT * dT);
    return (-0.5 + s2 / 2.) * E0 / EF * (EF + TXXF) / (E0 + TXX0);
  }

  double Asv() const { return 0.5 * E0 / EF * (EF + TXXF) / (E0 + TXX0); }

  double ScalingVariableInverse() const {
    return 3. / 4. * EtaByS * Entropy(EF) * M_HBARC / (EF * tauF);
  }

  double Entropy(double e) const {
    double T = pow(30. / (M_PI * M_PI * KoMPoSTParameters::NuG) * e, 0.25);
    return 4. / 3. * e / T;
  }
};

} // GreensFunctions
