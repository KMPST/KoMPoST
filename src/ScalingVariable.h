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
#ifndef CMPOSER_ScalingVariable_h
#define CMPOSER_ScalingVariable_h
#include "gsl/gsl_integration.h"
#ifndef M_HBARC
#define M_HBARC 0.197326979
#endif

class ScalingVariable {

private:
  double fEtaOverS0;
  double fLambda;

  bool EstimateIsSet; // True if the first pass is completed.
  double fT00BG;      // The background energy density from the first pass
  double fDT00;       // The energy perturbation from the first pass
  double fGx;         // The x-momentum perturbation from the first pass
  double fGy;         // The y-momentum perturbation from the first pass

public:
  ScalingVariable(double etaovers0, double lambda);
  ScalingVariable(const ScalingVariable &sv);
  ~ScalingVariable() { ; }

  void ClearEstimate() { EstimateIsSet = false; }

  //! Sets the estimate of the size of the initial perturbations from
  //! the first pass.  This estimate is used to locally modify the scaling
  //! varibale in the second pass.
  void SetEstimate(double t00bg, double dt00, double gx, double gy) {
    EstimateIsSet = true;
    fT00BG = t00bg;
    fDT00 = dt00;
    fGx = gx;
    fGy = gy;
  }

  double GetEtaOverS0() { return fEtaOverS0; }

  double ScalingVar(const double &tout_fm, const double &kvalue) {
    return ScalingVar(0., tout_fm, kvalue);
  }

  //!  Returns the scaling variable:
  //!
  //!  If EstimateIsSet is false, then the naive scaling variable is returned
  //!  for a given K value. Specifically we take a form
  //!
  //!  2/3 Integral[ dtau /( eta/s(T) T ), tau1, tau2) ]
  //!
  //!  With T(tau) = K /(K * tau)^(1/3).  Similarly we take a modified form
  //!  of eta/s
  //!
  //!  eta/s (T) = (eta/s)_0 (1  + Lambda^2/T^2)
  //!
  //!  which has a non-conformal cutoff Lambda.
  //!
  //!  If EstimateIsSet is true (which indicates that the energy and momentum
  //!  perturbations were calculated from the first pass) then  the naive
  //!  scaling variable is reduced (driving the system towards free streaming)
  //!  based on the size of the perturbations.
  double ScalingVar(const double &tin_fm, const double &tout_fm,
                    const double &kvalue);
};

#endif // CMPOSER_ScalingVariable_h
