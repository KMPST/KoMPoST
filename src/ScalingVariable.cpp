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
////////////////////////////////////////////////////////////////////////

#include "gsl/gsl_integration.h"
#include "gsl/gsl_math.h"
#include <cstdio>
#ifndef M_HBARC
#define M_HBARC 0.197326979
#endif
#include "ScalingVariable.h"

ScalingVariable::ScalingVariable(double etaovers0, double lambda) {
  fEtaOverS0 = etaovers0;
  fLambda = lambda;
  ClearEstimate();
}

ScalingVariable::ScalingVariable(const ScalingVariable &sv) {

  fEtaOverS0 = sv.fEtaOverS0;
  fLambda = sv.fLambda;

  EstimateIsSet = sv.EstimateIsSet;
  fT00BG = sv.fT00BG;
  fDT00 = sv.fDT00;
  fGx = sv.fGx;
  fGy = sv.fGy;
}

//! Returns a smooth (compact function) interpolating between one and zero,
//! i.e. for  xstart < x < x + w the function transition from one to zero
double InterpolatingFunction(const double &x, const double &xstart,
                             const double &width) {
  double u = (x - xstart) / width;
  if (u < 0.) {
    return 1.;
  }
  if (u > 1.) {
    return 0.;
  }
  if (u < 0.5) {
    return 1 - 2. * u * u;
  } else {
    return 2. * (1 - u) * (1 - u);
  }
}

//! Returns an estimate of the scaling variable
double ScalingVariable::ScalingVar(const double &tin_fm, const double &tout_fm,
                                   const double &K) {

  double TimeIn = tin_fm / M_HBARC;
  double TimeOut = tout_fm / M_HBARC;

  double result = 0.;

  // Evaluates the scaling variable with T(t) = K /(K*tau)^1/3:
  //
  // result =   3/2 * \int_tau1^(tau2) dt/(eta/(s(T) T)) .
  //
  // The shear viscosity to entropy is
  //
  // \eta/s (T) =   fEtaOverS0 (1 + Lambda^2/T^2)
  {
    double t2 = TimeOut * K;
    double t1 = TimeIn * K;
    double L2 = GSL_MAX(pow(fLambda / K, 2), GSL_DBL_EPSILON);
    result = 1 / (fEtaOverS0 * L2) *
             gsl_log1p((pow(t2, 2. / 3.) * L2 - pow(t1, 2. / 3.) * L2) /
                       (1. + pow(t1, 2. / 3.) * L2));
  }

  // If this is the second pass, we then modify (reduce) the scaling variable
  // stored in result, based on the local conditions as measured by x1 and x2
  // which
  // should be small if non-linear effects are small.
  if (EstimateIsSet) {
    // Given the estimated results evaluate a discriminant z,
    // which gives a good indication if the inversion will have trouble
    double x1 = sqrt(fGx * fGx + fGy * fGy) / (fT00BG * 4. / 3.);
    double x2 = fDT00 / fT00BG;
    const double zw = x1 - 2. / 3. * x2;

    const double zw_trouble = 0.4;
    const double width = 0.3;
    // Reduce the scaling variable (result) if the
    // if a the criterion measured by z = x1 - 2/3 x2 is greater than
    // 0.4.  By 0.7 the interpolating (zw_trouble + width) the interpolating
    // function goes to zero, i.e. we are evaluating the
    result *= InterpolatingFunction(zw, zw_trouble, width);
    // If the scaling variable is too low the code crashes.
    result = (result < 0.1 ? 0.1 : result);

    return result;
  } else {
    // This is the first pass return the scaling variable as it comes.
    return result;
  }
}
