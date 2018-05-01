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

#ifndef M_HBARC
#define M_HBARC 0.197326979
#endif


//! This computes the coordinate space response
//! functions in the low k limit.  To use

namespace CoordinateSpace {

double Gs(double dX, double dT, const GreensFunctions::LowKArguments &in) {
  return   -0.5 * in.Ass() * (in.DSOverR(dX, dT) + in.DDS(dX, dT));
}

double Gv(double dX, double dT, const GreensFunctions::LowKArguments &in) {
  return -1.0 * in.Asv() * in.DS(dX, dT);
}

// Gij=Gd \delta^{ij} + Gr r^{i}r^{j}/r^2 //
double Gd(double dX, double dT, const GreensFunctions::LowKArguments &in) {
  double ass = in.Ass();
  double asv = in.Asv();
  double u = in.ScalingVariableInverse();
  return (in.DSOverR(dX, dT) + in.DDS(dX, dT)) *
             (-1.0 * ass / 6. + 2. / 3 * u * (-1.0 *ass / 2. - asv)) +
         in.DSOverR(dX, dT) * asv * 2. * u;
}

double Gr(double dX, double dT, const GreensFunctions::LowKArguments &in) {
  double asv = in.Asv();
  double u = in.ScalingVariableInverse();
  return (-in.DSOverR(dX, dT) + in.DDS(dX, dT)) * asv * 2. * u;
}

} // CoordinateSpace
