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

#ifndef BACKGROUNDEVOLUTION_H
#define BACKGROUNDEVOLUTION_H

class ScalingVariable;

namespace BackgroundEvolution {

namespace KineticTheory {

double DetermineScalingFactor(double EIn, double tIn, ScalingVariable &sv, double &ScalingVarIn);

void Propagate(double tOut, double KValue, double ScalingVarOut, double EtaOverS, double &T00, double &TXX, double &TYY, double &TZZ);

} // Kinetic theory

namespace FreeStreaming {

void Propagate(double T00In, double TXXIn, double TYYIn, double TZZIn,
               double tIn, double tOut, double &T00, double &TXX, double &TYY,
               double &TZZ);
} // Free streaming

} // Background 

#endif
