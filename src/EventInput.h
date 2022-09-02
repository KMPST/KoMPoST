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
#ifndef EVENTINPUT_H
#define EVENTINPUT_H
#ifndef M_HBARC
#define M_HBARC 0.197326979
#endif
#include <string>

class INIReader;

namespace EventInput {

extern double afm;
extern int Ns;

extern int xSTART;
extern int xEND;

extern int ySTART;
extern int yEND;

void Setup(INIReader &reader);
}

namespace KoMPoSTParameters {
extern std::string KineticTheory;

extern double EtaOverS;
extern double EtaOverSTemperatureScale;
extern double Sigma;

extern std::string Regulator;

extern int EVOLUTION_MODE;
extern int ENERGY_PERTURBATIONS;
extern int MOMENTUM_PERTURBATIONS;
extern int DECOMPOSITION_METHOD;

void Setup(INIReader &reader);
}
#endif
