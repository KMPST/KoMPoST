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
#ifndef EVENTINPUT_CPP
#define EVENTINPUT_CPP
#include "EventInput.h"
#include <gsl/gsl_math.h>
#include <iostream>

#include "INIReader.h"

namespace EventInput {

double afm = 0.08;
int Ns = 512;

int xSTART = 0;
int xEND = Ns - 1;

int ySTART = 0;
int yEND = Ns - 1;

void Setup(INIReader &reader) {
  afm = reader.GetReal("EventInput", "afm", afm);
  Ns = reader.GetInteger("EventInput", "Ns", Ns);

  xSTART = reader.GetInteger("EventInput", "xSTART", xSTART);
  xEND = reader.GetInteger("EventInput", "xEND", xEND);

  ySTART = reader.GetInteger("EventInput", "ySTART", ySTART);
  yEND = reader.GetInteger("EventInput", "yEND", yEND);

  std::cerr << "** EventInput ** Initialized a grid layout:\n" 
            << "  afm    = " << afm << "\n"
            << "  Ns     = " << Ns << "\n"
            << "  xSTART = " << xSTART << "\n"
            << "  xEND   = " << xEND << "\n"
            << "  ySTART = " << ySTART << "\n"
            << "  yEND   = " << yEND <<  std::endl;
}
}

namespace KoMPoSTParameters {

double EtaOverS = 2. / (4 * M_PI);
double EtaOverSTemperatureScale = 0.1; // GeV cuttoff temperature scale

std::string Regulator("TwoPass");

double Sigma;

int EVOLUTION_MODE = 1;
int ENERGY_PERTURBATIONS = 1;
int MOMENTUM_PERTURBATIONS = 1;

void Setup(INIReader &reader) {

  // Determines the scaling variable
  EtaOverS = reader.GetReal("KoMPoSTParameters", "EtaOverS", EtaOverS);

  EtaOverSTemperatureScale = reader.GetReal(
      "KoMPoSTParameters", "EtaOverSTemperatureScale", EtaOverSTemperatureScale);

  Regulator = reader.GetString("KoMPoSTParameters", "Regulator", Regulator);

  //Sigma = reader.GetReal("KoMPoSTParameters", "Sigma", Sigma);

  EVOLUTION_MODE =
      reader.GetInteger("KoMPoSTParameters", "EVOLUTION_MODE", EVOLUTION_MODE);

  ENERGY_PERTURBATIONS = reader.GetInteger(
      "KoMPoSTParameters", "ENERGY_PERTURBATIONS", ENERGY_PERTURBATIONS);

  MOMENTUM_PERTURBATIONS = reader.GetInteger(
      "KoMPoSTParameters", "MOMENTUM_PERTURBATIONS", MOMENTUM_PERTURBATIONS);

  std::cerr << "** EventInput ** Initialized KoMPoST parameters:\n" 
            << "  EtaOverS                 = " << EtaOverS << "\n"
            << "  EtaOverSTemperatureScale = " << EtaOverSTemperatureScale << "\n"
            << "  Regulator                = " << Regulator << "\n"
            << "  EVOLUTION_MODE           = " << EVOLUTION_MODE << " ; 0 -- free streaming, 1 -- EKT\n"
            << "  ENERGY_PERTURBATIONS     = " << ENERGY_PERTURBATIONS << "\n"
            << "  MOMENTUM_PERTURBATIONS   = " << MOMENTUM_PERTURBATIONS << "\n"
            << "  Sigma                    = " << Sigma <<  std::endl;
}
}
#endif
