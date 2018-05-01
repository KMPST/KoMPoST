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
#include "KineticEvolution.h"
#include "ScalingVariable.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>


// ENERGY CUT-OFF //
static double ENERGY_CUTOFF = 1E-5;

// EVENT INPUT PARAMETERS //
#include "EventInput.h"

#include "EnergyMomentumTensor.h"

// BACKGROUND EVOLUTION //
#include "BackgroundEvolution.h"

// GREENS FUNCTIONS FOR ENERGY-MOMENTUM PERTRUBATIONS //
#include "GreensFunctions.h"

namespace KoMPoST {

//! Store the perturbations from the first pass, so that the they can be used
//! for the second pass.  This is only relevant if KoMPoSTParameters::Regulator is
//! "TwoPass". The relevant information for each cite is stored in
//! TOutBG->CellData  CellData[9 -- 12].
//!
//! The program runs in two passes. In the first pass we simply record the
//! values of the unregulated output. This is used to modify the scaling
//! variable in the second pass
void PrepareKoMPoSTEstimate(EnergyMomentumTensorMap *TOutBG,
                        EnergyMomentumTensorMap *TOutFull) {
  using namespace EventInput;

  for (int yS = xSTART; yS <= xEND; yS++) {
    for (int xS = xSTART; xS <= xEND; xS++) {
      double t00bg = TOutBG->Get(0, 0, xS, yS);
      double t00 = TOutFull->Get(0, 0, xS, yS);
      double dt00 = t00 - t00bg;
      double t0x = TOutFull->Get(0, 1, xS, yS);
      double t0y = TOutFull->Get(0, 2, xS, yS);

      // Store the initial estimate for the second pass.
      TOutBG->SetCellData(9, xS, yS, t00bg);
      TOutBG->SetCellData(10, xS, yS, dt00);
      TOutBG->SetCellData(11, xS, yS, t0x);
      TOutBG->SetCellData(12, xS, yS, t0y);
    }
  }
}

//! Regulate the size of the perturbations based on how large they are.  This
//! is relevant if KoMPoSTParameters::Regulator = "KoMPoSTAddition"
//!
//! The routine takes the unregulated result (*TOutFull) and the background
//! (*TOutBg), and modifies (*TOutFull) by modifying the size of the
//! perturbations according to the regulating procedure. In this case we simply
//! increase T00 (and also TXX, TYY, TYY) until a discriminant is satisfied.
void RegulateKoMPoSTAddition(EnergyMomentumTensorMap *TOutBG,
                         EnergyMomentumTensorMap *TOutFull) {
  using namespace EventInput;

  for (int yS = xSTART; yS <= xEND; yS++) {
    for (int xS = xSTART; xS <= xEND; xS++) {
      // Extract the stress
      double t00bg = TOutBG->Get(0, 0, xS, yS);
      double t00 = TOutFull->Get(0, 0, xS, yS);
      double txx = TOutFull->Get(1, 1, xS, yS);
      double tyy = TOutFull->Get(2, 2, xS, yS);
      double tzz = TOutFull->Get(3, 3, xS, yS);

      double dt00 = t00 - t00bg;
      double gx = TOutFull->Get(0, 1, xS, yS);
      double gy = TOutFull->Get(0, 2, xS, yS);

      // Store the initial values for the stress tensor
      TOutBG->SetCellData(9, xS, yS, t00bg);
      TOutBG->SetCellData(10, xS, yS, dt00);
      TOutBG->SetCellData(11, xS, yS, gx);
      TOutBG->SetCellData(12, xS, yS, gy);

      // Construct the discriminant
      double x1 = sqrt(gx * gx + gy * gy) / (t00bg * 4. / 3.);
      double x2 = dt00 / t00bg;
      double z = x1 - 2. / 3. * x2;

      // Implementing a "Breaking Function" so the regulated
      // z satisfies z < zstart + width
      const double zstart = 0.32;
      const double width = 0.15;
      double znew = z;
      if (z > zstart) {
        znew = zstart + width * std::tanh((z - zstart) / width);
      }
      double delta_t00 = 3. / 2. * (x1 - znew) * t00bg - dt00;

      TOutFull->SetComponent(0, 0, xS, yS, t00bg + dt00 + delta_t00);
      TOutFull->SetComponent(1, 1, xS, yS, txx + 1. / 3. * delta_t00);
      TOutFull->SetComponent(2, 2, xS, yS, tyy + 1. / 3. * delta_t00);
      TOutFull->SetComponent(3, 3, xS, yS, tzz + 1. / 3. * delta_t00);
    }
  }
}

// COMPUTE EVOLUTION OF THE BACKGROUND ENERGY-MOMENTUM TENSOR.
//
// Given the initial condition *TIn we compute an average initial energy
// density in a causal circle set by SigmaBG.
//
// This is used to determine the scaling variable x at the initial time and the
// final time at each point in the  transverse plane. The scaling variable is
// computed by the Scaler class, which knows about eta/s etc.
//
// Finally the average vale of the stress tensor at each point is stored in
// TOutBG at a time tOut.
//
// On output the TOutBG structure contains the background stress tensor.
//
// It also contains additional information about each grid point which is
// stored in the CellData of the TOutBG structure, CellsData(0...8), which can
// be used for diagnostics, but is also essential for propagating the
// perturbations.   In particular the essential information is the average
// input energy density and TXX, and the scaling variable.
//
// The detailed information about what information is stored is documented
// in the code below -- see SetCellData.
void ComputeBackground(bool IsFirstPass, EnergyMomentumTensorMap *TIn,
                       EnergyMomentumTensorMap *TOutBG,
                       const ScalingVariable &ScalerIn, double SigmaBG,
                       int EVOLUTION_MODE) {

  using namespace EventInput;

  // GET RELEVANT TIMES //
  double tIn = TIn->tau;
  double tOut = TOutBG->tau;

  // COMMANDLINE OUTPUT //
  std::cerr << "#COMPUTING EVOLUTION FROM " << tIn << " fm/c TO " << tOut
            << " fm/c" << std::endl;

#pragma omp parallel for
  for (int yS = ySTART; yS <= yEND; yS++) {
    // Class responsible for computing the scaling varaible.  We create a copy
    // for independently running  parallel process.
    ScalingVariable Scaler(ScalerIn);

    // COMMANDLINE OUTPUT -- PROGRESS MONITOR //
    if (omp_get_thread_num() == 0) {
      std::cerr << "#PROGRESS IS " << int(100*( (xEND - xSTART + 1) * (yS - ySTART) )/((double) (xEND - xSTART + 1) * (yEND - ySTART + 1) /
                       omp_get_max_threads() )) << "\%"
                << std::endl;
    }

    for (int xS = xSTART; xS <= xEND; xS++) {

      //////////////////////////
      // BACKGROUND EVOLUTION //
      //////////////////////////

      // ENERGY AVERAGED WITH GAUSSIAN PROFILE  //
      double T00InAvg = 0.0;

      // Minus Seond derivative of the averaged energy density background
      // divided by  (T00 + T^xx) ) =  - d_i d^i T00/(T00 + TXX)  . This is
      // needed for the low-k limit.
      double ddT00InAvg = 0.0;
      // The derivative of averaged T00 in the x direction / (T00 + TXX)
      double dxT00InAvg = 0.0;
      // The derivative of averaged T00 in the y direction / (T00 + TXX)
      double dyT00InAvg = 0.0;

      double TXXInAvg = 0.0;
      double TYYInAvg = 0.0;
      double TZZInAvg = 0.0;

      double sigma = SigmaBG;
      const double nsigma = 4.0;
      int range = int(nsigma * sigma / afm) + 1;

      // LIMITS OF THE GAUSSIAN PROFILE //
      int xstart = std::max(xS - range, 0);
      int xend = std::min(xS + range, Ns);
      int ystart = std::max(yS - range, 0);
      int yend = std::min(yS + range, Ns);

      double Normalization = 0.;
      for (int yE = ystart; yE < yend; yE++) {
        for (int xE = xstart; xE < xend; xE++) {
          // Compute an average background energy density in a causal patch
          double rsquare =
              ((xE - xS) * (xE - xS) + (yE - yS) * (yE - yS)) * afm * afm;

          double weight = afm * afm * exp(-rsquare / (sigma * sigma)) /
                          (M_PI * sigma * sigma);

          double dxweight = weight * -2. * (xE - xS) * afm / (sigma * sigma);
          double dyweight = weight * -2. * (yE - yS) * afm / (sigma * sigma);

          double ddweight =
              4. * weight / pow(sigma, 4) * (rsquare - pow(sigma, 2));

          T00InAvg += weight * TIn->Get(0, 0, xE, yE);

          dxT00InAvg += dxweight * TIn->Get(0, 0, xE, yE);
          dyT00InAvg += dyweight * TIn->Get(0, 0, xE, yE);
          ddT00InAvg += ddweight * TIn->Get(0, 0, xE, yE);

          TXXInAvg += weight * TIn->Get(1, 1, xE, yE);
          TYYInAvg += weight * TIn->Get(2, 2, xE, yE);
          TZZInAvg += weight * TIn->Get(3, 3, xE, yE);

          Normalization += weight;
        }
      }
      // NORMALIZE //
      T00InAvg /= Normalization;
      TXXInAvg /= Normalization;
      TYYInAvg /= Normalization;
      TZZInAvg /= Normalization;

      ddT00InAvg /= (Normalization * (T00InAvg + TXXInAvg));
      dxT00InAvg /= (Normalization * (T00InAvg + TXXInAvg));
      dyT00InAvg /= (Normalization * (T00InAvg + TXXInAvg));

      // EVOLUTION OF BACKGROUND ENERGY DENSITY  //
      double T00BG = 0.0;
      double TXXBG = 0.0;
      double TYYBG = 0.0;
      double TZZBG = 0.0;
      double BackgroundKValue = 0.0;

      // CHECK CUT-OFF CRITERION //
      if (T00InAvg < ENERGY_CUTOFF) {
        T00InAvg = ENERGY_CUTOFF;
        TXXInAvg = 0.5 * ENERGY_CUTOFF;
        TYYInAvg = 0.5 * ENERGY_CUTOFF;
        TZZInAvg = 0.;
      }

      // KINETIC-THEORY //
      if (EVOLUTION_MODE == 1 || EVOLUTION_MODE == 2) {

        if (IsFirstPass) {
          Scaler.ClearEstimate();
        } else {
          // This is the second pass. We use the results from the first pass to
          // modify the second pass.  The relevant results from the first  pass
          // were stored in CellData(9...12) by PrepareKoMPoSTAddition.
          double t00bg = TOutBG->GetCellData(9, xS, yS);
          double dt00 = TOutBG->GetCellData(10, xS, yS);
          double t0x = TOutBG->GetCellData(11, xS, yS);
          double t0y = TOutBG->GetCellData(12, xS, yS);
          Scaler.SetEstimate(t00bg, dt00, t0x, t0y);
        }

        double ScalingVarIn;
        BackgroundKValue =
            BackgroundEvolution::KineticTheory::DetermineScalingFactor(
                T00InAvg, tIn, Scaler, ScalingVarIn);

        // ScalingVarBG is an abosolution scaling variable with tIn = 0
        double ScalingVarBG = Scaler.ScalingVar(tOut, BackgroundKValue);
        double EtaByS0 = Scaler.GetEtaOverS0();

        BackgroundEvolution::KineticTheory::Propagate(
            tOut, BackgroundKValue, ScalingVarBG, EtaByS0, T00BG, TXXBG, TYYBG,
            TZZBG);

        // The perturbations are evolved with a difference in scaling variables
        // tOut**2/3 - tI**2/3 which is used for the perturbations, and not for
        // the background.
        double ScalingVarOut = Scaler.ScalingVar(tIn, tOut, BackgroundKValue);

        // Store information about the lattice cite in TOutBG's CellData.
        TOutBG->SetCellData(0, xS, yS, T00InAvg);

        // Store the K value and scaling var in CellData
        TOutBG->SetCellData(1, xS, yS, BackgroundKValue);

        // Store the scaling variable on output
        TOutBG->SetCellData(2, xS, yS, ScalingVarOut);

        // The average TXX pressure on in input
        TOutBG->SetCellData(3, xS, yS, TXXInAvg);

        TOutBG->SetCellData(4, xS, yS, EtaByS0);
        TOutBG->SetCellData(5, xS, yS, SigmaBG);

        // Minus Seond derivative of the averaged energy density background
        // divided by  (T00 + T^xx) ) =  - d_i d^i T00/(T00 + TXX)  . This and
        // the other derivatives are needed  for the low-k limit.
        TOutBG->SetCellData(6, xS, yS, ddT00InAvg);
        // The derivative of averaged T00 in the x direction / (T00 + TXX). T
        TOutBG->SetCellData(7, xS, yS, dxT00InAvg);
        // The derivative of averaged T00 in the y direction / (T00 + TXX)
        TOutBG->SetCellData(8, xS, yS, dyT00InAvg);

      }

      // FREE-STREAMING //
      else if (EVOLUTION_MODE == 0) {

        // COMPUTE EVOLUTION AND CHECK MATCHING EFFICICENCY //
        BackgroundEvolution::FreeStreaming::Propagate(
            T00InAvg, TXXInAvg, TYYInAvg, TZZInAvg, tIn, tOut, T00BG, TXXBG,
            TYYBG, TZZBG);

        // Store information about the lattice cite in TOutBG's CellData
        // structure -- see above. Some of this information is not relevant for
        // the free streaming case, in which case it is set to zero.
        TOutBG->SetCellData(0, xS, yS, T00InAvg);

        // Kvalue and scaling var make no sense for free sreaming
        TOutBG->SetCellData(1, xS, yS, 0.);
        TOutBG->SetCellData(2, xS, yS, 0.);
        TOutBG->SetCellData(3, xS, yS, TXXInAvg);
        TOutBG->SetCellData(4, xS, yS, 0.);
        TOutBG->SetCellData(5, xS, yS, SigmaBG);
        TOutBG->SetCellData(6, xS, yS, ddT00InAvg);
        TOutBG->SetCellData(7, xS, yS, dxT00InAvg);
        TOutBG->SetCellData(8, xS, yS, dyT00InAvg);

      }

      else {
        std::cerr << "#ERROR -- EVOLUTION MODE NOT SPECIFICED CORRECTLY"
                  << std::endl;
        exit(0);
      }

      TOutBG->Set(xS, yS, T00BG, TXXBG, TYYBG, TZZBG, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0);

    } // Loop over x
  }   // Loop over y

  std::cerr<< "#DONE" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

//! Propagate the perturbations.
//!
//! For a given backgound TOutBG propagate the perturbations and add the
//! perturbations as the come out to TOutBG to produce TOutFull.
void ComputePerturbations(EnergyMomentumTensorMap *TIn,
                          EnergyMomentumTensorMap *TOutBG,
                          EnergyMomentumTensorMap *TOutFull,
                          int ENERGY_PERTURBATIONS, int MOMENTUM_PERTURBATIONS,
                          int EVOLUTION_MODE, double CircleRadius) {

  using namespace EventInput;

  // GET RELEVANT TIMES //
  double tIn = TIn->tau;
  double tOut = TOutFull->tau;

  // COMMANDLINE OUTPUT //
  std::cerr << "#COMPUTING EVOLUTION FROM " << tIn << " fm/c TO " << tOut
            << " fm/c" << std::endl;
  // EVOLUTION TIME AND RADIUS OF CAUSAL CIRCLE //
  double EvolutionTime = (tOut - tIn);

#pragma omp parallel for
  for (int yS = ySTART; yS <= yEND; yS++) {

    // COMMANDLINE OUTPUT -- PROGRESS MONITOR //
    if (omp_get_thread_num() == 0) {
      std::cerr << "#PROGRESS IS " << int(100*( (xEND - xSTART + 1) * (yS - ySTART) )/((double) (xEND - xSTART + 1) * (yEND - ySTART + 1) /
                       omp_get_max_threads() )) << "\%"
                << std::endl;
    }

    for (int xS = xSTART; xS <= xEND; xS++) {

      // Get information about the background
      double T00BG = TOutBG->Get(0, 0, xS, yS);
      double TXXBG = TOutBG->Get(1, 1, xS, yS);
      double TYYBG = TOutBG->Get(2, 2, xS, yS);
      double TZZBG = TOutBG->Get(3, 3, xS, yS);

      // Extract extra-information about the lattice point stored in the
      // background structure TOutBG.
      double T00InAvg = TOutBG->GetCellData(0, xS, yS);
      double TXXInAvg = TOutBG->GetCellData(3, xS, yS);
      double ScalingVarOut = TOutBG->GetCellData(2, xS, yS);

      // The values of the perturbations
      double T00Pert = 0.0;
      double TXXPert = 0.0;
      double TYYPert = 0.0;
      double TZZPert = 0.0;
      double T0XPert = 0.0;
      double T0YPert = 0.0;
      double TXYPert = 0.0;

      // CHECK CUT-OFF CRITERION //
      if (T00InAvg <= ENERGY_CUTOFF) {
        // If below the cutoff ignore the perturbations
        goto PERTURBATION_FINISHUP;
      }

      // CHECK WETHER ENERGY/MOMENTUM PERTURBATIONS ARE INCLUDED //
      if (!(ENERGY_PERTURBATIONS || MOMENTUM_PERTURBATIONS)) {
        // Nothing to do ignore the perturbations
        goto PERTURBATION_FINISHUP;
      }

      for (int yE = 0; yE < Ns; yE++) {
        for (int xE = 0; xE < Ns; xE++) {

          // GET COORDINATES RELATIVE TO POINT OF INTEREST //
          double DeltaX = (xS - xE) * afm;
          double DeltaY = (yS - yE) * afm;
          double Distance = std::sqrt(DeltaX * DeltaX + DeltaY * DeltaY);

          double rX = DeltaX / Distance;
          double rY = DeltaY / Distance;

          // CHECK THAT DISTANCES ARE RELEVANT FOR EVOLUTION //
          if (Distance >= CircleRadius) {
            // We are outside the causal circle
            continue;
          }

          // COMPUTE AMPLITUDE OF INITIAL PERTURBATIONS //
          double dT00In = (TIn->Get(0, 0, xE, yE) - T00InAvg) / T00InAvg;
          double dT0XIn = TIn->Get(0, 1, xE, yE) / T00InAvg;
          double dT0YIn = TIn->Get(0, 2, xE, yE) / T00InAvg;

          //////////////////////////////////////////////////////
          // COMPUTE ENERGY-MOMENTUM PERTURBATION PROPAGATORS //
          //////////////////////////////////////////////////////

          // GREENS-FUNCTIONS FOR ENERGY PERTURBATIONS //
          double G00_00 = 0.0;
          double G0X_00 = 0.0;
          double G0Y_00 = 0.0;
          double GXX_00 = 0.0;
          double GYY_00 = 0.0;
          double GZZ_00 = 0.0;
          double GXY_00 = 0.0;

          if (ENERGY_PERTURBATIONS) {

            // GREENS-FUNCTIONS IN TENSOR BASIS //
            double Gs, Gv, Gd, Gr;

            // COMPUTE GREENS FUCNTIONS //

            // Kinetic Theory
            if (EVOLUTION_MODE == 1) {

              Gs = GreensFunctions::EnergyPerturbations::KineticTheory::
                  CoordinateSpace::Gs(Distance, EvolutionTime, ScalingVarOut);
              Gv = GreensFunctions::EnergyPerturbations::KineticTheory::
                  CoordinateSpace::Gv(Distance, EvolutionTime, ScalingVarOut);
              Gd = GreensFunctions::EnergyPerturbations::KineticTheory::
                  CoordinateSpace::Gd(Distance, EvolutionTime, ScalingVarOut);
              Gr = GreensFunctions::EnergyPerturbations::KineticTheory::
                  CoordinateSpace::Gr(Distance, EvolutionTime, ScalingVarOut);

              // LowK limit of Kinetic Theory
            } else if (EVOLUTION_MODE == 2) {

              GreensFunctions::LowKArguments Lowk(
                  tOut, tIn, T00InAvg, T00BG, TXXInAvg, TXXBG,
                  TOutBG->GetCellData(4, xS, yS),
                  TOutBG->GetCellData(5, xS, yS));

              Gs = GreensFunctions::EnergyPerturbations::LowKLimit::
                  CoordinateSpace::Gs(Distance, EvolutionTime, Lowk);
              Gv = GreensFunctions::EnergyPerturbations::LowKLimit::
                  CoordinateSpace::Gv(Distance, EvolutionTime, Lowk);
              Gd = GreensFunctions::EnergyPerturbations::LowKLimit::
                  CoordinateSpace::Gd(Distance, EvolutionTime, Lowk);
              Gr = GreensFunctions::EnergyPerturbations::LowKLimit::
                  CoordinateSpace::Gr(Distance, EvolutionTime, Lowk);

              // Free Streaming
            } else if (EVOLUTION_MODE == 0) {

              Gs = GreensFunctions::EnergyPerturbations::FreeStreaming::
                  CoordinateSpace::Gs(Distance, EvolutionTime);
              Gv = GreensFunctions::EnergyPerturbations::FreeStreaming::
                  CoordinateSpace::Gv(Distance, EvolutionTime);
              Gd = GreensFunctions::EnergyPerturbations::FreeStreaming::
                  CoordinateSpace::Gd(Distance, EvolutionTime);
              Gr = GreensFunctions::EnergyPerturbations::FreeStreaming::
                  CoordinateSpace::Gr(Distance, EvolutionTime);

            } else {
              std::cerr << "#ERROR -- EVOLUTION MODE NOT SPECIFICED CORRECTLY"
                        << std::endl;
              exit(0);
            }

            // ENERGY-MOMENTUM TENSOR RESPONSE //
            if (Distance == 0) {

              G00_00 = Gs;

              G0X_00 = 0.0;
              G0Y_00 = 0.0;

              GXX_00 = Gd + 0.5 * Gr;
              GYY_00 = Gd + 0.5 * Gr;
              GZZ_00 = (G00_00 - GXX_00 - GYY_00);

              GXY_00 = 0.0;

            } else {

              G00_00 = Gs;

              G0X_00 = rX * Gv;
              G0Y_00 = rY * Gv;

              GXX_00 = Gd + rX * rX * Gr;
              GYY_00 = Gd + rY * rY * Gr;
              GZZ_00 = (G00_00 - GXX_00 - GYY_00);

              GXY_00 = rX * rY * Gr;
            }
          } // If energy perturbations

          // GREENS-FUNCTIONS FOR MOMENTUM PERTRUBATIONS  //
          double G00_0X = 0.0;
          double G0X_0X = 0.0;
          double G0Y_0X = 0.0;
          double GXX_0X = 0.0;
          double GYY_0X = 0.0;
          double GZZ_0X = 0.0;
          double GXY_0X = 0.0;
          double G00_0Y = 0.0;
          double G0X_0Y = 0.0;
          double G0Y_0Y = 0.0;
          double GXX_0Y = 0.0;
          double GYY_0Y = 0.0;
          double GZZ_0Y = 0.0;
          double GXY_0Y = 0.0;

          if (MOMENTUM_PERTURBATIONS) {

            // GREENS FUNCTIONS FOR MOMENTUM PERTURBATIONS //
            double Hv, Hd, Hr, Htd, Htm, Htr;

            // COMPUTE GREENS FUCNTIONS //
            if (EVOLUTION_MODE == 1) {

              Hv = GreensFunctions::MomentumPerturbations::KineticTheory::
                  CoordinateSpace::Hv(Distance, EvolutionTime, ScalingVarOut);
              Hd = GreensFunctions::MomentumPerturbations::KineticTheory::
                  CoordinateSpace::Hd(Distance, EvolutionTime, ScalingVarOut);
              Hr = GreensFunctions::MomentumPerturbations::KineticTheory::
                  CoordinateSpace::Hr(Distance, EvolutionTime, ScalingVarOut);
              Htd = GreensFunctions::MomentumPerturbations::KineticTheory::
                  CoordinateSpace::Htd(Distance, EvolutionTime, ScalingVarOut);
              Htm = GreensFunctions::MomentumPerturbations::KineticTheory::
                  CoordinateSpace::Htm(Distance, EvolutionTime, ScalingVarOut);
              Htr = GreensFunctions::MomentumPerturbations::KineticTheory::
                  CoordinateSpace::Htr(Distance, EvolutionTime, ScalingVarOut);

            }

            else if (EVOLUTION_MODE == 0) {

              Hv = GreensFunctions::MomentumPerturbations::FreeStreaming::
                  CoordinateSpace::Hv(Distance, EvolutionTime);
              Hd = GreensFunctions::MomentumPerturbations::FreeStreaming::
                  CoordinateSpace::Hd(Distance, EvolutionTime);
              Hr = GreensFunctions::MomentumPerturbations::FreeStreaming::
                  CoordinateSpace::Hr(Distance, EvolutionTime);
              Htd = GreensFunctions::MomentumPerturbations::FreeStreaming::
                  CoordinateSpace::Htd(Distance, EvolutionTime);
              Htm = GreensFunctions::MomentumPerturbations::FreeStreaming::
                  CoordinateSpace::Htm(Distance, EvolutionTime);
              Htr = GreensFunctions::MomentumPerturbations::FreeStreaming::
                  CoordinateSpace::Htr(Distance, EvolutionTime);

            }

            else {
              std::cerr << "#ERROR -- EVOLUTION MODE NOT SPECIFICED CORRECTLY"
                        << std::endl;
              exit(0);
            }

            // ENERGY-MOMENTUM TENSOR RESPONSE //
            if (Distance == 0) {

              G00_0X = 0.0;
              G00_0Y = 0.0;

              G0X_0X = Hd + 0.5 * Hr;
              G0Y_0X = 0.0;
              G0X_0Y = 0.0;
              G0Y_0Y = Hd + 0.5 * Hr;

              GXX_0X = 0.0;
              GYY_0X = 0.0;
              GXX_0Y = 0.0;
              GYY_0Y = 0.0;
              GZZ_0X = 0.0;
              GZZ_0Y = 0.0;

              GXY_0X = 0.0;
              GXY_0Y = 0.0;

            }

            else {

              G00_0X = rX * Hv;
              G00_0Y = rY * Hv;

              G0X_0X = Hd + rX * rX * Hr;
              G0Y_0X = rY * rX * Hr;
              G0X_0Y = rX * rY * Hr;
              G0Y_0Y = Hd + rY * rY * Hr;

              GXX_0X = rX * Htd + rX * Htm + rX * rX * rX * Htr;
              GYY_0X = rX * Htd + rY * rY * rX * Htr;
              GXX_0Y = rY * Htd + rX * rX * rY * Htr;
              GYY_0Y = rY * Htd + rY * Htm + rY * rY * rY * Htr;
              GZZ_0X = (G00_0X - GXX_0X - GYY_0X);
              GZZ_0Y = (G00_0Y - GXX_0Y - GYY_0Y);

              GXY_0X = 0.5 * (rY)*Htm + rX * rY * rX * Htr;
              GXY_0Y = 0.5 * (rX)*Htm + rX * rY * rY * Htr;
            }
          } // If momentum perturbations

          // COMPUTE CONTRIBUTION TO ENERGY-MOMENTUM TENSOR AT POINT OF
          // INTEREST //
          T00Pert += (afm * afm) *
                     (+G00_00 * dT00In - G00_0X * dT0XIn - G00_0Y * dT0YIn) *
                     T00BG;
          TXXPert += (afm * afm) *
                     (+GXX_00 * dT00In - GXX_0X * dT0XIn - GXX_0Y * dT0YIn) *
                     T00BG;
          TYYPert += (afm * afm) *
                     (+GYY_00 * dT00In - GYY_0X * dT0XIn - GYY_0Y * dT0YIn) *
                     T00BG;
          TZZPert += (afm * afm) *
                     (+GZZ_00 * dT00In - GZZ_0X * dT0XIn - GZZ_0Y * dT0YIn) *
                     T00BG;

          T0XPert += (afm * afm) *
                     (-G0X_00 * dT00In + G0X_0X * dT0XIn + G0X_0Y * dT0YIn) *
                     T00BG;
          T0YPert += (afm * afm) *
                     (-G0Y_00 * dT00In + G0Y_0X * dT0XIn + G0Y_0Y * dT0YIn) *
                     T00BG;
          TXYPert += (afm * afm) *
                     (-GXY_00 * dT00In + GXY_0X * dT0XIn + GXY_0Y * dT0YIn) *
                     T00BG;
        } // Loop over the  x-coordinate of causal patch
      }   // Loop over the y-coordiante of causal patch

    PERTURBATION_FINISHUP:

      TOutFull->Set(xS, yS, T00BG + T00Pert, TXXBG + TXXPert, TYYBG + TYYPert,
                    TZZBG + TZZPert, T0XPert, T0YPert, 0.0, TXYPert, 0.0, 0.0);

    } // Loop over the x-coordinate of grid
  }   // Loop over the y-coordinate of grid

  // COMMANDLINE OUTPUT //
  std::cerr<< "#DONE" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

void Setup() {

  using namespace KoMPoSTParameters;

  // SETUP GREENS FUNCTIONS FOR ENERGY-MOMENTUM PERTURBATIONS //
  const int NumberOfPoints = 256;
  GreensFunctions::Setup(Sigma, NumberOfPoints, ENERGY_PERTURBATIONS,
                         MOMENTUM_PERTURBATIONS);

  // (optional) OUTPUT COORDINATE SPACE RESPONSE FUNCTIONS //
  GreensFunctions::Output(ENERGY_PERTURBATIONS, MOMENTUM_PERTURBATIONS);
}

////////////////////////////////////////////////////////////////////////////////

void Run(EnergyMomentumTensorMap *TIn, EnergyMomentumTensorMap *TOutBG,
         EnergyMomentumTensorMap *TOutFull) {

  using namespace KoMPoSTParameters;

  double EvolutionTime = TOutFull->tau - TIn->tau;
  double SigmaBG = 0.;

  // The Scaler is responsible for computing the scaling variable for a given
  // eta/s. The parameters EtaOverS and EtaOverSTemperature scale are passed
  // from the common block KoMPoSTParameters.
  ScalingVariable Scaler(EtaOverS, EtaOverSTemperatureScale);

  // Some error checking to make sure that the system can evolve the low k limit
  if (KoMPoSTParameters::Regulator == "TwoPass" && EVOLUTION_MODE == 2) {
    std::cerr << "*** KoMPoST::Run *** Regulator TwoPass is not appropriate for "
                 "the low k limit, selected by EVOLUTION_MODE == 2. Aborting!"
              << std::endl;
    exit(0);
  }

  // If the two pass regulator is used, then we run through the program twice
  // The first time is just to estimate the size of the perturbations.  In the
  // second pass we regulate the perturbations
  int npass = 1;
  if (KoMPoSTParameters::Regulator == "TwoPass") {
    npass = 2;
  }
  bool IsFirstPass = true;

  // Start the run
  for (int ipass = 0; ipass < npass; ipass++) {

    if (EVOLUTION_MODE != 2) {
      SigmaBG = EvolutionTime / sqrt(2.);
    } else {
      SigmaBG = EvolutionTime / sqrt(2.);
    }

    // Evolve the backround, compute the scaling variable, K etc
    ComputeBackground(IsFirstPass, TIn, TOutBG, Scaler, SigmaBG,
                      EVOLUTION_MODE);

    // The peturbations are propagated from anywhere within a circle of radius
    // CicleRadius.  Typically this is just the causal circle (plus a little
    // bit) But for the low k limit (EVOLUTION_MODE==2) we take something
    // different.
    double CircleRadius = 0.;
    if (EVOLUTION_MODE != 2) {
      // Causal circle
      CircleRadius = EvolutionTime * (1. + 3. * Sigma);
    } else {
      CircleRadius = 4 * SigmaBG; // Take a circle four times the BG
    }

    // Evolve the petrubations. The scaling variable and K are
    // stored in the CellData structure of TOutFull.
    ComputePerturbations(TIn, TOutBG, TOutFull, ENERGY_PERTURBATIONS,
                         MOMENTUM_PERTURBATIONS, EVOLUTION_MODE, CircleRadius);

    // Regulate the perturbations so that the inversion problem is well posed.
    if (KoMPoSTParameters::Regulator == "TwoPass") {
      if (IsFirstPass) {
        std::cerr << "#Preparing for second pass ... " << std::endl;
        PrepareKoMPoSTEstimate(TOutBG, TOutFull);
        IsFirstPass = false;
      }
    } else if (KoMPoSTParameters::Regulator == "KoMPoSTAddition") {
      std::cerr<< "#Regulation with KoMPoSTAddition ... " << std::endl;
      RegulateKoMPoSTAddition(TOutBG, TOutFull);
    } else if (KoMPoSTParameters::Regulator == "NoRegulator") {
      std::cerr<< "#Returning unregulated output ... " << std::endl;
    } else {
      std::cerr << "#KoMPoSTParameters::Regulator string does not match any of the "
                   "expected choices TwoPass/KoMPoSTAddition/NoRegulator! Aborting"
                << std::endl;
    }
  }
}

} // Namepsace KoMPoST
