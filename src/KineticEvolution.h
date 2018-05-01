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
#ifndef KINETICEVOLUTION_H
#define KINETICEVOLUTION_H

class EnergyMomentumTensorMap;
class ScalingVariable;

namespace KoMPoST {

//! Setup Green Functions based on the parameters in KoMPoSTParameters
void Setup();

//! Run the KoMPoST evolution based on  the parameters in KoMPoSTParameters
void Run(EnergyMomentumTensorMap *TIn, EnergyMomentumTensorMap *TOutBG,
         EnergyMomentumTensorMap *TOutFull);

void PrepareKoMPoSTEstimate(EnergyMomentumTensorMap *TOutBG,
                        EnergyMomentumTensorMap *TOutFull);

void RegulateKoMPoSTAddition(EnergyMomentumTensorMap *TOutBG,
                         EnergyMomentumTensorMap *TOutFull);

void ComputeBackground(bool IsFirstPass, EnergyMomentumTensorMap *TIn,
                       EnergyMomentumTensorMap *TOutBG,
                       const ScalingVariable &ScalerIn, double SigmaBG,
                       int EVOLUTION_MODE);

void ComputePerturbations(EnergyMomentumTensorMap *TIn,
                          EnergyMomentumTensorMap *TOutBG,
                          EnergyMomentumTensorMap *TOutFull,
                          int ENERGY_PERTURBATIONS, int MOMENTUM_PERTURBATIONS,
                          int EVOLUTION_MODE, double CircleRadius);
}

#endif
