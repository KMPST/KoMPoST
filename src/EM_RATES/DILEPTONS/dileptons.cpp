#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <complex>
#include <gsl/gsl_math.h>


using namespace std;

#include "../em_rates_parameters.h"
#include "ekt_dilepton.h"
#include "ekt_yields_dilepton.h"

namespace Dileptons{
    
    void GetDilepton(EnergyMomentumTensorMap *TOutFull) {
        
        std::cout << "#DILEPTON PRODUCTION STARTED" << std::endl;

        // Call relevant constants coming from KoMPoST evolution
        double NuEff=KoMPoSTParameters::NuG;
        double etaS=KoMPoSTParameters::EtaOverS;
        
        // Tconf: Minimal temperature used for photon production
        double Tconf=EMParameters::TConf;

        // Get hydrodynamization time from final energy-momentum tensor
        double tauHydro=TOutFull->tau;

        // Set mass grid in M
        double *MArr = new double[NK];
        for (int ik = 0; ik <  NK; ik++){  MArr[ik]=Kmin+ ik*dK;}

        // Import constants, evolution parameters and compute scaling curve from EKT simulation
        EKTDilepton Dilepton;

        // Compute pre-equilibrium dilepton yield
        EKTYieldsDilepton::ekt_yield_event(TOutFull,NuEff,MArr,tauHydro,etaS,Tconf,&Dilepton);

        std::cout << "#DILEPTON PRODUCTION DONE" << std::endl;


    }

}

