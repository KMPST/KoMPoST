#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <complex>
#include <gsl/gsl_math.h>


using namespace std;

#include "../em_rates_parameters.h"
#include "ekt_photon.h"
#include "ekt_yields_photon.h"

namespace Photons{
    
    void GetPhoton(EnergyMomentumTensorMap *TOutFull) {
        
        std::cout << "#PHOTON PRODUCTION STARTED" << std::endl;

        // Call relevant constants coming from KoMPoST evolution
        double NuEff=KoMPoSTParameters::NuG;
        double etaS=KoMPoSTParameters::EtaOverS;
        
        // CIdeal: Normalization constant for scaling curve; Tconf: Minimal temperature used for photon production
        double CIdeal=EMParameters::CIdeal;
        double Tconf=EMParameters::TConf;

        // Get hydrodynamization time from final energy-momentum tensor
        double tauHydro=TOutFull->tau;

        // Set momentum grid in pT
        double *PTArr = new double[NK];
        for (int ik = 0; ik <  NK; ik++){  PTArr[ik]=Kmin+ ik*dK;}

        // Import constants, evolution parameters and compute scaling curve from EKT simulation
        EKTPhoton Photon;

        // Compute pre-equilibrium photon yield
        EKTYieldsPhoton::ekt_yield_event(TOutFull,NuEff,PTArr,tauHydro,etaS,CIdeal,Tconf,&Photon);

        std::cout << "#PHOTON PRODUCTION DONE" << std::endl;


    }

}

