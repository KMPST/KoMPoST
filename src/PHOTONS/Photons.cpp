#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <complex>
#include <gsl/gsl_math.h>


using namespace std;

#include "parameters.h"
#include "ektphoton.h"
#include "ekt_yields.h"

namespace Photons{
    
    void GetPhoton(EnergyMomentumTensorMap *TOutFull) {
        
        std::cout << "#PHOTON PRODUCTION STARTED" << std::endl;

        // Call relevant constants coming from KoMPoST evolution
        double NuEff=KoMPoSTParameters::NuG;
        double etaS=KoMPoSTParameters::EtaOverS;
        
        // CIdeal: Normalization constant for scaling curve; Tconf: Miniml temperature used for photon production
        double CIdeal=PhotonParameters::CIdeal;
        double Tconf=PhotonParameters::TConf;

        // Get hydrodynamization time from final energy-momentum tensor
        double tauHydro=TOutFull->tau;

        // Set momentum grid in pT
        double *PTArr = new double[NK];
        for (int ik = 0; ik <  NK; ik++){  PTArr[ik]=Kmin+ ik*dK;}

        // Import constants, evolution parameters and compute scaling curve from EKT simulation
        EKTPhoton Photon;

        // Compute pre-equilibrium photon yield
        EKTYields::ekt_yield_event(TOutFull,NuEff,PTArr,tauHydro,etaS,CIdeal,Tconf,&Photon);

        std::cout << "#PHOTON PRODUCTION DONE" << std::endl;


    }

}

