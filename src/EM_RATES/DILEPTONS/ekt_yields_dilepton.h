#ifndef EKTYieldsDILEPTON_H
#define EKTYieldsDILEPTON_H

#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>

namespace EKTYieldsDilepton{
    
    // Modify volume element in case the cell is on the boundary of the grid
    double AreaWeight(int ix, int iy){
        if(ix==EventInput::xSTART || ix==EventInput::xEND-1){
            if(iy==EventInput::ySTART || iy==EventInput::yEND-1){return 1./4.;}
            else{return 1./2.;}
        }
        else{
            return 1.;
        }
    }

    // Compute dilepton spectrum for a given event for eta/s
	int ekt_yield_event(EnergyMomentumTensorMap *TOutFull,double NuEff,double * Marr,double thydro, double eta_o_s_in, double Tconf, EKTDilepton * Dilepton){   
        
        double * ektRate = new double[NK];
        double * ektRate_t = new double[NK];
        double wtilde_ev=0;
        double Ncells_sampled=0;
        
        // Clean the arrays
        for (int ik = 0; ik < NK; ik++){ektRate[ik]=0; ektRate_t[ik]=0;}

        // Set transverse area
        double dAT=EventInput::afm*EventInput::afm;

        // Perform integration over transverse plane
        for (int xS = EventInput::xSTART; xS <= EventInput::xEND; xS++) {
            for (int yS = EventInput::ySTART; yS <= EventInput::yEND; yS++) {
                
                // Extract energy density from KoMPoST evolution in GeV^4
                double Ed = TOutFull->Get(0,0,xS,yS);

                // Compute the temperature
                double Temp = pow(30./(M_PI*M_PI*NuEff)*Ed,1./4.);

                // In case the temperature is larger than the minimal temperature, compute pre-equilibrium dileptons
                if(Temp>Tconf)
                {
                    Dilepton->compute_dileptons(thydro*fmtoGeVm1, Temp,eta_o_s_in, Marr,  ektRate_t ,NK);

                    wtilde_ev += Dilepton->find_wtilde(thydro*fmtoGeVm1,Temp, eta_o_s_in);
                    Ncells_sampled += 1.0;

                    for (int ik = 0; ik < NK; ik++){ektRate[ik] += dAT * AreaWeight(xS,yS) * ektRate_t[ik]/GeVm2tofm2;}  
                } 

            }
        }
        

        // Create output
        ofstream data_plot;
        ostringstream data_plot_ll_name;
        
        data_plot_ll_name << "EKTDilepton_Nf_"<< nfamy << "_TauHydro_"<< thydro << "_eta_o_s_"<< eta_o_s_in << "_Tmin_"<< Tconf << ".txt";
        data_plot.open(data_plot_ll_name.str().c_str());
        data_plot << "#TauHydro = " << thydro << "fm" << std::endl;
        data_plot << "#[1]M [GeV] [2]Rate [GeV^{-2}]" << std::endl;
        for (int ik = 0; ik < NK; ik++){ data_plot << Marr[ik] << "\t"<< ektRate[ik] << std::endl;}
        data_plot.close();
        cout << "EKT Yields for Dileptons done with <wtilde> = " <<  wtilde_ev/Ncells_sampled << std::endl;
        return 0;
	}

}
	
#endif