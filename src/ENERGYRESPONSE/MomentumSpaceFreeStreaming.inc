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
/////////////////////////////////
// MOMENTUM SPACE EXPRESSIONS  //
/////////////////////////////////

namespace MomentumSpace{
    
    ////////////////////
    // FREE STREAMING //
    ////////////////////
    
    //  G00=GKs //
    double GKs(double ks){
        return BesselJ0(ks);
    }
    
    // G0i=GKv k^{i}/k //
    double GKv(double ks){
        return BesselJ1(ks);
    }
    
    // Gij=GKd \delta^{ij} + GKk k^{i}k^{j}/k^2 //
    double GKd(double ks){
        return BesselJ1(ks)/(ks);
    }
    
    double GKk(double ks){
        return -BesselJ2(ks);
    }
    
    
    
    ///////////////////
    // CREATE OUTPUT //
    ///////////////////
    
    void Output(std::string fname,int Nk,double k0,double k1){
        
        std::ofstream OutStream;
        OutStream.open(fname.c_str());
        
        for(int kIndex=0;kIndex<1000;kIndex++){
            
            double kVal=k0+(kIndex+0.5)*(k1-k0)/(Nk);
            
            OutStream << kVal << " " << GKs(kVal)  << " " << GKv(kVal) << " " << GKd(kVal) << " " << GKk(kVal) << std::endl;
            
        }
        
        OutStream.close();
        
        
    }
    
    
}
