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
    
    
    //  H00k=HKv k^{k}/k //
    double HKv(double ks){
        return 2.0*BesselJ1(ks);
    }
    
    // H0i,j=HKd \delta^{ij} + HKk k^{i}k^{j}/k^2 //
    double HKd(double ks){
        return 2.0*BesselJ1(ks)/(ks);
    }
    
    double HKk(double ks){
        return -2.0*BesselJ2(ks);
    }
    
    
    // Hij,k=HKtd\delta^{ij} k^{k}/k + HKtm 1/2(k^{i}/k \delta^{jk} + k^{j}/k \delta^{jk}) + HKtk k^{i}k^{j}k^{k}/k^3
    double HKtd(double ks){
        return 2.0*BesselJ2(ks)/(ks);
    }
    
    double HKtm(double ks){
        return 4.0*BesselJ2(ks)/(ks);
    }
    double HKtk(double ks){
        return -2.0*BesselJ3(ks);
    }
    
    
    
    ///////////////////
    // CREATE OUTPUT //
    ///////////////////
    
    void Output(std::string fname,int Nk,double k0,double k1){
        
        std::ofstream OutStream;
        OutStream.open(fname.c_str());
        
        for(int kIndex=0;kIndex<1000;kIndex++){
            
            double kVal=k0+(kIndex+0.5)*(k1-k0)/(Nk);
            
            OutStream << kVal  << " " << HKv(kVal) << " " << HKd(kVal) << " " << HKk(kVal) << " " << HKtd(kVal) << " " << HKtm(kVal) << " " << HKtk(kVal) << std::endl;
            
        }
        
        OutStream.close();
        
        
    }
    
    
}
