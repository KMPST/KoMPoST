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
    // KINETIC THEORY //
    ////////////////////
    
    // MIN/MAX MOMENTUM //
    double kMin; double kMax;
    
    // MOMENTUM VALUES //
    double *kValues;
    
    // GREENS FUNCTION VALUES //
    double *GkSValues;
    double *GkVValues;
    double *GkDValues;
    double *GkKValues;
    
    // GSL INTERPOLATION OBJECTS //
    gsl_interp_accel *GkSAcc,*GkVAcc,*GkDAcc,*GkKAcc;
    gsl_spline *GkSInt,*GkVInt,*GkDInt,*GkKInt;
    
    // ASYMPTOTIC COEFFICIENTS //
    double Zs,as,ss;
    double Zv,av,sv;
    double Zd,ad,sd;
    double Zk,ak,sk;
    
    // DEFINITION OF RESPONSE FUNCTIONS BASED ON GSL INTERPOLATORS //
    double GKs(double k){
        
        if(k<kMin){
            std::cerr << "#ERROR -- ACCESSING k=" << k << " WHICH IS OUT OF RANGE" << std::endl;
            exit(0);
        }
        else if(k<kMax){
            return gsl_spline_eval(GkSInt,k,GkSAcc);
        }
        else{
            return Zs*BesselJ0(as*k+ss);
        }
    }
    
    double GKv(double k){
        
        if(k<kMin){
            std::cerr << "#ERROR -- ACCESSING k=" << k << " WHICH IS OUT OF RANGE" << std::endl;
            exit(0);
        }
        else if(k<kMax){
            return gsl_spline_eval(GkVInt,k,GkVAcc);
        }
        else{
            return Zv*BesselJ1(av*k+sv);
        }
    }
    
    double GKd(double k){
        
        if(k<kMin){
            std::cerr << "#ERROR -- ACCESSING k=" << k << " WHICH IS OUT OF RANGE" << std::endl;
            exit(0);
        }
        else if(k<kMax){
            return gsl_spline_eval(GkDInt,k,GkDAcc);
        }
        else{
            return Zd*BesselJ1(ad*k+sd)/(ad*k+sd);
        }
    }
    
    double GKk(double k){
        
        if(k<kMin){
            std::cerr << "#ERROR -- ACCESSING k=" << k << " WHICH IS OUT OF RANGE" << std::endl;
            exit(0);
        }
        else if(k<kMax){
            return gsl_spline_eval(GkKInt,k,GkKAcc);
        }
        else{
            return -Zk*BesselJ2(ak*k+sk);
        }
    }
    
    
    /////////////////////////////
    // GET EKT GREENS FUNCTION //
    /////////////////////////////
    
    void GetGreensFunction(std::string fname,int NumberOfMomenta){
        
            std::cerr << "#IMPORTING " << fname << std::endl;
        
        // OPEN INPUT FILE //
        std::ifstream InStream;
        InStream.open(fname.c_str());
        
        std::string InputLine;
        
        // ENTRIES IN FILE //
        double ks,GkS,GkV,GkD,GkK;
        
        int InputCount=0;
        
        // GET ASMYPTOTIC COEFFICIENTS //
        if(InStream.good()){
            
            // GET FIRST LINE //
            getline(InStream,InputLine);
            
            InputLine.erase(InputLine.find("#"),1);
            std::stringstream InputValues(InputLine);
            
            // SET ASYMPTOTIC COEFFICIENTS //
            InputValues >> Zs; InputValues >> as; InputValues >> ss;
            InputValues >> Zv; InputValues >> av; InputValues >> sv;
            InputValues >> Zd; InputValues >> ad; InputValues >> sd;
            InputValues >> Zk; InputValues >> ak; InputValues >> sk;
            
           // Print free streaming tail fit parameters 
           //    std::cerr << as << " " << av << " " << ad << " " << ak  << " " << std::endl;
            
            
        }
        
        // GET DATA LINE BY LINE //
        while(InStream.good()){
            
            // GET LINE //
            getline(InStream,InputLine);
            
            if(!(InputLine.empty())){
                
                std::stringstream InputValues(InputLine);
                
                // GET INPUT MOMENTUM AND TIMES //
                InputValues >> ks;
                
                // GET GREENS FUNCTION //
                InputValues >> GkS; InputValues >> GkV;
                InputValues >> GkD; InputValues >> GkK;
                
                // SET VALUES //
                kValues[InputCount]=ks;
                
                GkSValues[InputCount]=GkS;
                GkVValues[InputCount]=GkV;
                GkDValues[InputCount]=GkD;
                GkKValues[InputCount]=GkK;
                
                // COUNTER //
                InputCount++;
                
                
            }
            
        }
        
        // CHECK //
        if(InputCount!=NumberOfMomenta){
            std::cerr << "#ERROR -- COULD NOT LOAD INPUT FILE " << fname << " CORRECTLY" << std::endl;
            exit(0);
        }
        
        // CLOSE //
        InStream.close();
        
    }
    
    ///////////////////////////
    // SETUP GREENS FUNCTION //
    ///////////////////////////
    
    void Setup(std::string fname,int NumberOfMomenta){
        
        // MOMENTUM VALUES //
        kValues=new double[NumberOfMomenta];
        
        // GREENS FUNCTION VALUES //
        GkSValues=new double[NumberOfMomenta];
        GkVValues=new double[NumberOfMomenta];
        GkDValues=new double[NumberOfMomenta];
        GkKValues=new double[NumberOfMomenta];
        
        // GET EKT GREENS FUNCTIONS FROM INPUT FILE //
        GetGreensFunction(fname,NumberOfMomenta);
        
        // SET BOUNDARY //
        kMin=kValues[0]; kMax=kValues[NumberOfMomenta-1];
        
        // SETUP GSL INTERPOLATION //
        GkSAcc = gsl_interp_accel_alloc ();
        GkVAcc = gsl_interp_accel_alloc ();
        GkDAcc = gsl_interp_accel_alloc ();
        GkKAcc = gsl_interp_accel_alloc ();
        
        GkSInt=gsl_spline_alloc(gsl_interp_akima,NumberOfMomenta);
        GkVInt=gsl_spline_alloc(gsl_interp_akima,NumberOfMomenta);
        GkDInt=gsl_spline_alloc(gsl_interp_akima,NumberOfMomenta);
        GkKInt=gsl_spline_alloc(gsl_interp_akima,NumberOfMomenta);
        
        gsl_spline_init(GkSInt,kValues,GkSValues,NumberOfMomenta);
        gsl_spline_init(GkVInt,kValues,GkVValues,NumberOfMomenta);
        gsl_spline_init(GkDInt,kValues,GkDValues,NumberOfMomenta);
        gsl_spline_init(GkKInt,kValues,GkKValues,NumberOfMomenta);
        
    }
    
    
    void Reset(){
        
        // GSL CLEAN-UP //
        gsl_spline_free (GkSInt);
        gsl_spline_free (GkVInt);
        gsl_spline_free (GkDInt);
        gsl_spline_free (GkKInt);
        
        gsl_interp_accel_free (GkSAcc);
        gsl_interp_accel_free (GkVAcc);
        gsl_interp_accel_free (GkDAcc);
        gsl_interp_accel_free (GkKAcc);
        
        
        // CLEAN-UP //
        delete[] kValues;
        delete[] GkSValues;
        delete[] GkVValues;
        delete[] GkDValues;
        delete[] GkKValues;
        
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
