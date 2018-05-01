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
    double *HkVValues;
    double *HkDValues;
    double *HkKValues;
    double *HkTDValues;
    double *HkTMValues;
    double *HkTKValues;
    
    
    // GSL INTERPOLATION OBJECTS //
    gsl_interp_accel *HkVAcc,*HkDAcc,*HkKAcc,*HkTDAcc,*HkTMAcc,*HkTKAcc;
    gsl_spline *HkVInt,*HkDInt,*HkKInt,*HkTDInt,*HkTMInt,*HkTKInt;
    
    // ASYMPTOTIC COEFFICIENTS //
    double Zv,av,sv;
    double Zd,ad,sd;
    double Zk,ak,sk;
    double Ztd,atd,std;
    double Ztm,atm,stm;
    double Ztk,atk,stk;
    
    // DEFINITION OF RESPONSE FUNCTIONS BASED ON GSL INTERPOLATORS //
    
    double HKv(double k){
        
        if(k<kMin){
            std::cerr << "#ERROR -- ACCESSING k=" << k << " WHICH IS OUT OF RANGE" << std::endl;
            exit(0);
        }
        else if(k<kMax){
            return gsl_spline_eval(HkVInt,k,HkVAcc);
        }
        else{
            return 2.0*Zv*BesselJ1(av*k+sv);
        }
    }
    
    double HKd(double k){
        
        if(k<kMin){
            std::cerr << "#ERROR -- ACCESSING k=" << k << " WHICH IS OUT OF RANGE" << std::endl;
            exit(0);
        }
        else if(k<kMax){
            return gsl_spline_eval(HkDInt,k,HkDAcc);
        }
        else{
            return 2.0*Zd*BesselJ1(ad*k+sd)/(ad*k+sd);
        }
        
    }
    
    double HKk(double k){
        
        if(k<kMin){
            std::cerr << "#ERROR -- ACCESSING k=" << k << " WHICH IS OUT OF RANGE" << std::endl;
            exit(0);
        }
        else if(k<kMax){
            return gsl_spline_eval(HkKInt,k,HkKAcc);
        }
        else{
            return -2.0*Zk*BesselJ2(ak*k+sk);
        }
    }
    
    double HKtd(double k){
        
        if(k<kMin){
            std::cerr << "#ERROR -- ACCESSING k=" << k << " WHICH IS OUT OF RANGE" << std::endl;
            exit(0);
        }
        else if(k<kMax){
            return gsl_spline_eval(HkTDInt,k,HkTDAcc);
        }
        else{
            return 2.0*Ztd*BesselJ2(atd*k+std)/(atd*k+std);
        }
    }
    
    double HKtm(double k){
        
        if(k<kMin){
            std::cerr << "#ERROR -- ACCESSING k=" << k << " WHICH IS OUT OF RANGE" << std::endl;
            exit(0);
        }
        else if(k<kMax){
            return gsl_spline_eval(HkTMInt,k,HkTMAcc);
        }
        else{
            return 4.0*Ztm*BesselJ2(atm*k+stm)/(atm*k+stm);
        }
    }
    
    double HKtk(double k){
        
        if(k<kMin){
            std::cerr << "#ERROR -- ACCESSING k=" << k << " WHICH IS OUT OF RANGE" << std::endl;
            exit(0);
        }
        else if(k<kMax){
            return gsl_spline_eval(HkTKInt,k,HkTKAcc);
        }
        else{
            return -2.0*Ztk*BesselJ3(atk*k+stk);
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
        double ks,HkV,HkD,HkK,HkTD,HkTM,HkTK;
        
        int InputCount=0;
        
        // GET ASMYPTOTIC COEFFICIENTS //
        if(InStream.good()){
            
            // GET FIRST LINE //
            getline(InStream,InputLine);
            
            InputLine.erase(InputLine.find("#"),1);
            std::stringstream InputValues(InputLine);
            
            // SET ASYMPTOTIC COEFFICIENTS //
            InputValues >> Zv; InputValues >> av; InputValues >> sv;
            InputValues >> Zd; InputValues >> ad; InputValues >> sd;
            InputValues >> Zk; InputValues >> ak; InputValues >> sk;
            
            InputValues >> Ztd; InputValues >> atd; InputValues >> std;
            InputValues >> Ztm; InputValues >> atm; InputValues >> stm;
            InputValues >> Ztk; InputValues >> atk; InputValues >> stk;

           // Print Free streaming tail fit parameters 
            //    std::cerr << av << " " << ad << " " << ak  << " " << atd << " " << atm << " " << atk << std::endl;
            
            
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
                InputValues >> HkV; InputValues >> HkD; InputValues >> HkK;
                InputValues >> HkTD; InputValues >> HkTM; InputValues >> HkTK;
                
                // SET VALUES //
                kValues[InputCount]=ks;
                
                HkVValues[InputCount]=HkV;
                HkDValues[InputCount]=HkD;
                HkKValues[InputCount]=HkK;
                
                HkTDValues[InputCount]=HkTD;
                HkTMValues[InputCount]=HkTM;
                HkTKValues[InputCount]=HkTK;
                
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
        HkVValues=new double[NumberOfMomenta];
        HkDValues=new double[NumberOfMomenta];
        HkKValues=new double[NumberOfMomenta];
        
        HkTDValues=new double[NumberOfMomenta];
        HkTMValues=new double[NumberOfMomenta];
        HkTKValues=new double[NumberOfMomenta];
        
        // GET EKT GREENS FUNCTIONS FROM INPUT FILE //
        GetGreensFunction(fname,NumberOfMomenta);
        
        // SET BOUNDARY //
        kMin=kValues[0]; kMax=kValues[NumberOfMomenta-1];
        
        // SETUP GSL INTERPOLATION //
        HkVAcc = gsl_interp_accel_alloc ();
        HkDAcc = gsl_interp_accel_alloc ();
        HkKAcc = gsl_interp_accel_alloc ();
        HkTDAcc = gsl_interp_accel_alloc ();
        HkTMAcc = gsl_interp_accel_alloc ();
        HkTKAcc = gsl_interp_accel_alloc ();
        
        HkVInt=gsl_spline_alloc(gsl_interp_akima,NumberOfMomenta);
        HkDInt=gsl_spline_alloc(gsl_interp_akima,NumberOfMomenta);
        HkKInt=gsl_spline_alloc(gsl_interp_akima,NumberOfMomenta);
        HkTDInt=gsl_spline_alloc(gsl_interp_akima,NumberOfMomenta);
        HkTMInt=gsl_spline_alloc(gsl_interp_akima,NumberOfMomenta);
        HkTKInt=gsl_spline_alloc(gsl_interp_akima,NumberOfMomenta);
        
        gsl_spline_init(HkVInt,kValues,HkVValues,NumberOfMomenta);
        gsl_spline_init(HkDInt,kValues,HkDValues,NumberOfMomenta);
        gsl_spline_init(HkKInt,kValues,HkKValues,NumberOfMomenta);
        gsl_spline_init(HkTDInt,kValues,HkTDValues,NumberOfMomenta);
        gsl_spline_init(HkTMInt,kValues,HkTMValues,NumberOfMomenta);
        gsl_spline_init(HkTKInt,kValues,HkTKValues,NumberOfMomenta);

        
    }
    
    
    void Reset(){
        
        // GSL CLEAN-UP //
        gsl_spline_free (HkVInt);
        gsl_spline_free (HkDInt);
        gsl_spline_free (HkKInt);
        gsl_spline_free (HkTDInt);
        gsl_spline_free (HkTMInt);
        gsl_spline_free (HkTKInt);

        
        gsl_interp_accel_free (HkVAcc);
        gsl_interp_accel_free (HkDAcc);
        gsl_interp_accel_free (HkKAcc);
        gsl_interp_accel_free (HkTDAcc);
        gsl_interp_accel_free (HkTMAcc);
        gsl_interp_accel_free (HkTKAcc);
        
        
        // CLEAN-UP //
        delete[] kValues;
        delete[] HkVValues;
        delete[] HkDValues;
        delete[] HkKValues;
        delete[] HkTDValues;
        delete[] HkTMValues;
        delete[] HkTKValues;
        
    }
    
    
    ///////////////////
    // CREATE OUTPUT //
    ///////////////////
    
    void Output(std::string fname,int Nk,double k0,double k1){
        
        std::ofstream OutStream;
        OutStream.open(fname.c_str());
        
        for(int kIndex=0;kIndex<1000;kIndex++){
            
            double kVal=k0+(kIndex+0.5)*(k1-k0)/(Nk);
            
            OutStream << kVal << " " << HKv(kVal) << " " << HKd(kVal) << " " << HKk(kVal)  << " " << HKtd(kVal) << " " << HKtm(kVal) << " " << HKtk(kVal) << std::endl;
            
        }
        
        OutStream.close();
        
        
    }
    
    
}
