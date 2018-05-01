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

/////////////////////////////////////////////////////////////////////
// COMPUTE BESSEL TRANSFORMS OF MOMENTUM SPACE RESPONSE FUNCTIONS  //
/////////////////////////////////////////////////////////////////////

namespace BesselTransform{
    
    // NUMBER OF POINTS USED IN QAG INTEGRATION //
    static const int NumberOfIntegrationPoints=64000;
    
    // DEFINITIONS OF MOMENTUM SPACE GREENS FUNCTIONS //
    using namespace MomentumSpace;
    
    // DEFINITION OF INTEGRAND FUNCTIONS //
    double GsIntegrand(double k,void *params){
        
        GSLVariables Variables=*(GSLVariables *) params;
        double r = Variables.r; double s=Variables.s;
        
        return 1.0/(2.0*M_PI)*k*BesselJ0(k*r)*GKs(k)*Smearing::KernelK(k,s);
    }
    
    double GvIntegrand(double k,void *params){
    
        GSLVariables Variables=*(GSLVariables *) params;
        double r = Variables.r; double s=Variables.s;
        
        return 1.0/(2.0*M_PI)*k*BesselJ1(k*r)*GKv(k)*Smearing::KernelK(k,s);
    }
    
    double GdIntegrand(double k,void *params){
        
        GSLVariables Variables=*(GSLVariables *) params;
        double r = Variables.r; double s=Variables.s;
        
        return 1.0/(2.0*M_PI)*k*(BesselJ0(k*r)*GKd(k)+BesselJ1(k*r)/(k*r)*GKk(k))*Smearing::KernelK(k,s);
    }
    
    double GrIntegrand(double k,void *params){

        GSLVariables Variables=*(GSLVariables *) params;
        double r = Variables.r; double s=Variables.s;
        
        return -1.0/(2.0*M_PI)*k*BesselJ2(k*r)*GKk(k)*Smearing::KernelK(k,s);
    }
    
    
    
    // COMPUTE INTEGRALS //
    double ComputeGs(double r,double s){
        GSLVariables Variables; Variables.r=r; Variables.s=s;
        EVALUATE_GSL_INTEGRAL(GsIntegrand,Variables);
    }
    
    double ComputeGv(double r,double s){
        GSLVariables Variables; Variables.r=r; Variables.s=s;
        EVALUATE_GSL_INTEGRAL(GvIntegrand,Variables);
    }
    
    double ComputeGd(double r,double s){
        GSLVariables Variables; Variables.r=r; Variables.s=s;
        EVALUATE_GSL_INTEGRAL(GdIntegrand,Variables);
    }
    
    double ComputeGr(double r,double s){
        GSLVariables Variables; Variables.r=r; Variables.s=s;
        EVALUATE_GSL_INTEGRAL(GrIntegrand,Variables);
    }
    
    
}
