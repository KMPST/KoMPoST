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
    double HvIntegrand(double k,void *params){
        
        GSLVariables Variables=*(GSLVariables *) params;
        double r = Variables.r; double s=Variables.s;
        
        return 1.0/(2.0*M_PI)*k*BesselJ1(k*r)*HKv(k)*Smearing::KernelK(k,s);
    }
    
    double HdIntegrand(double k,void *params){
        
        GSLVariables Variables=*(GSLVariables *) params;
        double r = Variables.r; double s=Variables.s;
        
        return 1.0/(2.0*M_PI)*k*(BesselJ0(k*r)*HKd(k)+BesselJ1(k*r)/(k*r)*HKk(k))*Smearing::KernelK(k,s);
    }
    
    double HrIntegrand(double k,void *params){
        
        GSLVariables Variables=*(GSLVariables *) params;
        double r = Variables.r; double s=Variables.s;
        
        return -1.0/(2.0*M_PI)*k*BesselJ2(k*r)*HKk(k)*Smearing::KernelK(k,s);
    }
    
    
    double HtdIntegrand(double k,void *params){
        
        GSLVariables Variables=*(GSLVariables *) params;
        double r = Variables.r; double s=Variables.s;
        
        return 1.0/(2.0*M_PI)*k*(BesselJ1(k*r)*HKtd(k)+BesselJ2(k*r)/(k*r)*HKtk(k))*Smearing::KernelK(k,s);
    }
    
    double HtmIntegrand(double k,void *params){
        
        GSLVariables Variables=*(GSLVariables *) params;
        double r = Variables.r; double s=Variables.s;
        
        return 1.0/(2.0*M_PI)*k*(BesselJ1(k*r)*HKtm(k)+2.0*BesselJ2(k*r)/(k*r)*HKtk(k))*Smearing::KernelK(k,s);
    }
    
    double HtrIntegrand(double k,void *params){
        
        GSLVariables Variables=*(GSLVariables *) params;
        double r = Variables.r; double s=Variables.s;
        
        return -1.0/(2.0*M_PI)*k*BesselJ3(k*r)*HKtk(k)*Smearing::KernelK(k,s);
    }
    
    
    // COMPUTE INTEGRALS //
    double ComputeHv(double r,double s){
        GSLVariables Variables; Variables.r=r; Variables.s=s;
        EVALUATE_GSL_INTEGRAL(HvIntegrand,Variables);
    }
    
    double ComputeHd(double r,double s){
        GSLVariables Variables; Variables.r=r; Variables.s=s;
        EVALUATE_GSL_INTEGRAL(HdIntegrand,Variables);
    }
    
    double ComputeHr(double r,double s){
        GSLVariables Variables; Variables.r=r; Variables.s=s;
        EVALUATE_GSL_INTEGRAL(HrIntegrand,Variables);
    }
    
    double ComputeHtd(double r,double s){
        GSLVariables Variables; Variables.r=r; Variables.s=s;
        EVALUATE_GSL_INTEGRAL(HtdIntegrand,Variables);
    }
    
    double ComputeHtm(double r,double s){
        GSLVariables Variables; Variables.r=r; Variables.s=s;
        EVALUATE_GSL_INTEGRAL(HtmIntegrand,Variables);
    }
    
    double ComputeHtr(double r,double s){
        GSLVariables Variables; Variables.r=r; Variables.s=s;
        EVALUATE_GSL_INTEGRAL(HtrIntegrand,Variables);
    }
    
    
}
