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


/////////////////////////////////////////
// COORDINATE SPACE RESPONSE FUNCTIONS //
/////////////////////////////////////////

namespace CoordinateSpace{
    
    //////////////////////////////////////////////
    // EVALUTE RESPONSE KERNEL BY INTERPOLATION //
    //////////////////////////////////////////////
    
    // EVALUATE BETWENN rMin,rMax AND tMin,tMax //
    double rMin,rMax,tMin,tMax;
    
    // GSL INTERPOLATION OBJECTS //
    gsl_interp_accel **HvTAcc,**HdTAcc,**HrTAcc,**HtdTAcc,**HtmTAcc,**HtrTAcc;
    gsl_interp_accel **HvXAcc,**HdXAcc,**HrXAcc,**HtdXAcc,**HtmXAcc,**HtrXAcc;

    gsl_spline2d *HvInt,*HdInt,*HrInt,*HtdInt,*HtmInt,*HtrInt;
    
    // H00,k=Hv r^{k}/r //
    double HvScalingCurve(double r,double t){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(HvInt,r,t,HvXAcc[tID],HvTAcc[tID],rMin,rMax,tMin,tMax);
    }
    double Hv(double dX,double dT,double ScalingVarOut){
        double tVal=pow(ScalingVarOut,3.0/2.0);
        return HvScalingCurve(dX/dT,tVal)/(dT*dT);
    }
    

    // H0i,j=Hd \delta^{ij} + Hr r^{i}r^{j}/r^2 //
    double HdScalingCurve(double r,double t){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(HdInt,r,t,HdXAcc[tID],HdTAcc[tID],rMin,rMax,tMin,tMax);
    }
    
    double Hd(double dX,double dT,double ScalingVarOut){
        double tVal=pow(ScalingVarOut,3.0/2.0);
        return HdScalingCurve(dX/dT,tVal)/(dT*dT);
    }
    
    
    double HrScalingCurve(double r,double t){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(HrInt,r,t,HrXAcc[tID],HrTAcc[tID],rMin,rMax,tMin,tMax);
    }
    
    double Hr(double dX,double dT,double ScalingVarOut){
        double tVal=pow(ScalingVarOut,3.0/2.0);
        return HrScalingCurve(dX/dT,tVal)/(dT*dT);
    }
    
    // Hij,k=Htd\delta^{ij} r^{k}/r + Htm 1/2(r^{i}/r \delta^{jk} + r^{j}/r \delta^{jk}) + Htr r^{i}r^{j}kr^{k}/r^3 //
    
    double HtdScalingCurve(double r,double t){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(HtdInt,r,t,HtdXAcc[tID],HtdTAcc[tID],rMin,rMax,tMin,tMax);
    }
    
    double Htd(double dX,double dT,double ScalingVarOut){
        double tVal=pow(ScalingVarOut,3.0/2.0);
        return HtdScalingCurve(dX/dT,tVal)/(dT*dT);
    }
    
    double HtmScalingCurve(double r,double t){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(HtmInt,r,t,HtmXAcc[tID],HtmTAcc[tID],rMin,rMax,tMin,tMax);
    }
    
    double Htm(double dX,double dT,double ScalingVarOut){
        double tVal=pow(ScalingVarOut,3.0/2.0);
        return HtmScalingCurve(dX/dT,tVal)/(dT*dT);
    }
    
    
    double HtrScalingCurve(double r,double t){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(HtrInt,r,t,HtrXAcc[tID],HtrTAcc[tID],rMin,rMax,tMin,tMax);
    }
    
    double Htr(double dX,double dT,double ScalingVarOut){
        double tVal=pow(ScalingVarOut,3.0/2.0);
        return HtrScalingCurve(dX/dT,tVal)/(dT*dT);
    }
    
    /////////////////////////////////////////////////////////
    // COMPUTE SMEARED COORDINATE SPACE RESPONSE FUNCTIONS //
    /////////////////////////////////////////////////////////
    
    double *rValues,*tValues;
    
    double *HvValues;
    double *HdValues;
    double *HrValues;
    double *HtdValues;
    double *HtmValues;
    double *HtrValues;
    
    
    void SetupGrid(int NumberOfPoints,int NumberOfTimes,double rMinVal,double rMaxVal){
        
        // ALLOCATE //
        rValues=new double[NumberOfPoints];
        tValues=new double[NumberOfTimes];
        
        HvValues=new double[NumberOfTimes*NumberOfPoints];
        HdValues=new double[NumberOfTimes*NumberOfPoints];
        HrValues=new double[NumberOfTimes*NumberOfPoints];
        HtdValues=new double[NumberOfTimes*NumberOfPoints];
        HtmValues=new double[NumberOfTimes*NumberOfPoints];
        HtrValues=new double[NumberOfTimes*NumberOfPoints];
        
        // SET r GRID //
        for(int rIndex=0;rIndex<NumberOfPoints;rIndex++){
            double r=rMinVal+(rMaxVal-rMinVal)*(rIndex+0.5)/(NumberOfPoints);  rValues[rIndex]=r;
        }
        
        // SETUP GSL INTERPOLATION //
        int NumberOfOpenMPThreads=omp_get_max_threads();
        
        HvXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        HdXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        HrXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        HtdXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        HtmXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        HtrXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        
        
        HvTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        HdTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        HrTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        HtdTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        HtmTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        HtrTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];

        
        #pragma omp parallel for
        for(int i=0;i<NumberOfOpenMPThreads;i++){
            
            HvXAcc[i] = gsl_interp_accel_alloc ();
            HdXAcc[i] = gsl_interp_accel_alloc ();
            HrXAcc[i] = gsl_interp_accel_alloc ();
            HtdXAcc[i] = gsl_interp_accel_alloc ();
            HtmXAcc[i] = gsl_interp_accel_alloc ();
            HtrXAcc[i] = gsl_interp_accel_alloc ();
            
            HvTAcc[i] = gsl_interp_accel_alloc ();
            HdTAcc[i] = gsl_interp_accel_alloc ();
            HrTAcc[i] = gsl_interp_accel_alloc ();
            HtdTAcc[i] = gsl_interp_accel_alloc ();
            HtmTAcc[i] = gsl_interp_accel_alloc ();
            HtrTAcc[i] = gsl_interp_accel_alloc ();
        }
        
        HvInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints,NumberOfTimes);
        HdInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints,NumberOfTimes);
        HrInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints,NumberOfTimes);
        HtdInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints,NumberOfTimes);
        HtmInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints,NumberOfTimes);
        HtrInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints,NumberOfTimes);
        
    }
    
    void SetValues(int NumberOfPoints,int tIndex,double tValue,double Sigma){
        
        // SET t VALUE //
        tValues[tIndex]=tValue;
        
        // SETUP VALUES OF RESPONSE FUNCTIONS //
        //#pragma omp parallel for
        for(int rIndex=0;rIndex<NumberOfPoints;rIndex++){
            
            double r=rValues[rIndex];
            
            gsl_spline2d_set(HvInt,HvValues,rIndex,tIndex,BesselTransform::ComputeHv(r,Sigma));
            gsl_spline2d_set(HdInt,HdValues,rIndex,tIndex,BesselTransform::ComputeHd(r,Sigma));
            gsl_spline2d_set(HrInt,HrValues,rIndex,tIndex,BesselTransform::ComputeHr(r,Sigma));
            gsl_spline2d_set(HtdInt,HtdValues,rIndex,tIndex,BesselTransform::ComputeHtd(r,Sigma));
            gsl_spline2d_set(HtmInt,HtmValues,rIndex,tIndex,BesselTransform::ComputeHtm(r,Sigma));
            gsl_spline2d_set(HtrInt,HtrValues,rIndex,tIndex,BesselTransform::ComputeHtr(r,Sigma));
            
            
        }
        
    }
    
    /////////////////////////
    // SETUP INTERPOLATORS //
    /////////////////////////
    
    void SetupInterpolators(int NumberOfPoints,int NumberOfTimes){

        
        // SET BOUNDARIES //
        rMin=rValues[0]; rMax=rValues[NumberOfPoints-1];
        tMin=tValues[0]; tMax=tValues[NumberOfTimes-1];
        
        gsl_spline2d_init(HvInt,rValues,tValues,HvValues,NumberOfPoints,NumberOfTimes);
        gsl_spline2d_init(HdInt,rValues,tValues,HdValues,NumberOfPoints,NumberOfTimes);
        gsl_spline2d_init(HrInt,rValues,tValues,HrValues,NumberOfPoints,NumberOfTimes);
        gsl_spline2d_init(HtdInt,rValues,tValues,HtdValues,NumberOfPoints,NumberOfTimes);
        gsl_spline2d_init(HtmInt,rValues,tValues,HtmValues,NumberOfPoints,NumberOfTimes);
        gsl_spline2d_init(HtrInt,rValues,tValues,HtrValues,NumberOfPoints,NumberOfTimes);
        
        // CLEAN-UP //
        delete[] rValues;
        delete[] tValues;
        
        delete[] HvValues;
        delete[] HdValues;
        delete[] HrValues;
        delete[] HtdValues;
        delete[] HtmValues;
        delete[] HtrValues;
        
    }
    
    
    ///////////////////
    // CREATE OUTPUT //
    ///////////////////
    
    void Output(std::string fname,int NdX,int NdT,double t0,double t1,double K,double EtaOverS0){
        
        std::ofstream Outstream;
        Outstream.open(fname.c_str());
        Outstream << "# tau_R(tau=1fm) = " << pow(K/M_HBARC,2.0/3.0)/EtaOverS0 <<std::endl;
        Outstream << "# 1:|x-x_0| 2:|tau-tau0| 3:Gvs 4:Gvvd 5:Gvvr 6:Gvtd 7:Gvtn 8:Gvtr" << std::endl;
        
        for(int tIndex=0;tIndex<NdT;tIndex++){
            
            double dT=t0+(tIndex+0.5)*(t1-t0)/(NdT);
            
            for(int xIndex=0;xIndex<NdX;xIndex++){
                
                double dX=2.0*(xIndex+0.5)*dT/(NdX);

                double ScalingVarOut=pow(K*dT/M_HBARC,2.0/3.0)/EtaOverS0 ;
                
                Outstream << dX << " " << dT  << " "
                   << Hv(dX,dT,ScalingVarOut) << " "
                   << Hd(dX,dT,ScalingVarOut) << " "
                   << Hr(dX,dT,ScalingVarOut) << " "
                   << Htd(dX,dT,ScalingVarOut) << " "
                   << Htm(dX,dT,ScalingVarOut) << " "
                   << Htr(dX,dT,ScalingVarOut) << std::endl;
                
            }
            
            Outstream << std::endl;
            Outstream << std::endl;
            
        }
        
        Outstream.close();
    }
    
    
}
