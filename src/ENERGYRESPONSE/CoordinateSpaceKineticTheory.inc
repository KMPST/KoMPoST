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
#ifndef M_HBARC
#define M_HBARC  0.197326979
#endif


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
    gsl_interp_accel **GsTAcc,**GvTAcc,**GdTAcc,**GrTAcc;
    gsl_interp_accel **GsXAcc,**GvXAcc,**GdXAcc,**GrXAcc;

    gsl_spline2d *GsInt,*GvInt,*GdInt,*GrInt;
    
    // G00=Gs //
    double GsScalingCurve(double r,double t){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GsInt,r,t,GsXAcc[tID],GsTAcc[tID],rMin,rMax,tMin,tMax);
    }
    
    double Gs(double dX,double dT,double ScalingVarOut){
        double tVal=pow(ScalingVarOut,3.0/2.0);
        return GsScalingCurve(dX/dT,tVal)/(dT*dT);
    }
    
    // G0i=r^{i}/r Gv //
    double GvScalingCurve(double r,double t){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GvInt,r,t,GvXAcc[tID],GvTAcc[tID],rMin,rMax,tMin,tMax);
    }
    double Gv(double dX,double dT,double ScalingVarOut){
        double tVal=pow(ScalingVarOut,3.0/2.0);
        return GvScalingCurve(dX/dT,tVal)/(dT*dT);
    }
    
    // Gij=Gd \delta^{ij} + Gr r^{i}r^{j}/r^2 //
    double GdScalingCurve(double r,double t){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GdInt,r,t,GdXAcc[tID],GdTAcc[tID],rMin,rMax,tMin,tMax);
    }
    
    double Gd(double dX,double dT,double ScalingVarOut){
        double tVal=pow(ScalingVarOut,3.0/2.0);
        return GdScalingCurve(dX/dT,tVal)/(dT*dT);
    }
    
    
    double GrScalingCurve(double r,double t){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR_2D(GrInt,r,t,GrXAcc[tID],GrTAcc[tID],rMin,rMax,tMin,tMax);
    }
    
    double Gr(double dX,double dT,double ScalingVarOut){
        double tVal=pow(ScalingVarOut,3.0/2.0);
        return GrScalingCurve(dX/dT,tVal)/(dT*dT);
    }
    
    /////////////////////////////////////////////////////////
    // COMPUTE SMEARED COORDINATE SPACE RESPONSE FUNCTIONS //
    /////////////////////////////////////////////////////////
    
    double *rValues,*tValues;
    
    double *GsValues;
    double *GvValues;
    double *GdValues;
    double *GrValues;
    
    void SetupGrid(int NumberOfPoints,int NumberOfTimes,double rMinVal,double rMaxVal){
        
        // ALLOCATE //
        rValues=new double[NumberOfPoints];
        tValues=new double[NumberOfTimes];
        
        GsValues=new double[NumberOfTimes*NumberOfPoints];
        GvValues=new double[NumberOfTimes*NumberOfPoints];
        GdValues=new double[NumberOfTimes*NumberOfPoints];
        GrValues=new double[NumberOfTimes*NumberOfPoints];
        
        // SET r GRID //
        for(int rIndex=0;rIndex<NumberOfPoints;rIndex++){
            double r=rMinVal+(rMaxVal-rMinVal)*(rIndex+0.5)/(NumberOfPoints);  rValues[rIndex]=r;
        }
        
        // SETUP GSL INTERPOLATION //
        int NumberOfOpenMPThreads=omp_get_max_threads();
        
        GsXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GvXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GdXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GrXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        
        GsTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GvTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GdTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GrTAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        
        #pragma omp parallel for
        for(int i=0;i<NumberOfOpenMPThreads;i++){
            
            GsXAcc[i] = gsl_interp_accel_alloc ();
            GvXAcc[i] = gsl_interp_accel_alloc ();
            GdXAcc[i] = gsl_interp_accel_alloc ();
            GrXAcc[i] = gsl_interp_accel_alloc ();
            
            GsTAcc[i] = gsl_interp_accel_alloc ();
            GvTAcc[i] = gsl_interp_accel_alloc ();
            GdTAcc[i] = gsl_interp_accel_alloc ();
            GrTAcc[i] = gsl_interp_accel_alloc ();
        }
        
        GsInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints,NumberOfTimes);
        GvInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints,NumberOfTimes);
        GdInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints,NumberOfTimes);
        GrInt=gsl_spline2d_alloc(gsl_interp2d_bilinear,NumberOfPoints,NumberOfTimes);
        
    }
    
    void SetValues(int NumberOfPoints,int tIndex,double tValue,double Sigma){
        
        // SET t VALUE //
        tValues[tIndex]=tValue;
        
        // SETUP VALUES OF RESPONSE FUNCTIONS //
        #pragma omp parallel for
        for(int rIndex=0;rIndex<NumberOfPoints;rIndex++){
            
            double r=rValues[rIndex];
            
            gsl_spline2d_set(GsInt,GsValues,rIndex,tIndex,BesselTransform::ComputeGs(r,Sigma));
            gsl_spline2d_set(GvInt,GvValues,rIndex,tIndex,BesselTransform::ComputeGv(r,Sigma));
            gsl_spline2d_set(GdInt,GdValues,rIndex,tIndex,BesselTransform::ComputeGd(r,Sigma));
            gsl_spline2d_set(GrInt,GrValues,rIndex,tIndex,BesselTransform::ComputeGr(r,Sigma));

            
        }
        
    }
    
    /////////////////////////
    // SETUP INTERPOLATORS //
    /////////////////////////
    
    void SetupInterpolators(int NumberOfPoints,int NumberOfTimes){

        
        // SET BOUNDARIES //
        rMin=rValues[0]; rMax=rValues[NumberOfPoints-1];
        tMin=tValues[0]; tMax=tValues[NumberOfTimes-1];
        
        gsl_spline2d_init(GsInt,rValues,tValues,GsValues,NumberOfPoints,NumberOfTimes);
        gsl_spline2d_init(GvInt,rValues,tValues,GvValues,NumberOfPoints,NumberOfTimes);
        gsl_spline2d_init(GdInt,rValues,tValues,GdValues,NumberOfPoints,NumberOfTimes);
        gsl_spline2d_init(GrInt,rValues,tValues,GrValues,NumberOfPoints,NumberOfTimes);
        
        // CLEAN-UP //
        delete[] rValues;
        delete[] tValues;
        
        delete[] GsValues;
        delete[] GvValues;
        delete[] GdValues;
        delete[] GrValues;
        
    }
    
    
    ///////////////////
    // CREATE OUTPUT //
    ///////////////////
    
    void Output(std::string fname,int NdX,int NdT,double t0,double t1,double K,double EtaOverS0){
        
        std::ofstream Outstream;
        Outstream.open(fname.c_str());
        Outstream << "# tau_R(tau=1fm) = " << pow(K/M_HBARC,2.0/3.0)/EtaOverS0 <<std::endl;
        Outstream << "# 1:|x-x_0| 2:|tau-tau0| 3:Gss 4:Gsv 5:Gsd 6:Gsr" << std::endl;
        
        for(int tIndex=0;tIndex<NdT;tIndex++){
            
            double dT=t0+(tIndex+0.5)*(t1-t0)/(NdT);
            
            for(int xIndex=0;xIndex<NdX;xIndex++){
                
                double dX=2.0*(xIndex+0.5)*dT/(NdX);

                double ScalingVarOut =pow(K*dT/M_HBARC,2.0/3.0)/EtaOverS0 ;
                
                Outstream << dX << " " << dT  << " "
                   << Gs(dX,dT,ScalingVarOut) << " "
                   << Gv(dX,dT,ScalingVarOut) << " "
                   << Gd(dX,dT,ScalingVarOut) << " "
                   << Gr(dX,dT,ScalingVarOut) << std::endl;
                
            }
            
            Outstream << std::endl;
            Outstream << std::endl;
            
        }
        
        Outstream.close();
    }
    
    
}
