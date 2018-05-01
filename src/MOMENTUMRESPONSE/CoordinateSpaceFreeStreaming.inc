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
    
    // EVALUATE BETWENN rMin AND rMax //
    double rMin,rMax;
    
    // GSL INTERPOLATION OBJECTS //
    gsl_interp_accel **HvAcc,**HdAcc,**HrAcc,**HtdAcc,**HtmAcc,**HtrAcc;
    gsl_spline *HvInt,*HdInt,*HrInt,*HtdInt,*HtmInt,*HtrInt;
    
    // H00,k=Hv r^{k}/r //
    double HvScalingCurve(double r){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR(HvInt,r,HvAcc[tID],rMin,rMax);
    }
    
    double Hv(double dX,double dT){
        return HvScalingCurve(dX/dT)/(dT*dT);
    }
    
    // H0i,j=Hd \delta^{ij} + Hr r^{i}r^{j}/r^2 //
    double HdScalingCurve(double r){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR(HdInt,r,HdAcc[tID],rMin,rMax);
    }
    
    double Hd(double dX,double dT){
        return HdScalingCurve(dX/dT)/(dT*dT);
    }
    
    double HrScalingCurve(double r){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR(HrInt,r,HrAcc[tID],rMin,rMax);
    }
    
    double Hr(double dX,double dT){
        return HrScalingCurve(dX/dT)/(dT*dT);
    }
    
    // Hij,k=Htd\delta^{ij} r^{k}/r + Htm 1/2(r^{i}/r \delta^{jk} + r^{j}/r \delta^{jk}) + Htr r^{i}r^{j}kr^{k}/r^3
    double HtdScalingCurve(double r){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR(HtdInt,r,HtdAcc[tID],rMin,rMax);
    }
    
    double Htd(double dX,double dT){
        return HtdScalingCurve(dX/dT)/(dT*dT);
    }
    
    double HtmScalingCurve(double r){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR(HtmInt,r,HtmAcc[tID],rMin,rMax);
    }
    
    double Htm(double dX,double dT){
        return HtmScalingCurve(dX/dT)/(dT*dT);
    }
    
    double HtrScalingCurve(double r){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR(HtrInt,r,HtrAcc[tID],rMin,rMax);
    }
    
    double Htr(double dX,double dT){
        return HtrScalingCurve(dX/dT)/(dT*dT);
    }
    
    
    /////////////////////////////////////////////////////////
    // COMPUTE SMEARED COORDINATE SPACE RESPONSE FUNCTIONS //
    /////////////////////////////////////////////////////////
    
    void Setup(int NumberOfPoints,double rMinVal,double rMaxVal,double Sigma){
     
        
        // SETUP DISCRETE r GRID //
        double *rValues=new double[NumberOfPoints];
        double *HvValues=new double[NumberOfPoints];
        double *HdValues=new double[NumberOfPoints];
        double *HrValues=new double[NumberOfPoints];
        double *HtdValues=new double[NumberOfPoints];
        double *HtmValues=new double[NumberOfPoints];
        double *HtrValues=new double[NumberOfPoints];
        
        
        #pragma omp parallel for
        for(int i=0;i<NumberOfPoints;i++){
            
            double r=rMinVal+(rMaxVal-rMinVal)*(i+0.5)/(NumberOfPoints);  rValues[i]=r;
            
            HvValues[i]=BesselTransform::ComputeHv(r,Sigma);
            HdValues[i]=BesselTransform::ComputeHd(r,Sigma);
            HrValues[i]=BesselTransform::ComputeHr(r,Sigma);
            HtdValues[i]=BesselTransform::ComputeHtd(r,Sigma);
            HtmValues[i]=BesselTransform::ComputeHtm(r,Sigma);
            HtrValues[i]=BesselTransform::ComputeHtr(r,Sigma);
            
        }
        
        /////////////////////////
        // SETUP INTERPOLATORS //
        /////////////////////////
        
        // SET BOUNDARY //
        rMin=rValues[0]; rMax=rValues[NumberOfPoints-1];
        
        // SETUP GSL INTERPOLATION //
        int NumberOfOpenMPThreads=omp_get_max_threads();
        
        HvAcc =new gsl_interp_accel*[NumberOfOpenMPThreads];
        HdAcc =new gsl_interp_accel*[NumberOfOpenMPThreads];
        HrAcc =new gsl_interp_accel*[NumberOfOpenMPThreads];
        HtdAcc =new gsl_interp_accel*[NumberOfOpenMPThreads];
        HtmAcc =new gsl_interp_accel*[NumberOfOpenMPThreads];
        HtrAcc =new gsl_interp_accel*[NumberOfOpenMPThreads];
        
        #pragma omp parallel for
        for(int i=0;i<NumberOfOpenMPThreads;i++){
            
            HvAcc[i] = gsl_interp_accel_alloc ();
            HdAcc[i] = gsl_interp_accel_alloc ();
            HrAcc[i] = gsl_interp_accel_alloc ();
            HtdAcc[i] = gsl_interp_accel_alloc ();
            HtmAcc[i] = gsl_interp_accel_alloc ();
            HtrAcc[i] = gsl_interp_accel_alloc ();
            
        }
        
        
        HvInt=gsl_spline_alloc(gsl_interp_akima,NumberOfPoints);
        HdInt=gsl_spline_alloc(gsl_interp_akima,NumberOfPoints);
        HrInt=gsl_spline_alloc(gsl_interp_akima,NumberOfPoints);
        HtdInt=gsl_spline_alloc(gsl_interp_akima,NumberOfPoints);
        HtmInt=gsl_spline_alloc(gsl_interp_akima,NumberOfPoints);
        HtrInt=gsl_spline_alloc(gsl_interp_akima,NumberOfPoints);
        
        gsl_spline_init(HvInt,rValues,HvValues,NumberOfPoints);
        gsl_spline_init(HdInt,rValues,HdValues,NumberOfPoints);
        gsl_spline_init(HrInt,rValues,HrValues,NumberOfPoints);
        gsl_spline_init(HtdInt,rValues,HtdValues,NumberOfPoints);
        gsl_spline_init(HtmInt,rValues,HtmValues,NumberOfPoints);
        gsl_spline_init(HtrInt,rValues,HtrValues,NumberOfPoints);
        
        // CLEAN-UP //
        delete[] rValues;
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
    
    void Output(std::string fname,int NdX,int NdT,double t0,double t1){
        
        std::ofstream Outstream;
        Outstream.open(fname.c_str());
        Outstream << "# 1:|x-x_0| 2:|tau-tau0| 3:Gvs 4:Gvvd 5:Gvvr 6:Gvtd 7:Gvtn 8:Gvtr" << std::endl;
        
        for(int tIndex=0;tIndex<NdT;tIndex++){
            
            double dT=t0+(tIndex+0.5)*(t1-t0)/(NdT);
            
            for(int xIndex=0;xIndex<NdX;xIndex++){
                
                double dX=2.0*(xIndex+0.5)*dT/(NdX);
                
                Outstream << dX << " " << dT  << " " << Hv(dX,dT) << " " << Hd(dX,dT) << " " << Hr(dX,dT) << " " << Htd(dX,dT) << " " << Htm(dX,dT) << " " << Htr(dX,dT)  << std::endl;
                
            }
            
            Outstream << std::endl;
            Outstream << std::endl;
            
        }
        
        Outstream.close();
    }
    
}


