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
    gsl_interp_accel **GsAcc,**GvAcc,**GdAcc,**GrAcc;
    gsl_spline *GsInt,*GvInt,*GdInt,*GrInt;
    
    // G00=Gs //
    double GsScalingCurve(double r){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR(GsInt,r,GsAcc[tID],rMin,rMax);
    }
    
    double Gs(double dX,double dT){
        return GsScalingCurve(dX/dT)/(dT*dT);
    }
    
    // G0i=r^{i}/r Gv //
    double GvScalingCurve(double r){

        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR(GvInt,r,GvAcc[tID],rMin,rMax);
    }
    
    double Gv(double dX,double dT){
        return GvScalingCurve(dX/dT)/(dT*dT);
    }
    
    // Gij=Gd \delta^{ij} + Gr r^{i}r^{j}/r^2 //
    double GdScalingCurve(double r){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR(GdInt,r,GdAcc[tID],rMin,rMax);
    }
    
    double Gd(double dX,double dT){
        return GdScalingCurve(dX/dT)/(dT*dT);
    }
    
    double GrScalingCurve(double r){
        
        int tID=omp_get_thread_num();
        EVALUATE_GSL_INTERPOLATOR(GrInt,r,GrAcc[tID],rMin,rMax);
    }
    
    double Gr(double dX,double dT){
        return GrScalingCurve(dX/dT)/(dT*dT);
    }
    
    /////////////////////////////////////////////////////////
    // COMPUTE SMEARED COORDINATE SPACE RESPONSE FUNCTIONS //
    /////////////////////////////////////////////////////////
    
    void Setup(int NumberOfPoints,double rMinVal,double rMaxVal,double Sigma){
        
        // SETUP DISCRETE r GRID //
        double *rValues=new double[NumberOfPoints];
        double *GsValues=new double[NumberOfPoints];
        double *GvValues=new double[NumberOfPoints];
        double *GdValues=new double[NumberOfPoints];
        double *GrValues=new double[NumberOfPoints];
        
        #pragma omp parallel for
        for(int i=0;i<NumberOfPoints;i++){
            
            double r=rMinVal+(rMaxVal-rMinVal)*(i+0.5)/(NumberOfPoints);  rValues[i]=r;
            
            GsValues[i]=BesselTransform::ComputeGs(r,Sigma);
            GvValues[i]=BesselTransform::ComputeGv(r,Sigma);
            GdValues[i]=BesselTransform::ComputeGd(r,Sigma);
            GrValues[i]=BesselTransform::ComputeGr(r,Sigma);
            
        }
        
        /////////////////////////
        // SETUP INTERPOLATORS //
        /////////////////////////
        
        // SET BOUNDARY //
        rMin=rValues[0]; rMax=rValues[NumberOfPoints-1];
        
        // SETUP GSL INTERPOLATION //
        int NumberOfOpenMPThreads=omp_get_max_threads();
        
        GsAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GvAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GdAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        GrAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        
        #pragma omp parallel for
        for(int i=0;i<NumberOfOpenMPThreads;i++){
            
            GsAcc[i] = gsl_interp_accel_alloc ();
            GvAcc[i] = gsl_interp_accel_alloc ();
            GdAcc[i] = gsl_interp_accel_alloc ();
            GrAcc[i] = gsl_interp_accel_alloc ();
        }
        
        GsInt=gsl_spline_alloc(gsl_interp_akima,NumberOfPoints);
        GvInt=gsl_spline_alloc(gsl_interp_akima,NumberOfPoints);
        GdInt=gsl_spline_alloc(gsl_interp_akima,NumberOfPoints);
        GrInt=gsl_spline_alloc(gsl_interp_akima,NumberOfPoints);
        
        gsl_spline_init(GsInt,rValues,GsValues,NumberOfPoints);
        gsl_spline_init(GvInt,rValues,GvValues,NumberOfPoints);
        gsl_spline_init(GdInt,rValues,GdValues,NumberOfPoints);
        gsl_spline_init(GrInt,rValues,GrValues,NumberOfPoints);
        
        // CLEAN-UP //
        delete[] rValues;
        delete[] GsValues;
        delete[] GvValues;
        delete[] GdValues;
        delete[] GrValues;
        
    }
    
    
    
    ///////////////////
    // CREATE OUTPUT //
    ///////////////////
    
    void Output(std::string fname,int NdX,int NdT,double t0,double t1){
        
        std::ofstream Outstream;
        Outstream.open(fname.c_str());
        
        Outstream << "# 1:|x-x_0| 2:|tau-tau0| 3:Gss 4:Gsv 5:Gsd 6:Gsr" << std::endl;
        for(int tIndex=0;tIndex<NdT;tIndex++){
            
            double dT=t0+(tIndex+0.5)*(t1-t0)/(NdT);
            
            for(int xIndex=0;xIndex<NdX;xIndex++){
                
                double dX=2.0*(xIndex+0.5)*dT/(NdX);
                
                Outstream << dX << " " << dT  << " " << Gs(dX,dT) << " " << Gv(dX,dT) << " " << Gd(dX,dT) << " " << Gr(dX,dT) << std::endl;
                
            }
            
            Outstream << std::endl;
            Outstream << std::endl;
            
        }
        
        Outstream.close();
    }
    
    
    
}
