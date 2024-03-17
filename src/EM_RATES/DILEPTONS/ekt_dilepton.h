#ifndef EKTDILEPTON_H
#define EKTDILEPTON_H

#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

const int NStepsDileptons=244;
const int NM=1000;

class EKTDilepton{
	
	public:
    EKTDilepton()
	{
        // Import data and constants from EKT simulations for lambda=10
        char input1[256];
        char buffer[1024];
        snprintf(input1,256,"EKT_DILEPTON/ConstantsLambda10.txt");
        FILE *tableConst = fopen(input1,"r");
        fgets(buffer, 1024,tableConst);
        fscanf(tableConst,"%lf %lf %lf", &lambda10, &tau13T_Infty10, &eta_o_s_10);
        fclose(tableConst);
        
        std::cout<< "#Set constants for import of EKT dilepton rate:"<< std::endl;
        std::cout<< "#Lambda = 10 => (tau^(1/3)T)_infinity = "<< tau13T_Infty10 << ", eta/s = " << eta_o_s_10 << std::endl;
        
        char input2[256];
        snprintf(input2,256,"EKT_DILEPTON/EvolutionParametersLambda10.txt");
        FILE *tableEvoPars = fopen(input2,"r");
        fgets(buffer, 1024,tableEvoPars);
        int sc10;
        double tau10;
        double t10;
        double wt10;
        
        for (int i=0; i<NStepsDileptons; i++) {
            fscanf(tableEvoPars,"%d %lf %lf %lf", &sc10, &tau10, &t10, &wt10);
            
            Steps10[i]=sc10;
            Tau10[i]=tau10;
            T10[i]=t10;
            wTilde10[i]=wt10;
        }
        fclose(tableEvoPars);

        for (int i=0; i<NStepsDileptons; i++) {
            import_dileptons(i);
        }

        std::cout<< "Dilepton imported and universal scaling function set up!" << std::endl;
        std::cout<< "Boundaries of rescaled mass:" << std::endl;
        std::cerr<< "lambda = 10 => Msc- = "<< Msc10min << ", Msc+ = "<< Msc10max << std::endl; 
		
	}
	
	virtual ~EKTDilepton ()
	{
		for(int ij =0; ij< NStepsDileptons; ij++ )
		{
			gsl_spline_free (RateInt10[ij]);
            gsl_interp_accel_free (acc10[ij]);
		}
	}   


    double get_alpha_s(int lambda){return lambda/(4.*M_PI*NC);}
    
    double find_t13Thydro(double tauhydro, double Thydro, double eta_over_s_input){
        return pow(tauhydro,1/3.)*(Thydro + 2*eta_over_s_input/(3.*tauhydro));
    }

    double find_wtilde(double tauhydro, double Thydro, double eta_over_s_input){
         return tauhydro*Thydro/(4*M_PI*eta_over_s_input);
    }
    
    int find_wtilde_index(double wtil, double &wTm, double &wTp){
        int index_t= 0;
        for (int i =0; i< NStepsDileptons; i++) {
            double wt1,wt2;
            index_t=i;
            wt1=wTilde10[i];wt2=wTilde10[i+1];
            if(wtil>=wt1 && wtil <wt2){
                wTm=wt1;
                wTp=wt2;
                break;
            }
        }
        return index_t;
    }

    // Compute pre-equilibrium dilepton rates based on universal scaling function in arXiv:2403.04846 for a given eta/s
    // WARNING: tauhydro NEEDS TO BE GIVEN IN GeV HERE
    void compute_dileptons(double tauhydro, double Thydro, double eta_over_s_input, double *M, double * Rate, int nM){
        
        double WTM, WTP;
        double scaling_x, scaling_y;
        double rate_m, rate_p;

        // Finding Wtilde 
        double wtilde = find_wtilde(tauhydro,Thydro, eta_over_s_input);
        int iw=find_wtilde_index(wtilde,WTM, WTP);

        // Establishing rescaling factors
        scaling_x=pow(eta_over_s_input,1/2.)*pow( find_t13Thydro( tauhydro, Thydro, eta_over_s_input) ,-3/2.); 
        scaling_y=pow(eta_over_s_input,2.);

        // Computing max and min rescaled mass 
        double Mscm_t, Mscp_t;
        Mscm_t = Msc10min;
        Mscp_t = Msc10max;
        
        // Compute pre-equilibrium spectra (Linear interpolation in the w direction)
        for (int j = 0; j < nM; j++)
            {
                if(M[j]*scaling_x>Mscm_t && M[j]*scaling_x<Mscp_t){
                rate_m=scaling_y*gsl_spline_eval(RateInt10[iw], M[j]*scaling_x, acc10[iw]);
                rate_p=scaling_y*gsl_spline_eval(RateInt10[iw+1], M[j]*scaling_x, acc10[iw+1]);
                Rate[j]=(wtilde-WTM)*(rate_p-rate_m)/(WTP-WTM) + rate_m;
            }
            else{
                Rate[j]=0; 
            }   
            
        }
        
        
        
        
    }


    // Create output
    int PrintDileptons(const char *out,double tauhydro, double Thydro, double eta_over_s_input, int nM, double Mmax, double Mmin)
	{   

		char fname[1024];
        double wtilde = find_wtilde(tauhydro,Thydro,eta_over_s_input);
        int NF= nfamy;
		snprintf(fname,1024,"%s/DileptonRate_NF_%d_tauhydro_%0.2ffm_Thydro_%0.2fGeV_etas_%0.2f.dat",out,NF,tauhydro*GeVm1tofm,Thydro,eta_over_s_input);
		printf("Writing wf to file %s\n",fname);
		FILE *wfout = fopen(fname,"w");
        fprintf(wfout,"# wtilde %10.5e\n",wtilde);

		double * M =new double[nM];
        double * Rate10 =new double[nM];

        double dM = (Mmax-Mmin)/(nM-1);
        for(int k = 0; k < nM; k++ ){M[k] = k *dM +Mmin;}
        compute_dileptons(tauhydro,  Thydro,  eta_over_s_input,M, Rate10, nM);
        
		for(int k = 0; k < nM; k++ )
		{   
			fprintf(wfout,"%10.3e\t%10.5e\n",M[k],Rate10[k]);
		}
		
	
		fclose(wfout);
		return 0;
    }
	
	

    private:
    

    
	gsl_interp_accel *acc10[NStepsDileptons];
    gsl_spline *RateInt10[NStepsDileptons];
    
    double lambda10;
    double tau13T_Infty10;
    double eta_o_s_10;

    double Msc10max,Msc10min;
    
    double Steps10[NStepsDileptons];
    double Tau10[NStepsDileptons];
    double T10[NStepsDileptons];
    double wTilde10[NStepsDileptons];

    
    // Import dilepton spectrum from EKT simulations for lambda=10
    void import_dileptons(int i){
        
        char inputllo[1024];
        char bufferllo[256];
        int step;
        double scaling_x, scaling_y;
        
        step=Steps10[i];
        scaling_x=pow(tau13T_Infty10,-3./2.)*pow(eta_o_s_10,1./2.);
        scaling_y=pow(eta_o_s_10,-2.);
        
        snprintf(inputllo,1024,"EKT_DILEPTON/LAMBDA_10/YIELD_DILEPTONS/Yield_dNMdM_LO_%d.txt", step);
        FILE *tableDilepton = fopen(inputllo,"r");
        fgets(bufferllo,256,tableDilepton);
        
        double M_t, rate_t;
        double M[NM];
        double rate[NM];
        
        for (int j = 0; j < NM; j++)
        {
            fscanf(tableDilepton,"%lf %lf", &M_t, &rate_t);
            M[j]=M_t*scaling_x;
            rate[j]=llPrefactor*rate_t*scaling_y;
        }

        Msc10min= M[0];
        Msc10max= M[NM-1];
        
        fclose(tableDilepton);

        acc10[i] = gsl_interp_accel_alloc();
        RateInt10[i] = gsl_spline_alloc(gsl_interp_cspline, NM);
        gsl_spline_init (RateInt10[i], &M[0], &rate[0], NM);
        

    }

    
	
};

#endif

