#ifndef EKTPHOTON_H
#define EKTPHOTON_H

#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


const int NStepsPhotons=244;
const int NP=4000;

class EKTPhoton{
	
	public:
    EKTPhoton()
	{
        // Import data and constants from EKT simulations for lambda=10
        char input1[256];
        char buffer[1024];
        snprintf(input1,256,"EKT_PHOTON/ConstantsLambda10.txt");
        FILE *tableConst = fopen(input1,"r");
        fgets(buffer, 1024,tableConst);
        fscanf(tableConst,"%lf %lf %lf %lf", &lambda10, &tau13T_Infty10, &eta_o_s_10, &intZ4Ceq4_10);
        fclose(tableConst);

        // Include charges in CIdealEKT
        CIdealEKT=intZ4Ceq4_10*sum_qf2;
        
        std::cout<< "#Set constants for import of EKT photon rate:"<< std::endl;
        std::cout<< "#Lambda = 10 => (tau^(1/3)T)_infinity = "<< tau13T_Infty10 << ", eta/s = " << eta_o_s_10 << ", C_(gamma)^(ideal) = " << CIdealEKT << std::endl;
        
        char input2[256];
        snprintf(input2,256,"EKT_PHOTON/EvolutionParametersLambda10.txt");
        FILE *tableEvoPars = fopen(input2,"r");
        fgets(buffer, 1024,tableEvoPars);
        int sc10;
        double tau10;
        double t10;
        double wt10;
        
        for (int i=0; i<NStepsPhotons; i++) {
            fscanf(tableEvoPars,"%d %lf %lf %lf", &sc10, &tau10, &t10, &wt10);
            
            Steps10[i]=sc10;
            Tau10[i]=tau10;
            T10[i]=t10;
            wTilde10[i]=wt10;
        }
        fclose(tableEvoPars);

        for (int i=0; i<NStepsPhotons; i++) {
            import_photons(i);
        }
        
        std::cout<< "Photons imported and universal scaling function set up!" << std::endl;
        std::cout<< "Boundaries of rescaled momentum:" << std::endl;
        std::cerr<< "lambda = 10 => psc- = "<< psc10min << ", psc+ = "<< psc10max << std::endl; 
		
	}
	
	virtual ~EKTPhoton ()
	{
		for(int ij =0; ij< NStepsPhotons; ij++ )
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
        for (int i =0; i< NStepsPhotons; i++) {
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

    // Compute pre-equilibrium photon rates based on universal scaling function in arxiv:2308.09747 for a given eta/s and CIdeal
    // WARNING: tauhydro NEEDS TO BE GIVEN IN GeV HERE
    void compute_photons(double tauhydro, double Thydro, double eta_over_s_input, double CIdeal_input, double *pT, double * Rate, int np){
        
        double WTM, WTP;
        double scaling_x, scaling_y;
        double rate_m, rate_p;

        // Finding Wtilde 
        double wtilde = find_wtilde(tauhydro,Thydro, eta_over_s_input);
        int iw=find_wtilde_index(wtilde,WTM, WTP);

        // Establishing rescaling factors
        scaling_x=pow(eta_over_s_input,1/2.)*pow( find_t13Thydro( tauhydro, Thydro, eta_over_s_input) ,-3/2.); 
        scaling_y=pow(eta_over_s_input,2.)*pow(CIdeal_input,1.);

        // Computing max and min rescaled momenta 
        double pscm_t, pscp_t;
        pscm_t = psc10min;
        pscp_t = psc10max;
        
        // Compute pre-equilibrium spectra (Linear interpolation in the w direction)
        for (int j = 0; j < np; j++)
            {
                if(pT[j]*scaling_x>pscm_t && pT[j]*scaling_x<pscp_t){
                rate_m=scaling_y*gsl_spline_eval(RateInt10[iw], pT[j]*scaling_x, acc10[iw]);
                rate_p=scaling_y*gsl_spline_eval(RateInt10[iw+1], pT[j]*scaling_x, acc10[iw+1]);
                Rate[j]=(wtilde-WTM)*(rate_p-rate_m)/(WTP-WTM) + rate_m;
            }
            else{
                Rate[j]=0; 
            }   
            
        }
        
        
        
        
    }


    // Create output
    int PrintPhotons(const char *out,double tauhydro, double Thydro, double eta_over_s_input, double CIdeal_input, int np, double pmax, double pmin)
	{   

		char fname[1024];
        double wtilde = find_wtilde(tauhydro,Thydro, eta_over_s_input);
        int NF= nfamy;
		snprintf(fname,1024,"%s/PhotonRate_NF_%d_tauhydro_%0.2ffm_Thydro_%0.2fGeV_etas_%0.2f.dat",out,NF,tauhydro*GeVm1tofm, Thydro,eta_over_s_input);
		printf("Writing wf to file %s\n",fname);
		FILE *wfout = fopen(fname,"w");
        fprintf(wfout,"# wtilde %10.5e\n",wtilde);

		double * pT =new double[np];
        double * Rate10 =new double[np];

        double dp = (pmax-pmin)/(np-1);
        for(int k = 0; k < np; k++ ){pT[k] = k *dp +pmin;}
        compute_photons(tauhydro,  Thydro,  eta_over_s_input, CIdeal_input,pT, Rate10, np);
        
		for(int k = 0; k < np; k++ )
		{   
			fprintf(wfout,"%10.3e\t%10.5e\n",pT[k],Rate10[k]);
		}
		
	
		fclose(wfout);
		return 0;
    }
	
	

    private:
    

    
	gsl_interp_accel *acc10[NStepsPhotons];
    gsl_spline *RateInt10[NStepsPhotons];
    
    double lambda10;
    double tau13T_Infty10;
    double eta_o_s_10;
    double intZ4Ceq4_10;
    double CIdealEKT;

    double psc10max,psc10min;
    
    double Steps10[NStepsPhotons];
    double Tau10[NStepsPhotons];
    double T10[NStepsPhotons];
    double wTilde10[NStepsPhotons];

    
    // Import photon spectrum from EKT simulations for lambda=10
    void import_photons(int i){
        
        char inputpho[1024];
        char bufferpho[256];
        int step;
        double scaling_x, scaling_y;
        
        step=Steps10[i];
        scaling_x=pow(tau13T_Infty10,-3./2.)*pow(eta_o_s_10,1./2.);
        scaling_y=pow(eta_o_s_10,-2.)*pow(CIdealEKT,-1.);
        
        snprintf(inputpho,1024,"EKT_PHOTON/LAMBDA_10/pT_SPECTRUM_LEADING_ORDER_RATES/pTSpectrum_LeadingOrderRatePhotons%d.txt", step);
        FILE *tablePhoton = fopen(inputpho,"r");
        fgets(bufferpho,256,tablePhoton);
        
        double pt_t, rate_t;
        double pT[NP];
        double rate[NP];
        
        for (int j = 0; j < NP; j++)
        {
            fscanf(tablePhoton,"%lf %lf", &pt_t, &rate_t);
            pT[j]=pt_t*scaling_x;
            rate[j]=GammaPrefactor*rate_t*scaling_y;
        }

        psc10min= pT[0];
        psc10max= pT[NP-1];
        
        fclose(tablePhoton);

        acc10[i] = gsl_interp_accel_alloc();
        RateInt10[i] = gsl_spline_alloc(gsl_interp_cspline, NP);
        gsl_spline_init (RateInt10[i], &pT[0], &rate[0], NP);
        

    }

    
	
};

#endif
