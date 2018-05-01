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
#ifndef ComputeSigma_cpp
#define ComputeSigma_cpp

#include <gsl/gsl_math.h>
#include <cmath>

#include "EventInput.h"
#include "EnergyMomentumTensor.h"

// Flux limiter for derivative: min-mod (kMinmod=1) or no limiter (kNolim=0)
class TLimiter
{
	private:
		int fKmethod;
		double minmod(const double &dpm, const double &dpp) const
		{ double r = dpm*dpp/GSL_MAX(dpm*dpm, GSL_DBL_MIN);
			return dpm*GSL_MAX(0.0, GSL_MIN(1.0,r)); }
		double nolim(const double &dpm, const double &dpp) const
		{ return (dpm + dpp)/2.;  } 

	public:
		enum EMethod: int { kNolim=0, kMinmod=1};
		TLimiter(): fKmethod(kMinmod){};
		TLimiter(int kmethod): fKmethod(kmethod){};
		double D(double p, double pm, double pp) const
		{
			switch (fKmethod){
				case EMethod::kNolim:
					return nolim( p-pm, pp-p );
					break;
				case EMethod::kMinmod:
					return minmod(p-pm, pp-p);
					break;
				default:
					std::cout << "** TLimiter::D ** Bad method selctor in switch!" << std::endl;
					return 0.;
			}
		}

};

// Spatial derivative of variable given by index "m" (index of "vector variable") in direction "direction"
double DIp(const int &direction, class EnergyMomentumTensorMap *TOut_Full, class EnergyMomentumTensorMap *TOut_BG,
		const int &i,
		const int &j,
		const int &k,
		const int &m)
{
	// Forward declaration
	double GetPData(class EnergyMomentumTensorMap *TOut_Full, class EnergyMomentumTensorMap *TOut_BG, int var_id, int xS, int yS, int zS);

	//TLimiter lim(TLimiter::kMinmod);
	TLimiter lim(TLimiter::kNolim);
	switch(direction){
		//x direction
		case 0:
			//return lim.D(TOut->GetPData(m,i,j,k) , TOut->GetPData(m,i-1,j,k),TOut->GetPData(m,i+1,j,k));
			return lim.D(GetPData(TOut_Full,TOut_BG,m,i,j,k), GetPData(TOut_Full,TOut_BG,m,i-1,j,k), GetPData(TOut_Full,TOut_BG,m,i+1,j,k));
			break;
			//y direction
		case 1:
			//return lim.D(TOut->GetPData(m,i,j,k) , TOut->GetPData(m,i,j-1,k),TOut->GetPData(m,i,j+1,k));
			return lim.D(GetPData(TOut_Full,TOut_BG,m,i,j,k), GetPData(TOut_Full,TOut_BG,m,i,j-1,k), GetPData(TOut_Full,TOut_BG,m,i,j+1,k));
			break;
			//z direction
		case 2:
			//return lim.D(TOut->GetPData(m,i,j,k) , TOut->GetPData(m,i,j,k-1),TOut->GetPData(m,i,j,k+1));
			return lim.D(GetPData(TOut_Full,TOut_BG,m,i,j,k), GetPData(TOut_Full,TOut_BG,m,i,j,k-1), GetPData(TOut_Full,TOut_BG,m,i,j,k+1));
			break;
		default:
			std::cerr << "** TBRSSSStep::DIp ** Unrecognized direction \n" << direction << std::endl ;
			abort();
	}
	return 0;
}


//! Implements 4-vector acceleration (see above) using the ideal equations of motion for time derivative
void accelerationIdeal( class EnergyMomentumTensorMap *TOut_Full, class EnergyMomentumTensorMap *TOut_BG,
		const int &i, const int &j, const int &k,
		const double (&dxI)[3],
		const double (&u)[4],
		const double (&DIuJ)[3][3],
		// const TNDArray4D &DaI,
		double (&a)[4],
		double &theta)
{
	// Forward declaration
	double GetPData(class EnergyMomentumTensorMap *TOut_Full, class EnergyMomentumTensorMap *TOut_BG, int var_id, int xS, int yS, int zS);

	// Evaluate DIuI, uIuJDIuJ and uIDIuJ
	double DIuI = 0; 
	double uIuJDIuJ = 0;
	double uIDIuJ[3] = {0., 0., 0.} ;
	for (int J=0; J<3; J++) {
		for (int I=0; I<3; I++) {
			uIDIuJ[J] += u[1+I]*DIuJ[I][J];
		}
		uIuJDIuJ += u[1+J]*uIDIuJ[J] ;
		DIuI += DIuJ[J][J];
	}

	// For extracting the time derivatives using the ideal EOM
	// we need spatial derivatives of energy density and the speed of sound
	double DIe[3];
	for (int I=0; I<3; I++){
		DIe[I] = DIp(I, TOut_Full, TOut_BG, i, j, k, 0 ) /dxI[I];
	}
	//double e = TOut->GetEd(i,j);
	double e = GetPData(TOut_Full,TOut_BG,0,i,j,k);
	double p = e/3.;
	double cs = sqrt(1./3.);
	double cs2 = cs*cs;

	// Construct the acceleration using the ideal EOM
	const double cs2ep = cs2/(e+p);
	const double uIDIe = (u[1]*DIe[0]+u[2]*DIe[1]+u[3]*DIe[2]);
	const double term2 = cs2/((1.-cs2)*u[0]*u[0]+cs2)*
		(cs2ep * uIDIe + uIuJDIuJ - u[0]*u[0]*DIuI);

	// We now can determine the acceleration,
	a[1] = -cs2ep*DIe[0] - u[1]*term2 ; //+ DaI(i,j,k,0);
	a[2] = -cs2ep*DIe[1] - u[2]*term2 ; //+ DaI(i,j,k,1);
	a[3] = -cs2ep*DIe[2] - u[3]*term2 ; //+ DaI(i,j,k,2);
	a[0] = (a[1]*u[1]+a[2]*u[2]+a[3]*u[3])/u[0];

	// // Store the a correction based on this approximate time derivative.
	// // updateAccelerationCorrection will modify this value after the stage
	// // to determine the difference between the actual time derivaive  of
	// // the step and the approximate derivative
	// DaI(i,j,k,0) = -(-cs2ep*DIe[0] - u[1]*term2 - uIDIuJ[0]);
	// DaI(i,j,k,1) = -(-cs2ep*DIe[1] - u[2]*term2 - uIDIuJ[1]);
	// DaI(i,j,k,2) = -(-cs2ep*DIe[2] - u[3]*term2 - uIDIuJ[2]);

	// Extract theta using D_{\tau} u^\tau, using a[0]
	theta = a[0]/u[0]-uIuJDIuJ/(u[0]*u[0])+DIuI;
}
const int NPi44=10;
const int NOmega44=6;

//Get \sigma^{\mu\nu}
void GetSigma(class EnergyMomentumTensorMap *TOut_Full, class EnergyMomentumTensorMap *TOut_BG,
		const int i,
		const int j,
		const int k,
		double sg[NPi44],
		double om[NOmega44]
	     ) {
	// Forward declaration
	double GetPData(class EnergyMomentumTensorMap *TOut_Full, class EnergyMomentumTensorMap *TOut_BG, int var_id, int xS, int yS, int zS);

	double tau = TOut_Full->tau; //fm
	double dx,dy,deta;
	dx = EventInput::afm;
	dy = EventInput::afm;
	deta = 1.0;
	double dxI[3];
	dxI[0] = dx;
	dxI[1] = dy;
	dxI[2] = tau*deta;

	// u is the flow, with index raise
	// u = (u^tau, u^x, u^y, tau u^eta)
	double u[4];
	int kU0=1;
	for (int m=0; m < 4; m++){
		//u[m] = TOut->GetPData(kU0+m,i,j,k) ;
		u[m] = GetPData(TOut_Full,TOut_BG,kU0+m,i,j,k);
	}

	//Find local derivatives for accelleration and sigma/omega
	// pn = flow
	double DIuJ[3][3];
	int kUI = 2;
	for (int J=0; J<3; J++){
		for (int I=0; I<3; I++){
			DIuJ[I][J] = DIp(I, TOut_Full, TOut_BG, i, j, k, kUI + J)/dxI[I] ;
		}
	}
	DIuJ[2][2] += u[0]/tau;


	// Extract u, DIuJ, Find acceleration, theta
	double a[4];
	double theta;
	//AM Checked numerically
	accelerationIdeal(TOut_Full, TOut_BG, i, j, k, dxI, u, DIuJ, a, theta);
	{
		//Reorder lines in the actual evaluations order
		sg[4]=2*DIuJ[0][0]-(2*theta*(1+pow(u[1],2)))/3.+2*a[1]*u[1];
		sg[5]=DIuJ[0][1]+DIuJ[1][0]+a[2]*u[1]+a[1]*u[2]-(2*theta*u[1]*u[2])/3.;
		sg[6]=DIuJ[0][2]+DIuJ[2][0]+a[3]*u[1]+a[1]*u[3]-(2*theta*u[1]*u[3])/3.;
		sg[7]=2*DIuJ[1][1]-(2*theta*(1+pow(u[2],2)))/3.+2*a[2]*u[2];
		sg[8]=DIuJ[1][2]+DIuJ[2][1]+a[3]*u[2]+a[2]*u[3]-(2*theta*u[2]*u[3])/3.;
		sg[9]=2*DIuJ[2][2]-(2*theta*(1+pow(u[3],2)))/3.+2*a[3]*u[3];
		sg[0]=sg[4]+sg[7]+sg[9];
		sg[1]=(sg[4]*u[1]+sg[5]*u[2]+sg[6]*u[3])/u[0];
		sg[2]=(sg[5]*u[1]+sg[7]*u[2]+sg[8]*u[3])/u[0];
		sg[3]=(sg[6]*u[1]+sg[8]*u[2]+sg[9]*u[3])/u[0];
		om[3]=(DIuJ[0][1]-DIuJ[1][0]+a[2]*u[1]-a[1]*u[2])/2.;
		om[4]=(DIuJ[0][2]-DIuJ[2][0]+a[3]*u[1]-a[1]*u[3])/2.;
		om[5]=(DIuJ[1][2]-DIuJ[2][1]+a[3]*u[2]-a[2]*u[3])/2.;
		om[0]=(om[3]*u[2]+om[4]*u[3])/u[0];
		om[1]=(-(om[3]*u[1])+om[5]*u[3])/u[0];
		om[2]=-((om[4]*u[1]+om[5]*u[2])/u[0]);
	}

}


// Get energy and flow at a given position
double GetPData(class EnergyMomentumTensorMap *TOut_Full, class EnergyMomentumTensorMap *TOut_BG, int var_id, int xS, int yS, int zS) {

	//
        void get_energy_flow_from_Tmunu(double Tmunu_raised[4][4], double tau_in_GeV, double * energy_density, double flow[], enum error_codes * error_code);

	//Safety
	if ((xS<0)||(xS>EventInput::Ns-1)||(yS<0)||(yS>EventInput::Ns-1)) return 0.0;

	double Tmunu_Full_raised[4][4];
	double Tmunu_BG_raised[4][4];

	double res=1./0.;

	TOut_Full->GetRaised(xS,yS,&Tmunu_Full_raised[0][0]);
	TOut_BG->GetRaised(xS,yS,&Tmunu_BG_raised[0][0]);

	double energy_density;
	double flow_raised[4];
        enum error_codes error_code;

	get_energy_flow_from_Tmunu(Tmunu_Full_raised, TOut_Full->tau/M_HBARC, & energy_density, flow_raised, &error_code);

	// Energy density
	switch(var_id) {

		//Energy density
		case 0:
			res=energy_density;	
			break;

			//u^\tau
		case 1:
			res=flow_raised[0];
			break;

			//u^x

		case 2:
			res=flow_raised[1];
			break;

			//u^y
		case 3:
			res=flow_raised[2];
			break;

			//tau*u^\eta
		case 4:
			res=0.0*0.0;
			break;
	}

	return res;


}

//Get \sigma^{\mu\nu} in a more convenient format
void get_sigma_proper(class EnergyMomentumTensorMap *TOut_Full, class EnergyMomentumTensorMap *TOut_BG, const int ix, const int iy, double sigma[4][4]) {

	// Forward declaration
	void GetSigma(class EnergyMomentumTensorMap *TOut_Full, class EnergyMomentumTensorMap *TOut_BG, const int i, const int j, const int k, double sg[NPi44], double om[NOmega44]);

	double sigma_array[10];
	double dummy[6];

	// Get sigma
	GetSigma(TOut_Full, TOut_BG, ix, iy, 0, sigma_array, dummy);

	//
	const double tau_in_GeV=TOut_Full->tau/M_HBARC;

	//hbar-c factor necessary because the units of sigma_array is 1/[dx] with dx in fermi
	//factor of 1/2 necessary because it's not included in the definition
	sigma[0][0]=-1./2.*M_HBARC*sigma_array[0];
	sigma[0][1]=-1./2.*M_HBARC*sigma_array[1];
	sigma[0][2]=-1./2.*M_HBARC*sigma_array[2];
	sigma[0][3]=-1./2.*M_HBARC/(tau_in_GeV)*sigma_array[3];
	sigma[1][1]=-1./2.*M_HBARC*sigma_array[4];
	sigma[1][2]=-1./2.*M_HBARC*sigma_array[5];
	sigma[1][3]=-1./2.*M_HBARC/(tau_in_GeV)*sigma_array[6];
	sigma[2][2]=-1./2.*M_HBARC*sigma_array[7];
	sigma[2][3]=-1./2.*M_HBARC/(tau_in_GeV)*sigma_array[8];
	sigma[3][3]=-1./2.*M_HBARC/(tau_in_GeV*tau_in_GeV)*sigma_array[9];

	//Set lower diagonal by symmetry
	sigma[1][0]=sigma[0][1];
	sigma[2][0]=sigma[0][2];
	sigma[3][0]=sigma[0][3];
	sigma[2][1]=sigma[1][2];
	sigma[3][1]=sigma[1][3];
	sigma[3][2]=sigma[2][3];

}

#endif
