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
#ifndef EnergyMomentumDecomposition_cpp
#define EnergyMomentumDecomposition_cpp

#include <iostream>
#include <cstdlib>
#include <vector>
#include <array>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_poly.h>
#include <typeinfo>
#include <sstream>

enum error_codes { decomp_ok=0, decomp_failed=2};

// Get shear tensor \pi^{\mu\nu} and bulk pressure \Pi given T^\mu_\nu, epsilon
// and u^\mu at one spacetime point. The energy-momentum tensor of that point is given by Tmunu_raised[4][4]
// Assume conformal EOS
void get_shear_and_bulk_from_Tmunu(double Tmunu_raised[4][4], double tau_in_GeV, double energy_density, double * flow, double pimunu_raised[4][4], double & bulk_pressure) {

	//Forward declarations

	//For convenience
	const double Ttautau=Tmunu_raised[0][0];
	const double Txx    =Tmunu_raised[1][1];
	const double Tyy    =Tmunu_raised[2][2];
	const double Tetaeta=Tmunu_raised[3][3];
	const double gmunu_raised_hyperb[4][4]={{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1./(tau_in_GeV*tau_in_GeV)}};

	//Assume a conformal system
	double pressure=energy_density/3.;
	//Trace of Tmunu
	double Tmunu_trace=Ttautau-Txx-Tyy-tau_in_GeV*tau_in_GeV*Tetaeta;
	//Bulk pressure
	double tmp_bulk_pressure=(energy_density/3.-pressure)-Tmunu_trace/3.;

	bulk_pressure=tmp_bulk_pressure;

	// flow lowered
	double flow_lowered[4]={flow[0],-1*flow[1],-1*flow[2],-tau_in_GeV*tau_in_GeV*flow[3]};

   // Shear tensor: pi^{\mu\nu}=T^{\mu\nu}+epsilon*u^\mu*u^\nu-T^{\mu \alpha}
   // u_{\alpha} u^\nu-T^{\nu \alpha} u_{\alpha} u^\mu-1/3*(g^{\mu\nu}-u^\mu
   // u^\nu)*(T^\alpha_\alpha - epsilon) or: pi^{\mu\nu}=T^{\mu\nu}-e u^\mu
   // u^\nu+(g^{\mu\nu}-u^\mu u^\nu) (P+Pi)
	for(int mu=0;mu<4;mu++) {
		for(int nu=0;nu<4;nu++) {
         // method 1 (insured orthogonality to flow):
         // pi^{\mu\nu}=T^{\mu\nu}+epsilon*u^\mu*u^\nu-T^{\mu \alpha}
         // u_{\alpha} u^\nu-T^{\nu \alpha} u_{\alpha}
         // u^\mu-1/3*(g^{\mu\nu}-u^\mu u^\nu)*(T^\alpha_\alpha - epsilon)

			double tmp=Tmunu_raised[mu][nu]+energy_density*flow[mu]*flow[nu]-1./3.*(gmunu_raised_hyperb[mu][nu]-flow[mu]*flow[nu])*(Tmunu_trace-energy_density);

			for(int alpha=0;alpha<4;alpha++) tmp-=Tmunu_raised[mu][alpha]*flow[nu]*flow_lowered[alpha];
			for(int alpha=0;alpha<4;alpha++) tmp-=Tmunu_raised[nu][alpha]*flow[mu]*flow_lowered[alpha];

         // method 2:
         // pimunu_raised[mu][nu]=Tmunu_raised[mu][nu]-energy_density*flow[mu]*flow[nu]+(gmunu_raised_hyperb[mu][nu]-flow[mu]*flow[nu])*(pressure+tmp_bulk_pressure);
			// const double tmp2=Tmunu_raised[mu][nu]-energy_density*flow[mu]*flow[nu]+(gmunu_raised_hyperb[mu][nu]-flow[mu]*flow[nu])*(pressure+tmp_bulk_pressure);

			pimunu_raised[mu][nu]=tmp;
		}
	}
}

//  Solves T^\mu_\nu u^\nu = epsilon u^\mu for one spacetime point
//  The only timelike eigenvalue is the flow u^\mu (when normalised to unity) and
//  to the corresponding eigenvalue is the energy density epsilon 
//  The flow is unitless and the energy density has the units of T^munu
//  Input are Tmunu (Tmunu_raised[4][4]) and time (tau_in_GeV)
//  Outputs are the energy_density, the flow and an error code (decomp_failed if no rest frame could be found, decomp_OK if epsilon and flow were found)
void get_energy_flow_from_Tmunu(double Tmunu_raised[4][4], double tau_in_GeV, double * energy_density, double flow[], enum error_codes * error_code) {

	//GSL stuff
	gsl_eigen_nonsymmv_workspace * eigen_ws = gsl_eigen_nonsymmv_alloc(4);
	const int balance_matrix=1;
	gsl_eigen_nonsymmv_params(balance_matrix, eigen_ws);
	gsl_vector_complex * eigenvalues= gsl_vector_complex_alloc (4);
	gsl_matrix_complex * eigenvectors=gsl_matrix_complex_alloc(4,4);


	// Make proper 4X4 matrix from Tmunu (depending on the raised and lower indices
	// / T^\tau_\tau   T^\tau_x    T^\tau_y    T^\tau_\eta \       / T^\tau\tau    -T^\tau x    -T^\tau y    -tau^2*T^\tau \eta
	// | T^x_\tau      T^x_x       T^x_y       T^x_\eta    | ==== |  T^x \tau      -T^x x       -T^x y       -tau^2*T^x \eta    |
	// | T^y_\tau      T^y_x       T^y_y       T^y_\eta    | ==== |  T^y \tau      -T^x y       -T^y y       -tau^2*T^y \eta    |
	// \ T^\eta_\tau   T^\eta_x    T^\eta_y    T^\eta_\eta /       \ T^\eta\tau    -T^\eta x    -T^\eta y    -tau^2*T^\eta\eta /


	const double Ttautau=Tmunu_raised[0][0];
	const double Ttaux  =Tmunu_raised[0][1];
	const double Ttauy  =Tmunu_raised[0][2];
	const double Ttaueta=Tmunu_raised[0][3];
	const double Txx    =Tmunu_raised[1][1];
	const double Txy    =Tmunu_raised[1][2];
	const double Txeta  =Tmunu_raised[1][3];
	const double Tyy    =Tmunu_raised[2][2];
	const double Tyeta  =Tmunu_raised[2][3];
	const double Tetaeta=Tmunu_raised[3][3];

	double tmpdata[]={  Ttautau,   -Ttaux,  -Ttauy, -tau_in_GeV*tau_in_GeV*Ttaueta,
	                    Ttaux,     -Txx,    -Txy,   -tau_in_GeV*tau_in_GeV*Txeta,
	                    Ttauy,     -Txy,    -Tyy,   -tau_in_GeV*tau_in_GeV*Tyeta,
	                    Ttaueta,   -Txeta,  -Tyeta, -tau_in_GeV*tau_in_GeV*Tetaeta}; 	                  
        // GSL magic
	gsl_matrix_view tmp_matrix = gsl_matrix_view_array(tmpdata, 4, 4);
	gsl_matrix * matrix = gsl_matrix_alloc (4, 4);
	gsl_matrix * matrix_orig = gsl_matrix_alloc (4, 4);
	*matrix=tmp_matrix.matrix;
	gsl_matrix_memcpy(matrix_orig,&tmp_matrix.matrix);

	
	// Solve eigensystem 
	//gsl_eigen_nonsymmv (gsl_matrix * A, gsl_vector_complex * eval, gsl_matrix_complex * evec, gsl_eigen_nonsymmv_workspace * w);
	gsl_eigen_nonsymmv(matrix, eigenvalues, eigenvectors, eigen_ws);

	
	//If a complex number has fabs(imag/real) smaller than this number, we consider if real
	//(to accommodate for possible tiny numerical error)
	const double real_enough_ratio=1e-30;

	// Find the only timelike eigenvector (flow) and corresponding eigenvalue (energy)
	bool parse_eigensystem(double tau_in_GeV, double im_tolerance, gsl_matrix * matrix_orig, gsl_vector_complex * eigenvalues, gsl_matrix_complex * eigenvectors, double * energy, double flow[]);
	bool decomp_success=parse_eigensystem(tau_in_GeV, real_enough_ratio, matrix_orig,eigenvalues,eigenvectors,energy_density, flow);
	enum error_codes local_error_code=decomp_failed;

	if (decomp_success) {
		local_error_code=decomp_ok;
	}

	*error_code=local_error_code;


	//GSL stuff
	gsl_matrix_free(matrix);
	gsl_matrix_free(matrix_orig);
	gsl_vector_complex_free(eigenvalues);
	gsl_matrix_complex_free(eigenvectors);
	gsl_eigen_nonsymmv_free(eigen_ws);

}


// Take eigenvalues and eigenvectors of a matrix and identify the eigenvalue with a positive real part, a small imaginary part given by parameter "im_tolerance", and the corresponding eigenvector
// An exact solution correspond to a single real positive eigenvalue and a timelike vector
// Need tau for scalar product, since all of this is happening with hyperbolic coordinates
// Return: true for success (a single solution is found), false for failure
bool parse_eigensystem(double tau_in_GeV, double im_tolerance, gsl_matrix * matrix_orig, gsl_vector_complex * eigenvalues, gsl_matrix_complex * eigenvectors, double * energy, double flow[]) {

	bool found_a_candidate=false;
	bool found_too_many_candidates=false;
	double tmp_energy=-1.;
	double tmp_flow[4]={0.0};
	double eigenvalue_candidate;
	double eigenvector_candidate[4]={0.0};
	double eigenvector_candidate_norm;

	//Loop over eigenvalues and eigenvectors
	for(int i=0;i<4;i++) {

		//Eigenvalue
		gsl_complex tmp_eigenval=gsl_vector_complex_get(eigenvalues,i);
		const double tmp_eigenval_real=GSL_REAL(tmp_eigenval);
		const double tmp_eigenval_imag=GSL_IMAG(tmp_eigenval);

		//Only continue if the eigenvalue is real and positive
		if ((tmp_eigenval_real > 0)&&(fabs(tmp_eigenval_imag/tmp_eigenval_real) < im_tolerance)) {

			//Loop over elements of the eigenvector
			bool all_real=true;
			double norm=0.0;
			double eigenvector[4];
			for(int j=0;j<4;j++) {

				//get eigenvector element
				gsl_complex tmp_eigenvec = gsl_matrix_complex_get(eigenvectors, j, i);
				const double tmp_eigenvec_real=GSL_REAL(tmp_eigenvec);
				const double tmp_eigenvec_imag=GSL_IMAG(tmp_eigenvec);

				//weight from g_{\mu\nu} in hyperbolic coordinates: (1,-1,-1,-tau^2)
				double weight;
				switch(j) {
					//time component
					case 0:
						weight=1;
						break;

						//x, y
					case 1:
					case 2:
						weight=-1;
						break;
						//eta
					case 3:
						weight=-tau_in_GeV*tau_in_GeV;
						break;
				}	

				//is it real?
				if ((0.0 == tmp_eigenvec_imag)||(fabs(tmp_eigenvec_imag/tmp_eigenvec_real) < im_tolerance)) {
					//compute (Minkowski) norm
					norm+=weight*tmp_eigenvec_real*tmp_eigenvec_real;
					eigenvector[j]=tmp_eigenvec_real;
				}	
				else {
					all_real=false;
					break;
				}
			}
			//If the vector is timelike and all its components are real, jackpot!...
			if ((norm > 0)&&(all_real)) {
				//...except if another real timelike candidate was already found, in which case we're in trouble
				if (found_a_candidate) { 
					const double similar_enough=1e-3;
					//... although if the other candidate differs only by the imaginary part, that's OK
                                        // I know this is not great coding
					if ((fabs(1-eigenvalue_candidate/tmp_eigenval_real) > similar_enough)&&(eigenvalue_candidate != tmp_eigenval_real)&&(fabs(1-eigenvector_candidate[0]/eigenvector[0]) > similar_enough)&&(eigenvector_candidate[0] != eigenvector[0])&&(fabs(1-eigenvector_candidate[1]/eigenvector[1]) > similar_enough)&&(eigenvector_candidate[1] != eigenvector[1])&&(fabs(1-eigenvector_candidate[2]/eigenvector[2]) > similar_enough)&&(eigenvector_candidate[2] != eigenvector[2])&&(fabs(1-eigenvector_candidate[3]/eigenvector[3]) > similar_enough)&&(eigenvector_candidate[3] != eigenvector[3])) {
						found_too_many_candidates=true;
					}
					//std::cout << "Too many solutions!\n";
				}

				eigenvalue_candidate=tmp_eigenval_real;
				eigenvector_candidate_norm=norm;

				for(int k=0;k<4;k++) eigenvector_candidate[k]=eigenvector[k];

				found_a_candidate=true;
			}
		}
	}

	tmp_energy=eigenvalue_candidate;

	//Assume the time-component is never zero. If it is, there's something weird going on
	const double time_component_sign=eigenvector_candidate[0]/fabs(eigenvector_candidate[0]);
	for(int k=0;k<4;k++) tmp_flow[k]=eigenvector_candidate[k]*time_component_sign/(sqrt(eigenvector_candidate_norm));

	bool success=(!found_too_many_candidates)&&(found_a_candidate);

	//If there's a problem, output debug info 
	if (success) {

			*energy=tmp_energy;

			for(int i=0;i<4;i++) flow[i]=tmp_flow[i];

	}

	return  success;

}
#endif
