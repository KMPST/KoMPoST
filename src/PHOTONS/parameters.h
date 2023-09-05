#ifndef PARAMETERS_H
#define PARAMETERS_H

// Photon degrees of freedom
const int NuGamma = 2;

// Units
const double hbarc=0.197;
const double fmtoGeVm1=1/hbarc;
const double fm2toGeVm2=fmtoGeVm1*fmtoGeVm1;
const double GeVm1tofm=hbarc;
const double GeVm2tofm2=GeVm1tofm*GeVm1tofm;

// Set momentum grid
int const NK=2*99+1;
double const Kmax = 12.0;
double const Kmin = 0.1;
double const dK = (Kmax-Kmin)/(NK-1);

double const alpha_em = 1. / 137.;
const int Nf = 3;
const int NC=3;
const double CF = (NC*NC - 1.0)/(2.0*NC);
const int CA = NC;
const int dF = NC;
const double dA = NC*NC - 1.0;

// Set electric charges of quarks
#define nfamy 3
double const qu=2/3.;
double const qd=-1/3.;
double const qs=-1/3.;

#if (nfamy == 2)
	double const  sum_qf2 =pow(qu,2)+pow(qd,2);
#else
	double const sum_qf2 = pow(qu,2)+pow(qd,2)+pow(qs,2);
#endif

// Include prefactors in order to get spectrum from raw data
double const GammaPrefactor = NuGamma*pow(2*M_PI,-3.)*sum_qf2;

#endif