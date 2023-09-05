This directory contains data from the EKT simulation described in https://arxiv.org/abs/2308.09747.

In "ConstantsLambda10.txt" you find constants regarding the EKT simulation for \lambda=10. Namely, it contains (\tau^(1/3)T)_\infty, eta/s and the normalization
constant CIdeal = 4.0*\int z^4 \overline{C}_eq(z). Note that CIdeal itself does not contain the sum of squared charges here and needs to be multiplied by it.

"EvolutionParametersLambda10.txt" provides a table with which the stepcounter can be translated into the eigentime \tau, the dimensionless time
variable \tilde{w} and the current temperature determined from Landau-matching for \nu_{eff}=\nu_G+2*7/8*N_f*\nu_Q (and N_f=3).

The files "LAMBDA10/pT_SPECTRUM_LEADING_ORDER_RATES/pTSpectrum_LeadingOrderRatePhotons*.txt" contain data about the pT spectrum of the EKT simulation for \lambda=10 (\lambda = g^2*Nc). Each file has
the following structure:
1:pT 2:(dN/d^2p_T d^2x_T dy)/(NuGamma/(2.0*pi)**3*(Qu**2+Qd**2+Qs**2))
It is important to note that the value in column 2 needs to be multiplied by NuGamma/(2.0*pi)**3*(Qu**2+Qd**2+Qs**2) in order to get the actual rate. pT is given 
in units of a Q_s. To translate this into physical units please see XXX.
The Number in the files (e.g. 550 in "LAMBDA10/pTSpectrum_LeadingOrderRatePhotons550.txt") corresponds to the stepcounter (i.e. time step).
Using "EvolutionParametersLambda10.txt" the stepcounter can be translated to other evolution parameters as described above.