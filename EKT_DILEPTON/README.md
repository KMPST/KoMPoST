This directory contains data from the EKT simulation described in https://arxiv.org/abs/2403.04846.

In "ConstantsLambda10.txt" you find constants regarding the EKT simulation for $\lambda=10$. Namely, it contains $(\tau^{1/3}T)_\infty$, $\eta/s$.

"EvolutionParametersLambda10.txt" provides a table with which the stepcounter can be translated into the eigentime $\tau$, the dimensionless time
variable $\tilde{w}$ and the current temperature determined from Landau-matching for $\nu_{eff}=\nu_G+2*7/8*N_f*\nu_Q$ (and $N_f=3$).

The files "LAMBDA10/YIELD_DILEPTONS/Yield_dNQdQ_LO_*.txt" contain data about the M spectrum of the EKT simulation for $\lambda=10$ ($\lambda = g^2*N_c$). Each file has
the following structure:
1:$M$ 2:$dN/MdM dy_q d2x_T/(Q_{\text{up}}^2+Q_{\text{down}}^2+Q_{\text{strange}}^2)$
It is important to note that the value in column 2 needs to be multiplied by $(Q_{\text{up}}^2+Q_{\text{down}}^2+Q_{\text{strange}}^2)$ in order to get the actual rate. $M$ is given 
in units of a $Q_s$. To translate this into physical units please see https://arxiv.org/abs/2403.04846.
The Number in the files (e.g. 550 in "LAMBDA10/YIELD_DILEPTONS/Yield_dNQdQ_LO_550.txt") corresponds to the stepcounter (i.e. time step).
Using "EvolutionParametersLambda10.txt" the stepcounter can be translated to other evolution parameters as described above.