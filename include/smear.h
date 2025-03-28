#ifndef _SMEAR_H
#define _SMEAR_H

#include "mp.h"

extern Real normalization;

//target functions prototypes
Real Z( Real s, Real Es);
Real Int_PGauss(Real Es, Real s, Real beta);
Real Target_F(Real E, Real Es, Real s, Real norm, Real beta);
// basis function prototype
Real K(Real omega, Real t, Real beta);
// DEFINITIONS OF THE FUNCTIONS USED FOR THE RESOLUTION METHOD
Real Delta_Smear(Real omega, PrecVec q, Real t_in[], Real beta, int Nt);
Real Target_Int(double ap, Real Estar, Real sigma, Real norm, Real beta);
PrecVec Coeff(PrecVec R, PrecMatr Winv, PrecVec f, Real l);
Real spectral(PrecVec q, PrecVec C, int Nt);
Real stat_unc(PrecVec q, PrecVec dC, int Nt);
// COMPUTATION A, R, f (INDEPENDENT FROM \lambda AND CORRELATORS)
// Computation A
PrecMatr A_Comp(Real t_a[], double Estar, double ap, Real beta, int Nt, Real E0);
// Computation R
PrecVec R_Comp(Real t_a[], int Nt, Real beta);
// Computation f 
Real f_NInt(Real infLimit, Real supLimit, Real ti, Real s, Real Es, double ap, Real beta, Real norm);
Real N(Real t, Real s, Real Es);
Real D(Real t, Real s, Real Es, Real E0);
PrecVec f_func(Real t_a[], Real s, Real Es, Real E0, double ap, int Nt, Real beta, Real norm);
//Computation W e coefficients
PrecVec g_comp(Real lambda, PrecMatr A, PrecMatr Cov, PrecVec f, PrecVec R, Real Corr0, int Nt);
//Computation functional a posteriori (global coefficients)
Real W_func_comp(Real lambda, Real Corr_zero, PrecMatr Cov, PrecVec g, double ap, PrecMatr A, PrecVec f, Real Estar, Real sigma, Real norm, Real beta);
// SOME OUTPUT FUNCTIONS
//Output Smearing function
void smear_reconstruction_output(char open_Delta_S[], Real t_a[], PrecVec g, double Es, Real s, Real beta, int Nt, Real norm);
void Residual_Study(char out_res_study[], Real Corr_Z, PrecMatr Cov, PrecMatr gl, Real t[], PrecVec f, PrecMatr A, Real Estar, Real sigma, double apar, double dLim, double uLim, Real norm, Real beta);
Real Sigma(Real Sigma_Stat2, Real Sigma_Stat1, Real M_1, Real M_2);
void Print_Study_Lambda(char output_file[], PrecVec Dens_Mu, PrecVec Dens_S, Real Corr_Z, PrecMatr Cov, Real Estar, Real sigma, PrecMatr gl, Real t[], Real SigmaF, int Nt, Real beta, Real E0, Real norm);

#endif