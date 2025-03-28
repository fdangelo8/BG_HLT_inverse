#ifndef PARAMS_H
#define PARAMS_H

#include "mp.h"

#define STD_STRING_LENGTH 100 // standarg length of unknown strings

// set desired precision in bits                                  
const int P = 1024;
struct Initer
{ 
  Initer()
  {
    mpfr_class::set_dprec(P);
  }
};

extern Initer initer;

typedef struct inv_params {

	// lattice dimensions
	int Ns;  //Ns
	int Nt; //Nt, estensione totale del reticolo

	// paths input files
	char boot_corr_file[STD_STRING_LENGTH];
	int Nboot;
	char mean_corr_file[STD_STRING_LENGTH];

	// parameters for inversion
	int n_points; //number of points used for the inversion
	int trash; //number of points excluded, starting from t=0
	double Estar; //point where inversion is performed
	Real sigma;
	double dLim;
	double uLim;
	double apar; //parametro alpha nel funzionale
	Real E0; //punto da cui si inizia a integrare nel funzionale

	// output file names
	char out_file[STD_STRING_LENGTH]; //per il file out_ens_sigma
	char smear_delta_rec_file[STD_STRING_LENGTH]; //per il file delta_smear
	char rho_vs_lambda_file[STD_STRING_LENGTH]; //per lo studio di rho vs lambda
	char rs_vs_lambda[STD_STRING_LENGTH]; //corrisponde al file RS

} inv_params;

// global parameters for lambda array
extern const int Nlambda;
extern PrecVec lambda;
extern Real lambda1, lambda2;
extern int ilambda1, ilambda2;

// parameters for numerical integration
#if defined(BG)
const Real infLimit=0;
#endif
#if defined(HLN)
extern Real infLimit;
#endif
const Real supLimit=inf;			

//other useful quantity (I don't know what it is, it is needed when one has the simple exponential basis, this is not our case)
static const Real alpha=0; 
                      
extern string ens;

void remove_white_lines_and_comments(FILE *input);
void read_input(char *in_file, inv_params *param);

#endif