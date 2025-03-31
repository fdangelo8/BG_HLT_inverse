#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <cmath>
#include <string>
#include <cstdlib>

#include "../include/mp.h"
#include "../include/smear.h"
#include "../include/params.h"
#include "../include/statistical.h"


int main(int argc, char *argv[])
{

  if (argc < 2) 
  {
#if defined(EXP)
    std::cerr << "Compiled with EXP basis" << std::endl;
#endif
#if defined(COS)
    std::cerr << "Compiled with COS basis" << std::endl;
#endif
#if defined(COS_SPHAL)
    std::cerr << "Compiled with COS_SPHAL basis" << std::endl;
#endif
#if defined(BASIS_TEST)
    std::cerr << "Compiled with BASIS_TEST basis" << std::endl;
#endif
#if defined(GAUSS)
    std::cerr << "Compiled with GAUSS target" << std::endl;
#endif
#if defined(PSEUDO_GAUSS)
    std::cerr << "Compiled with PSEUDO_GAUSS target" << std::endl;
#endif
#if defined(COV)
    std::cerr << "Compiled with COV" << std::endl;
#endif
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return 1;  // Exit with an error code if no file is provided
  } 
  cout << "start " << endl;
  char *input_file = argv[1];
  inv_params param;
  read_input(input_file, &param);

  // initialize suplim for integration
#if defined(HLN)
  infLimit=param.E0;
#endif
#if defined(PSEUDO_GAUSS)
  normalization=1.0/Int_PGauss(param.Estar, param.sigma, param.Nt);
#endif

  cout << "E0: " << param.E0 << " infLimit: " << infLimit << endl;
  Real t[param.n_points];
  for(int i=0; i<param.n_points; i++){
    t[i] = i+param.trash;
    cout << "t: " << t[i] << endl;
  }
  cout << "Nt: " << param.n_points << endl;
  
  cout << "********* " << param.Estar  << " *********" << endl;
  cout << "L: " << param.Ns << " #Points: " << param.n_points <<  " Nt: " << param.Nt << " dLim: " << param.dLim << " uLim: " << param.uLim << " trash: " << param.trash << " sigma: " << param.sigma  << " Estar: " << param.Estar << "apar:" << param.apar  << endl;

  //Defining lambda vector (logarithmic scan from 1e-15 to 0.99999)
  for(int ilambda=0; ilambda<Nlambda; ilambda++)
	{
		if(ilambda < 9) lambda(ilambda) = conv(to_string(ilambda+1))*Real("1e-15");
		else if(ilambda >= 9  and ilambda < 18) lambda(ilambda) = conv(to_string(ilambda-9+1))*Real("1e-14");
		else if(ilambda >= 18 and ilambda < 27) lambda(ilambda) = conv(to_string(ilambda-18+1))*Real("1e-13");
		else if(ilambda >= 27 and ilambda < 36) lambda(ilambda) = conv(to_string(ilambda-27+1))*Real("1e-12");
		else if(ilambda >= 36 and ilambda < 45) lambda(ilambda) = conv(to_string(ilambda-36+1))*Real("1e-11");
		else if(ilambda >= 45 and ilambda < 54) lambda(ilambda) = conv(to_string(ilambda-45+1))*Real("1e-10");
		else if(ilambda >= 54 and ilambda < 63) lambda(ilambda) = conv(to_string(ilambda-54+1))*Real("1e-9");
		else if(ilambda >= 63 and ilambda < 72) lambda(ilambda) = conv(to_string(ilambda-63+1))*Real("1e-8");
		else if(ilambda >= 72 and ilambda < 81) lambda(ilambda) = conv(to_string(ilambda-72+1))*Real("1e-7");
		else if(ilambda >= 81 and ilambda < 90) lambda(ilambda) = conv(to_string(ilambda-81+1))*Real("1e-6");
    else if(ilambda >= 90 and ilambda < 99) lambda(ilambda) = conv(to_string(ilambda-90+1))*Real("1e-5");
    else if(ilambda >= 99 and ilambda < 108) lambda(ilambda) = conv(to_string(ilambda-99+1))*Real("1e-4");
    else if(ilambda >= 108 and ilambda < 117) lambda(ilambda) = conv(to_string(ilambda-108+1))*Real("1e-3");
    else if(ilambda >= 117 and ilambda < 126) lambda(ilambda) = conv(to_string(ilambda-117+1))*Real("1e-2");
    else if(ilambda >= 126 and ilambda < 135) lambda(ilambda) = conv(to_string(ilambda-126+1))*Real("1e-1");
    else if(ilambda >= 135 and ilambda < 144) lambda(ilambda) = Real("0.9")+conv(to_string(ilambda-135+1))/100;
    else if(ilambda >= 144 and ilambda < 153) lambda(ilambda) = Real("0.99")+conv(to_string(ilambda-144+1))/1000;
    else if(ilambda >= 153 and ilambda < 162) lambda(ilambda) = Real("0.999")+conv(to_string(ilambda-153+1))/10000;
    else if(ilambda >= 162 and ilambda < 171) lambda(ilambda) = Real("0.9999")+conv(to_string(ilambda-162+1))/100000;
    else if(ilambda >= 171 and ilambda < 180) lambda(ilambda) = Real("0.99999")+conv(to_string(ilambda-171+1))/1000000;
	}

  //BOOTSTRAP
  FILE *Input_Corrs;
  if ((Input_Corrs = fopen(param.boot_corr_file, "r")) == NULL )
  {
    printf("Error opening the input file: %s\n", param.boot_corr_file);
    exit(EXIT_FAILURE);
  } 
  
  PrecMatr Corr_Boot(param.Nboot, param.n_points+param.trash);
  for(int iboot=0; iboot<param.Nboot; iboot++)
  {
    for(int i=0; i<param.n_points+param.trash; i++){     
      int App1, App2;
      char App3[1024], App4[1024], App5[1024];
      fscanf(Input_Corrs, "%d" "%s" "%s" "%s", &App1, App3, App4, App5);
      Corr_Boot(iboot,i) = conv(App4);
      cout << "iboot " << iboot << "  "  << i << "  " <<  Corr_Boot(iboot,i) <<  endl;
    }
  }
  fclose(Input_Corrs);

  //open file with mean and sigma values of correlatore and store it in corr_mu, corr_err
  FILE *Input_Corr_M;
  if ((Input_Corr_M = fopen(param.mean_corr_file, "r")) == NULL )
  {
    printf("Error opening the input file: %s\n", param.mean_corr_file);
    exit(EXIT_FAILURE);
  }
  cout << "Nt+trash " << param.n_points+param.trash << endl;
  PrecVec Corr_err(param.n_points+param.trash), Corr_Mu(param.n_points+param.trash);
  cout << "N: " << param.n_points+param.trash << "  " << param.n_points << "  " << param.trash <<  endl;
  for(int i=0; i<param.n_points+param.trash; i++)
  {
    double a;
    char App3[1024], App4[1024];
    fscanf(Input_Corr_M, "%lf" "%s" "%s", &a, App3, App4);
    Corr_Mu(i) = conv(App3);
    Corr_err(i) = conv(App4);
    cout << "i: " << i << endl; 
    cout << "KKK " << a << " " << Corr_Mu(i) << "  " <<  Corr_err(i) << endl;
  }
  fclose(Input_Corr_M);
  
  //Matrice di covarianza
  PrecMatr Cov(param.n_points, param.n_points);
  for(int i=0; i<param.n_points; i++)
  {
    for(int j=0; j<param.n_points; j++)
    {
      if(i==j)
      {
	      Cov(i,j)=Corr_err(i+param.trash)*Corr_err(i+param.trash);
      }
#if defined(COV)
      else Cov(i,j)=(Corr_Boot.col(i+param.trash).dot(Corr_Boot.col(j+param.trash))/param.Nboot)-Corr_Mu(i+param.trash)*Corr_Mu(j+param.trash);
#else      
      else Cov(i,j)=0;
#endif   
      cout << "i : " << i << " j : " << j << " Cov: " << Cov(i,j) << " Cov_norm: " << Cov(i,j)/(Corr_err(i+param.trash)*Corr_err(j+param.trash)) << endl;
    }
  }
   
  //A,R,f
  cout << "OO " << endl;
  cout << "EEE: " << param.Estar << "  " << param.apar << "  " << t[0] << endl;
  PrecMatr A = A_Comp(t, param.Estar, param.apar, param.Nt, param.n_points, param.E0);
  cout << "HERE " << endl;
  PrecVec R = R_Comp(t, param.n_points, param.Nt);
  PrecVec f;
#if defined(HLN)
  f = f_func(t, param.sigma, param.Estar, param.E0, param.apar, param.n_points, param.Nt, normalization);
  cout << "f(" << param.Estar << ")= " << f << " f(0) = " << f_func(t,param.sigma,0, param.E0, param.apar, param.n_points, param.Nt, normalization) << endl;
#endif
  
  //Coeff
  PrecMatr gl(Nlambda, param.n_points);
  for(int ilambda=0; ilambda<Nlambda; ilambda++)
  {
    gl.row(ilambda) = g_comp(lambda(ilambda), A, Cov, f, R, Corr_Mu(0), param.n_points);
    cout << "gl: " << gl.row(ilambda) << endl;
  }
  
  Real Ver;
  for(int i=0; i<param.n_points; i++)
  {
    Ver += gl(0,i)*Corr_Mu(i+param.trash);
    cout << "EEE: " << gl(0,i) << "  " << Corr_Mu(i+param.trash) << endl;
  }
  cout << "EEE: " << Pi*Ver << endl;
  
  
  //Calcolo integrale quadrato target per funzionale W
  Real TSq = Target_Int(param.apar, param.Estar, param.sigma, normalization, param.Nt);
  

  
  //Calcolo Densità
  PrecMatr Dens(param.Nboot, Nlambda);
  cout << "HERE" << endl;
  for(int iboot=0; iboot<param.Nboot; iboot++)
  {
    for(int ilambda=0; ilambda<Nlambda; ilambda++)
    {
      for(int it=0; it<param.n_points; it++)
      {
       	Dens(iboot, ilambda) += Corr_Boot(iboot, it+param.trash)*gl(ilambda,it);
      }
    }
  }

  PrecVec Dens_Mu(Nlambda), Dens_S(Nlambda);
  
  //Qua mettere un if Estar==0
  for(int ilambda=0; ilambda<Nlambda; ilambda++)
  {
#if defined(COS_SPHAL)
    Dens_Mu(ilambda) = (-1)*2*Pi*Boot_Mean(Dens.col(ilambda), param.Nboot)/(param.Nt*param.Nt);
    Dens_S(ilambda) = 2*Pi*Boot_Sigma(Dens.col(ilambda), param.Nboot)/(param.Nt*param.Nt);
#endif
#if defined(COS) || (BASIS_TEST)
    Dens_Mu(ilambda) = (-4)*Pi*Boot_Mean(Dens.col(ilambda), param.Nboot)/param.Nt;
    Dens_S(ilambda) = 4*Pi*Boot_Sigma(Dens.col(ilambda), param.Nboot)/param.Nt;
#endif
  }
 
  cout << "QUI" << endl;
  Residual_Study(param.rs_vs_lambda, Corr_Mu(0), Cov, gl, t, f, A, param.Estar, param.sigma, param.apar, param.dLim, param.uLim, normalization, param.Nt);
  cout << "SSS1 " << Dens_S(ilambda2) << "  " << Dens_S(ilambda1) << "  " << Dens_Mu(ilambda1) << "  " << Dens_Mu(ilambda2) << endl;
  Real SigmaF = Sigma(Dens_S(ilambda2), Dens_S(ilambda1), Dens_Mu(ilambda1), Dens_Mu(ilambda2));

  // output file reconstruction target function
  smear_reconstruction_output(param.smear_delta_rec_file, t, gl.row(ilambda2), param.Estar, param.sigma, param.Nt, param.n_points, normalization);
  // output file study of spectral density vs lambda
  Print_Study_Lambda(param.rho_vs_lambda_file, Dens_Mu, Dens_S, Corr_Mu(0), Cov, param.Estar, param.sigma, gl, t, SigmaF, param.n_points, param.Nt, param.E0, normalization);
  // final_result
  ofstream file; 
  file.open(param.out_file, std::ios::app); 
  if (!file.is_open()) 
  {
    std::cerr << "Errore: non è stato possibile aprire il file." << std::endl;
    return 1;
  }
  file << param.sigma << "\t" << Dens_Mu(ilambda2) << "\t" << SigmaF << endl;
  file.close();

  return 0;
}