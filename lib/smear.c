#include "../include/smear.h"
#include "../include/params.h"
#include "../include/statistical.h"

#include <string.h>

//Gaussian integral
Real Z(Real s, Real Es)
{
  return (1+erf(Es/(sqrt(2)*s)))/2;
}

//Pseudogauss integral
Real Int_PGauss(Real Es, Real s, Real beta)
{
  const auto f_cs=
    [=](const Real& E) -> Real
    {
#if defined(COS_SPHAL)
      return (E-Es)/(sinh((E-Es)/s));
#endif
#if defined(COS)
      return ((E-Es)/(sinh((E-Es)/s)))*(E/(2*(1-exp(-beta*E))));
#endif
    };
  return bq::gauss_kronrod<Real,61>::integrate(f_cs,0.0,supLimit,5,1e-16);

}

Real normalization = 1.0;

//Target function
Real Target_F(Real E, Real Es, Real s, Real norm, Real beta)
{
  Real ret;
#if defined(GAUSS)
  ret = exp(-pow(E-Es,2)/(2*s*s))/(sqrt(2*Pi)*s*Z(s,Es));
#endif
#if defined(PSEUDO_GAUSS)
#if defined(COS_SPHAL)
  ret = (E-Es)/(sinh((E-Es)/s)) * norm;
#endif
#if defined(COS)
  ret = ((E-Es)/(sinh((E-Es)/s)))*(E/(2*(1-exp(-beta*E))))*norm;
#endif
#endif
  return ret;
}

//Basis function
Real K(Real omega, Real t, Real beta)
{
  Real ret;
#if defined(EXP)
  ret=exp(-omega*t); //+ exp(-(beta-t)*omega);
#endif 
#if defined(COS)
  //Real A=cosh(omega*(t - beta/2));
  //Real B=sinh(beta*omega/2);
  //ret = A;
  ret = exp(-omega*t)+ exp(-(beta-t)*omega);
#endif
#if defined(COS_SPHAL)
  Real A=cosh(omega*(t-beta/2)); //basis for sphal rate 
  Real B=sinh(beta*omega/2);
  ret = A/B*omega;
#endif
#if defined(BASIS_TEST)
  Real A=(1-exp((-1)*beta*omega))*cosh(omega*(t - beta/2));  // for tests
  Real B=2.0*sinh(beta*omega/2);
  ret = A/B;
#endif
  return ret; 
}

//Smearing function
Real Delta_Smear(Real omega, PrecVec q, Real t_in[], Real beta, int Nt)
{  
  Real D;
  for(int i=0; i<Nt; i++)
  {
    //if(omega==0) D += q(i)*2/beta;
    if(omega==0) D += q(i);
    else D += q(i)*K(omega, t_in[i], beta);
  }
  return D; 
}

//Target funct integral
Real Target_Int(double ap, Real Estar, Real sigma, Real norm, Real beta)
{
  const auto f_cs=
    [=](const Real& E) -> Real
    {
      return Target_F(E, Estar, sigma, norm, beta)*Target_F(E, Estar, sigma, norm, beta)*exp(ap*E);
    };
  return bq::gauss_kronrod<Real,61>::integrate(f_cs,infLimit,supLimit,5,1e-16);
}

//Coefficients computation
PrecVec Coeff(PrecVec R, PrecMatr Winv, PrecVec f, Real l)
{
  Real den =  R.transpose()*Winv*R;
#if defined(HLN)
  Real numA = R.transpose()*Winv*(1-l)*f;
  Real num = 1-numA;
  return Winv*(1-l)*f+ Winv*R*num/den;
  //return Winv*(1-l)*f;
#endif
#if defined(BG)
  return Winv*R/den;
#endif
}

//Spectral function from coefficients
Real spectral(PrecVec q, PrecVec C, int Nt)
{
  Real rho=0;
  for(int i=0; i<Nt; i++)
  {
    rho += q(i)*C(i);
  }
  return Pi*rho;
}

//Naive statistical uncertainty
Real stat_unc(PrecVec q, PrecVec dC, int Nt)
{
  Real err=0;
  for(int i=0; i<Nt; i++)
  {
    err += q(i)*dC(i);
    cout << "q: " << q(i) << "  "  << dC(i) << "  " << err << endl;
  } 
  return Pi*err;
}

//Output Smearing function
void smear_reconstruction_output(char open_Delta_S[], Real t_a[], PrecVec g, double Es, Real sigma, Real beta, int Nt, Real norm)
{
  FILE *Delta_S;
  if ((Delta_S = fopen(open_Delta_S, "w")) == NULL)
  {
    printf("Error opening the input file: %s\n",open_Delta_S);
    exit(EXIT_FAILURE);
  }
#if defined(BG)
  for(double i=0; i<300; i++)
    fprintf(Delta_S, "%s\t%s\n", conv(0.0001 + i/100).c_str(), conv(Delta_Smear(0.0001 + i/100, g, t_a, Nt)).c_str());
#endif
#if defined(HLN)
  fprintf(Delta_S, "#omega\tdelta_rec\tTarget_F\tTarget_F-delta_rec\n");
  // Energy vector creation for the reconstruction plot, 500 equally spaced points between Estar-5*sigma and Estar+5*sigma
  double s = sigma.convert_to<double>();
  double half_width = 8*s;
  double lower_bound = Es - half_width;
  double upper_bound = Es + half_width;
  int N = 800;
  double step = (upper_bound - lower_bound) / (N - 1);
  for(int i = 0; i < N; i++)
  {
    Real energy    = lower_bound + i * step;
    Real target    = Target_F(energy, Es, s, norm, beta);
    Real smear_rec = Delta_Smear(energy, g, t_a, beta, Nt);
    Real df        = target - smear_rec;
    fprintf(Delta_S, "%s\t%s\t%s\t%s\n", conv(energy).c_str(), conv(smear_rec).c_str(), conv(target).c_str(), conv(df).c_str());
  }
#endif 
  fclose(Delta_S);
}
  
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Computation A, R, f (Independent from \lambda and correlators)
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//Computation A
PrecMatr A_Comp(Real t_a[], double Estar, double ap, Real beta, int Nt, Real E0)
{
  PrecMatr A(Nt,Nt);
  for(int i=0; i<Nt; i++)
  {
    for(int j=0; j<Nt; j++)
    {
    //With exp basis 
#if defined(EXP)
#if defined(HLN)
      A(i,j) = exp(-(t_a[i]+t_a[j]-alpha)*E0)/(t_a[i]+t_a[j]-alpha);
#endif
#if defined(BG)
      A(i,j) = -(-2+2*Estar*(t_a[i]+t_a[j])-pow(Estar,2)*pow(t_a[i]+t_a[j],2))/pow(t_a[i]+t_a[j],3);
#endif      
#endif
//With cos basis integrate
#if defined(COS) || (COS_SPHAL) || (BASIS_TEST)
      const auto f_cs=
	[=](const Real& E) -> Real
	{
#if defined(BG)
	  return K(E, t_a[i], beta)*(E - Estar)*(E - Estar)*K(E,t_a[j],beta);
#endif
#if defined(HLN)
	  return K(E, t_a[i], beta)*K(E,t_a[j],beta)*exp(ap*E); 
#endif
	};
      A(i,j) = bq::gauss_kronrod<Real,61>::integrate(f_cs,infLimit,supLimit,5,1e-16);
#endif
       
    }//j
  }//i
  cout << "A: " << A << endl;
  return A;
} 

//Computation R
PrecVec R_Comp(Real t_a[], int Nt, Real beta)
{
  PrecVec R(Nt);
  for(int i=0; i<Nt; i++)
  {
#if defined(EXP)
#if defined(HLN)
    R(i) = 1/(t_a[i])*exp(-E0*t_a[i]);
#endif
#if defined(BG)
    R(i) = 1/(t_a[i]);
#endif
#endif
#if defined(COS) || (COS_SPHAL) || (BASIS_TEST)
    const auto f_cs=
    [=](const Real& E) -> Real
    {
     return K(E, t_a[i], beta);
    };
    R(i) = bq::gauss_kronrod<Real,61>::integrate(f_cs,0.0,supLimit,5,1e-16);
#endif
  }
  return R;
}

//Computation f
Real f_NInt(Real infLimit, Real supLimit, Real ti, Real s, Real Es, double ap, Real beta, Real norm)
{
  const auto f_cs=
    [=](const Real& E) -> Real
    {
     return K(E, ti, beta)*Target_F(E,Es,s,norm,beta)*exp(ap*E);
    };
  return
    bq::gauss_kronrod<Real,61>::integrate(f_cs,infLimit,supLimit,5,1e-16);
  
}

Real N(Real t, Real s, Real Es)
{
  return 1/(2*Z(s,Es))*exp(((alpha-t)*((alpha-t)*pow(s,2)+2*Es))/2);
}

Real D(Real t, Real s, Real Es, Real E0)
{
  return  1+erf(((alpha-t)*pow(s,2)+Es-E0)/(sqrt(2)*s)); 
}

PrecVec f_func(Real t_a[], Real s, Real Es, Real E0, double ap, int Nt, Real beta, Real norm)
{
  PrecVec f(Nt);
  for(int i=0; i<Nt; i++)
  {
#if defined(EXP)
    f(i) = N(t_a[i], s, Es)*D(t_a[i], s, Es, E0);
#endif
#if defined(COS) || (COS_SPHAL) || (BASIS_TEST)
    f(i) = f_NInt(infLimit, supLimit, t_a[i], s, Es, ap, beta, norm);
#endif
  }
  return f;
}

//Computation W e coefficients
PrecVec g_comp(Real lambda, PrecMatr A, PrecMatr Cov, PrecVec f, PrecVec R, Real Corr0, int Nt)
{
  PrecMatr W_Mat(Nt,Nt);
  W_Mat = (1-lambda)*A + lambda*Cov/(pow(Corr0,2));
  //Inversion matrix W
  const auto Winv=W_Mat.inverse();
  //Computation g
  return Coeff(R,Winv,f,lambda);
}

//Computation functional a posteriori (global coefficients)
Real W_func_comp(Real lambda, Real Corr_zero, PrecMatr Cov, PrecVec g, double ap, PrecMatr A, PrecVec f, Real Estar, Real sigma, Real norm, Real beta)
{
  Real A1=g.transpose()*A*g;
  Real A2=-2*f.transpose()*g;
  Real A3=Target_Int(ap, Estar, sigma, norm, beta);
  Real B=lambda/(Corr_zero*Corr_zero)*g.transpose()*Cov*g;
  return (1-lambda)*(A1+A2+A3)+B;
}

void Residual_Study(char out_res_study[], Real Corr_Z, PrecMatr Cov, PrecMatr gl, Real t[], PrecVec f, PrecMatr A, Real Estar, Real sigma, double apar, double dLim, double uLim, Real norm, Real beta, int Nt, Real E0)
{
  FILE *RS;
  if ((RS = fopen(out_res_study, "w")) == NULL)
  {
    printf("Error opening the input file: %s\n", out_res_study);
    exit(EXIT_FAILURE);
  }

  double Estar_double = Estar.convert_to<double>();
  PrecMatr A0 = A_Comp(t, Estar_double, 0, beta, Nt, E0);
  PrecVec f0 = f_func(t, sigma, Estar, E0, 0, Nt, beta, norm);

  int a=0;
  for(int ilambda=0; ilambda<Nlambda; ilambda++)
  {
    fprintf(RS, "%s\t%s\t%s\t%s\n", conv(lambda(ilambda)).c_str(), conv(W_func_comp(0, Corr_Z, Cov, gl.row(ilambda), apar, A, f, Estar, sigma, norm, beta)/Target_Int(apar, Estar, sigma, norm, beta)).c_str(), conv(W_func_comp(1, Corr_Z, Cov, gl.row(ilambda), apar, A, f, Estar, sigma, norm, beta)).c_str(), conv((W_func_comp(1, Corr_Z, Cov, gl.row(ilambda), apar, A, f, Estar, sigma, norm, beta)/(W_func_comp(0, Corr_Z, Cov, gl.row(ilambda), apar, A, f, Estar, sigma, norm, beta))*Target_Int(apar, Estar, sigma, norm, beta))).c_str());
    if(dLim<sqrt(W_func_comp(0, Corr_Z, Cov, gl.row(ilambda), 0, A0, f0, Estar, sigma, norm, beta)/Target_Int(0, Estar, sigma, norm, beta)) and a==0)
    {
      lambda1 = lambda(ilambda-1);
      ilambda1=ilambda;
      a++;
    }
    if(uLim<sqrt(W_func_comp(0, Corr_Z, Cov, gl.row(ilambda), 0, A0, f0, Estar, sigma, norm, beta)/Target_Int(0, Estar, sigma, norm, beta)) and a==1)
    {
      lambda2 = lambda(ilambda-1);
      ilambda2=ilambda;
      a++;
    }
    cout << " ilambda : " << ilambda << "\tlambda : " << lambda(ilambda) <<endl;
  }

  cout << "lambda1: " << lambda1 << " lambda2: " << lambda2 << endl; 
  fclose(RS); 
}

Real Sigma(Real Sigma_Stat2, Real Sigma_Stat1, Real M_1, Real M_2)
{
  Real Ps = M_1-M_2/Sigma_Stat1;
  Real syst = abs(M_1-M_2)*erf(abs(Ps)/sqrt(2));
  return sqrt(Sigma_Stat2*Sigma_Stat2 + syst*syst); 
}

void Print_Study_Lambda(char output_file[], PrecVec Dens_Mu, PrecVec Dens_S, Real Corr_Z, PrecMatr Cov, Real Estar, Real sigma, PrecMatr gl, Real t[], Real SigmaF, int Nt, Real beta, Real E0, Real norm)
{ 
  FILE *Lambda_Shape_out;
  if ((Lambda_Shape_out = fopen(output_file, "w")) == NULL)
  {
    printf("Error opening the input file: %s\n", output_file);
    exit(EXIT_FAILURE);
  }   
  
  //Output rho in funzione di lambda
  double Estar_double = Estar.convert_to<double>();
  PrecMatr A0 = A_Comp(t, Estar_double, 0, beta, Nt, E0);
  PrecVec f0 = f_func(t, sigma, Estar, E0, 0, Nt, beta, norm);
  
  fprintf(Lambda_Shape_out, "#lambda\tA[g]/Tar_Int\td[g_t]\tdens_mu(lambda)\tdens_err(lambda)\n");
  for(int ilambda=0; ilambda<Nlambda; ilambda++)
  {
    cout << "lambda: " << lambda(ilambda) << "\tA[g]/Tar_Int: " << sqrt(W_func_comp(0,Corr_Z, Cov, gl.row(ilambda), 0, A0, f0, Estar, sigma, norm, beta)/Target_Int(0, Estar, sigma, norm, beta)) << "\td[g]: " << sqrt(W_func_comp(0,Corr_Z, Cov, gl.row(ilambda), 0, A0, f0, Estar, sigma, norm, beta)/W_func_comp(1, Corr_Z, Cov, gl.row(ilambda), 0, A0, f0, Estar, sigma, norm, beta)) << endl; 
    fprintf(Lambda_Shape_out, "%s\t%s\t%s\t%s\t%s\n", conv(lambda(ilambda)).c_str(), conv(sqrt(W_func_comp(0,Corr_Z, Cov, gl.row(ilambda), 0, A0, f0, Estar, sigma, norm, beta)/Target_Int(0, Estar, sigma, norm, beta))).c_str(), conv(sqrt(W_func_comp(0,Corr_Z, Cov, gl.row(ilambda), 0, A0, f0, Estar, sigma, norm, beta)/W_func_comp(1, Corr_Z, Cov, gl.row(ilambda), 0, A0, f0, Estar, sigma, norm, beta))).c_str(), conv(Dens_Mu(ilambda)).c_str(), conv(Dens_S(ilambda)).c_str());
  }

  fclose(Lambda_Shape_out);

  FILE *file_final_points;
  strcat(output_file, "_final_res");
  if ((file_final_points = fopen(output_file, "w")) == NULL)
  {
    printf("hey\n");
    printf("Error opening the input file: %s\n", output_file);
    exit(EXIT_FAILURE);
  }

  fprintf(file_final_points, "#lambda\tA[g]/Tar_Int\td[g_t]\tdens_mu(lambda)\tdens_err(lambda)\n");


  fprintf(file_final_points, "%s\t%s\t%s\t%s\t%s\n", conv(lambda(ilambda1)).c_str(), conv(sqrt(W_func_comp(0,Corr_Z, Cov, gl.row(ilambda1), 0, A0, f0, Estar, sigma, norm, beta)/Target_Int(0, Estar, sigma, norm, beta))).c_str(), conv(sqrt(W_func_comp(0,Corr_Z, Cov, gl.row(ilambda1), 0, A0, f0, Estar, sigma, norm, beta)/W_func_comp(1, Corr_Z, Cov, gl.row(ilambda1), 0, A0, f0, Estar, sigma, norm, beta))).c_str(), conv(Dens_Mu(ilambda1)).c_str(), conv(Dens_S(ilambda1)).c_str());

  fprintf(file_final_points, "%s\t%s\t%s\t%s\t%s\n", conv(lambda(ilambda2)).c_str(), conv(sqrt(W_func_comp(0,Corr_Z, Cov, gl.row(ilambda2), 0, A0, f0, Estar, sigma, norm, beta)/Target_Int(0, Estar, sigma, norm, beta))).c_str(), conv(sqrt(W_func_comp(0,Corr_Z, Cov, gl.row(ilambda2), 0, A0, f0, Estar, sigma, norm, beta)/W_func_comp(1, Corr_Z, Cov, gl.row(ilambda2), 0, A0, f0, Estar, sigma, norm, beta))).c_str(), conv(Dens_Mu(ilambda2)).c_str(), conv(Dens_S(ilambda2)).c_str());

  fprintf(file_final_points, "%s\t%s\t%s\t%s\t%s\n", conv(lambda(ilambda2)).c_str(), conv(sqrt(W_func_comp(0,Corr_Z, Cov, gl.row(ilambda2), 0, A0, f0, Estar, sigma, norm, beta)/Target_Int(0, Estar, sigma, norm, beta))).c_str(), conv(sqrt(W_func_comp(0,Corr_Z, Cov, gl.row(ilambda2), 0, A0, f0, Estar, sigma, norm, beta)/W_func_comp(1, Corr_Z, Cov, gl.row(ilambda2), 0, A0, f0, Estar, sigma, norm, beta))).c_str(), conv(Dens_Mu(ilambda2)).c_str(), conv(SigmaF).c_str());

  fclose(file_final_points);  
}
