////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    Copyright 2023 Amol Upadhye                                             //
//                                                                            //
//    This file is part of MuFLR-HDM.                                         //
//                                                                            //
//    MuFLR-HDM is free software: you can redistribute it and/or modify       //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation, either version 3 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    MuFLR-HDM is distributed in the hope that it will be useful,            //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with MuFLR-HDM.  If not, see <http://www.gnu.org/licenses/>.      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_spline.h>

///////////#include "AU_pcu.h"
#include "AU_ncint.h"
#include "AU_fftgrid.h"
#include "AU_cosmoparam.h"
#include "AU_cosmofunc.h"
#include "AU_fastpt_coord.h"

//////////////////////////// SWITCHES AND TOLERANCES ///////////////////////////

const int SWITCH_NU_SOURCE_NONLIN = 1; //source nu growth using nonlin CB

const double PARAM_DETA0 = 1e-3; //default starting step size in eta
const double PARAM_EABS = 0; //absolute error tolerance
const double PARAM_EREL = 1e-3; //relative error tolerance

////////////////////////////// PERTURBATIONS ///////////////////////////////////

//Perturbation array has dimensionality N_EQ
//
//neutrinos: 2 * N_tau * N_mu * NK 
//  y[(2*alpha + 0)*N_mu*NK + ell*NK + i] = delta_{alpha,ell}(k_i)
//  y[(2*alpha + 1)*N_mu*NK + ell*NK + i] = theta_{alpha,ell}(k_i)
//
//CDM+Baryons: 5*NK elements
//  y[2*N_tau*N_mu*NK + 0*NK + i] = delta_{CB}(k_i) (lin)
//  y[2*N_tau*N_mu*NK + 1*NK + i] = theta_{CB}(k_i) (lin)
//  y[2*N_tau*N_mu*NK + 2*NK + i] = delta_{CB}(k_i) (non-lin)
//  y[2*N_tau*N_mu*NK + 3*NK + i] = theta_{CB}(k_i) (non-lin)
//  y[2*N_tau*N_mu*NK + 4*NK + i] = r_{CB}(k_i) (non-lin)

inline double yhdm0(int alpha, int ell, int ik, const double *y){
  return y[(2*alpha+0)*N_mu*NK + ell*NK + ik];
}

inline double yhdm1(int alpha, int ell, int ik, const double *y){
  return y[(2*alpha+1)*N_mu*NK + ell*NK + ik];
}

double yhdm(int ny, int alpha, int ell, int ik, const double *y){
  if(ny<0 || alpha<0 || alpha>=N_tau || ell<0 || ell>=N_mu || ik<0 || ik>=NK)
    return 0;
  if(ny==0) return yhdm0(alpha,ell,ik,y);
  else if (ny==1) return yhdm1(alpha,ell,ik,y);
  return 1e100; //should not get here
}

inline double ycb0l(int ik, const double *y){
  return y[2*N_tau*N_mu*NK + 0*NK + ik];
}

inline double ycb1l(int ik, const double *y){
  return y[2*N_tau*N_mu*NK + 1*NK + ik];
}

inline double Pcb0n(int ik, const double *y){
  return exp( y[2*N_tau*N_mu*NK + 2*NK + ik] );
}

inline double Pcb1n(int ik, const double *y){
  return exp( y[2*N_tau*N_mu*NK + 3*NK + ik] );
}

inline double Pcb2n(int ik, const double *y){
  return exp( y[2*N_tau*N_mu*NK + 4*NK + ik] );
}

inline double Pcbn(int ab, int ik, const double *y){
  return exp( y[2*N_tau*N_mu*NK + (2+ab)*NK + ik] );
}

//hdm density monopole from perturbation array
double d_hdm_mono(int ik, double z, const double *y){
  
  double d_mono = 0, norm = 0, aeta = 1.0/(1.0+z), eta = log(aeta/aeta_in);
  
  for(int t=0; t<N_tau; t++){
    double E_m = 1.0;
    d_mono += y[2*t*N_mu*NK + 0*NK + ik] * E_m;
    norm += E_m;
  }

  return d_mono / norm;
}  

//functions for extrapolating Legendre moments, using first and second
//finite-difference approximations
double Lex(int type, int ny, int alpha, int ell, int ik, const double *y){

  //if no extrapolation needed, return actual element from y
  if(ell>=0 && ell < N_mu) return yhdm(ny,alpha,ell,ik,y);

  if(type==1 && ell<=N_mu){
    int lm = ell - 1;
    double fac0 = 3.0 * lm * (2.0*lm+3.0)
      / ( (1.0+lm)*(2.0*lm+1.0) );
    double fac1 = -3.0 * (-1.0+lm) * (2.0*lm+3.0)
      / ( (1.0+lm)*(2.0*lm-1.0) );
    double fac2 = 1.0 * (-2.0+lm) * (2.0*lm+3.0)
      / ( (1.0+lm) * (2.0*lm-3.0) );
    return fac0 * yhdm(ny,alpha,lm,ik,y)
      + fac1 * yhdm(ny,alpha,lm-1,ik,y)
      + fac2 * yhdm(ny,alpha,lm-2,ik,y);
  }

  //shouldn't get here
  printf("ERROR: Invalid moment or extrapolation type in Lex.  Quitting.\n");
  fflush(stdout);
  abort();
  return 0;
}

//Gaussian or Lorentzian filter for power spectrum input to mode-coupling integ
inline double filter(double k, double sigv2){
  return 1.0 / (1.0 + k*k*sigv2);
}

//function for extrapolating P_00, P01, P11 for padded k grid
double pad_power(double eta, const double *Pin, double *Ppad,
		 const struct cosmoparam C){

  //smoothing scale
  double kP[NK];
  double pre = 0.01 * exp(2.0*eta) / (6.0*M_PI*M_PI);
  for(int ik=0; ik<NK; ik++) kP[ik] = exp(LNKMIN+DLNK*ik)*Pin[0*NK+ik];
  double sigv2 = pre * ncint_cf(NK, DLNK, kP);
  
  //extrapolate left using P \propto k^ns
  for(int i=0; i<NSHIFT; i++){
    double lnk = LNK_PAD_MIN+DLNK*i, k=exp(lnk), kr=pow(k/KMIN,C.n_s);
    double Fk = filter(k,sigv2);
    Ppad[i]       = Pin[0] * kr * Fk;
    Ppad[i+NKP]   = Pin[NK] * kr * Fk;
    Ppad[i+2*NKP] = Pin[2*NK] * kr * Fk;
  }

  //no extrapolation needed
  for(int i=NSHIFT; i<NSHIFT+NK; i++){
    double lnk = LNK_PAD_MIN+DLNK*i, k=exp(lnk), Fk = filter(k,sigv2);
    Ppad[i]       = Pin[i-NSHIFT] * Fk;
    Ppad[i+NKP]   = Pin[NK+i-NSHIFT] * Fk;
    Ppad[i+2*NKP] = Pin[2*NK+i-NSHIFT] * Fk;
  }

  //extrapolate right using P \propto k^ns T^2
  double T1 = T_EH(KMAX,C);
  for(int i=NSHIFT+NK; i<NKP; i++){
    double lnk = LNK_PAD_MIN + DLNK*i, k=exp(lnk), kr=pow(k/KMAX,C.n_s);
    double T2 = sq(T_EH(k,C)/T1), kns_T2 = kr * T2, Fk = filter(k,sigv2);
    Ppad[i]   	  = Pin[NK-1] * kns_T2 * Fk;
    Ppad[i+NKP]   = Pin[2*NK-1]	* kns_T2 * Fk;
    Ppad[i+2*NKP] = Pin[3*NK-1]	* kns_T2 * Fk;
  }
  
  return sqrt(sigv2);
}

//Time-RG indices
const int nUI = 14; //number of unique indices of A, I
const int aU[] = {0,0,0,0,0,0,0,0, 1,1,1,1,1,1};
const int cU[] = {0,0,0,0,0,0,0,0, 1,1,1,1,1,1};
const int dU[] = {1,1,1,1,1,1,1,1, 1,1,1,1,1,1};
const int bU[] = {0,0,0,0,1,1,1,1, 0,0,0,1,1,1};
const int eU[] = {0,0,1,1,0,0,1,1, 0,0,1,0,0,1};
const int fU[] = {0,1,0,1,0,1,0,1, 0,1,1,0,1,1};
const int JU[] = {8,9,10,11,12,13,14,15,56,57,59,60,61,63};

inline int nAI(int a, int c, int d, int b, int e, int f){
  return 32*a + 16*c + 8*d + 4*b + 2*e + f; }

double Itrg(int a, int c, int d, int b, int e, int f, int i, const double *y){
  if(a==0 && c==1 && d==0) return Itrg(a,d,c,b,f,e,i,y);
  int nI = nAI(a,c,d,b,e,f);
  if(nI<8 || nI>63) return 0;
  else if(nI<16) return y[2*N_tau*N_mu*NK + (5+nI-8)*NK + i];
  else if(nI<56) return 0;
  return y[2*N_tau*N_mu*NK + (5+8+3*b+e+f)*NK + i];
}

////////////////////////////// PRINT RESULTS ///////////////////////////////////

//print growth factor D, growth rate f, and total nu growth 
int print_all_growth(double z, const double *w){
  for(int i=0; i<NK; i++){
    double k = KMIN * exp(DLNK * i);
    printf("%g %g %g %g %g\n", z, k,
           w[2*N_tau*N_mu*NK + 0*NK + i],
           w[2*N_tau*N_mu*NK + 1*NK + i] / w[2*N_tau*N_mu*NK + 0*NK + i],
           d_hdm_mono(i,z,w));
    fflush(stdout);
  }
  return 0;
}

//print linear and nonlinear cb power and total nu power
int print_all_Pcblin_Pcbnl_Pnutot(double z, const double *w){
  double a2_ain2 = 1.0 / (sq( (1.0+z) * aeta_in ));
  for(int i=0; i<NK; i++){
    double k = KMIN*exp(DLNK*i), dl = ycb0l(i,w), tl = ycb1l(i,w);
    printf("%e   %e %e %e   %e %e %e   %e\n", k, sq(dl), dl*tl, sq(tl),
	   Pcb0n(i,w)*a2_ain2, Pcb1n(i,w)*a2_ain2, Pcb2n(i,w)*a2_ain2,
	   sq(d_hdm_mono(i,z,w)) );
    fflush(stdout);
  }
  return 0;
}

//print dd, dt, tt monopole powers for cb(lin), cb(nl) and all nu fluids
int print_all_Pmono(double z, const double *w){
  double a2_ain2 = 1.0 / (sq( (1.0+z) * aeta_in ));
  for(int i=0; i<NK; i++){
    printf("%e   %e %e %e   %e %e %e",
           KMIN*exp(DLNK*i),
           sq(ycb0l(i,w)), ycb0l(i,w)*ycb1l(i,w), sq(ycb1l(i,w)),
	   Pcb0n(i,w)*a2_ain2, Pcb1n(i,w)*a2_ain2, Pcb2n(i,w)*a2_ain2);
    for(int t=0; t<N_tau; t++)
      printf("   %e %e %e", sq(yhdm0(t,0,i,w)),
	     yhdm0(t,0,i,w)*yhdm1(t,0,i,w), sq(yhdm1(t,0,i,w)));
    printf("\n");
    fflush(stdout);
  }
  return 0;
}

//print all nu and cb perturbations
int print_all_perturbations(double z, const double *w){
  for(int i=0; i<NK; i++){
    printf("%e", KMIN*exp(DLNK*i));
    for(int j=0; j<2*N_tau*N_mu+19; j++) printf(" %e", w[j*NK+i]);
    printf("\n");
    fflush(stdout);
  }
  return 0;
}  

//print user-determined information
int print_menu(int ptype, double z, const double *w){
  switch(ptype){
  case 0: return print_all_growth(z,w); 
  case 1: return print_all_Pcblin_Pcbnl_Pnutot(z,w);
  case 2: return print_all_Pmono(z,w);
  case 3: return print_all_perturbations(z,w);
  default: return 1;
  }
  return 1;
}

/////////////////////////////// DERIVATIVES ////////////////////////////////////

int der(double eta, const double *y, double *dy, void *par){

  //initialize
  struct cosmoparam *C = (struct cosmoparam *)par;
  double Hc2_Hc02 = Hc2_Hc02_eta(eta,*C), Hc2 = Hc2_Hc02*Hc0h2, Hc = sqrt(Hc2),
    dlnHc = dlnHc_eta(eta,*C), ee = exp(eta), aeta = aeta_in*ee, ze=1.0/aeta-1;
  double Xi[NK][2][2], fhc=C->f_hdm_0/C->f_cb_0; //linear evolution matrix
  for(int ieq=0; ieq<N_EQ; ieq++) dy[ieq] = 0;
  
  //initialize non-linear calculation
  double Ppad[3*NKP], Aacdbef[64*NK];
  if(C->switch_nonlinear && ze<C->z_nonlinear_initial){
    double Pin[3*NK];

#pragma omp parallel for schedule(dynamic)
    for(int i=0; i<NK; i++){
      Pin[0*NK+i] = Pcbn(0,i,y);
      Pin[1*NK+i] = Pcbn(1,i,y);
      Pin[2*NK+i] = Pcbn(2,i,y);
    }
      
    double sigv = pad_power(eta, Pin, Ppad, *C);
    compute_Aacdbef(eta,Ppad,Aacdbef);
  }
  
  //loop over wave numbers
#pragma omp parallel for schedule(dynamic)
  for(int i=0; i<NK; i++){
    
    double k = KMIN * exp(DLNK * i), k_H = k/Hc, k2_H2=k_H*k_H;
    double Phi_l = Poisson_lin(eta,i,y,*C), Phi = Phi_l;
    if(SWITCH_NU_SOURCE_NONLIN) Phi = Poisson_nonlin(eta,i,y,*C);
    
    //linear evolution matrix for CDM+Baryon Time-RG
    double dhdm = d_hdm_mono(i,ze,y), dcb = ee*sqrt(Pcb0n(i,y));
    Xi[i][0][0] = 1.0;
    Xi[i][0][1] = -1.0;
    Xi[i][1][0] = -1.5 * OF_eta(N_tau,eta,*C) * (1.0 + fhc*dhdm/dcb);
    Xi[i][1][1] = 2.0 + dlnHc;
    
    //neutrino stream perturbations
    for(int t=0; t<N_tau; t++){
      
      double vt = v_t_eta(t,eta,*C), kv_H = vt*k_H;
    
      //sum over Legendre moments of fluid equations
      for(int ell=0; ell<N_mu; ell++){
	dy[(2*t+0)*N_mu*NK + ell*NK + i]
	  = kv_H * ( yhdm(0,t,ell-1,i,y) * ell / (2*ell-1)
		     - Lex(1,0,t,ell+1,i,y) * (ell+1) / (2*ell+3) )
	  + yhdm(1,t,ell,i,y);
	
	dy[(2*t+1)*N_mu*NK + ell*NK + i]
	  = -(1.0 + dlnHc) * yhdm(1,t,ell,i,y)
	  - k2_H2 * (ell==0) * Phi
	  + kv_H * ( yhdm(1,t,ell-1,i,y) * ell / (2*ell-1)
		     - Lex(1,1,t,ell+1,i,y) * (ell+1) / (2*ell+3) );
      } //end for ell
      
    } //end for t

    //cdm perturbations: linear; always use Phi with linear delta
    dy[2*N_tau*N_mu*NK + 0*NK + i] = ycb1l(i,y);
    dy[2*N_tau*N_mu*NK + 1*NK + i] = -(1.0 + dlnHc)*ycb1l(i,y) - k2_H2*Phi_l;

    //linear evolution of power spectra
    double dP_i[3]={0,0,0}, P_i[3]={Pcbn(0,i,y),Pcbn(1,i,y),Pcbn(2,i,y)};
    dy[2*N_tau*N_mu*NK + 2*NK + i] = 0;
    dy[2*N_tau*N_mu*NK + 3*NK + i] = 0;
    dy[2*N_tau*N_mu*NK + 4*NK + i] = 0;

    for(int c=0; c<2; c++){
      dP_i[0] += -Xi[i][0][c]*P_i[c] - Xi[i][0][c]*P_i[c];
      dP_i[1] += -Xi[i][0][c]*P_i[c+1] - Xi[i][1][c]*P_i[c];
      dP_i[2] += -Xi[i][1][c]*P_i[c+1] - Xi[i][1][c]*P_i[c+1];
    }

    //non-linear evolution of power spectra
    if(C->switch_nonlinear && ze<C->z_nonlinear_initial){

      double pre = 4.0*M_PI*ee/k;
							  
      //non-linear evolution terms for power spectrum
      for(int c=0; c<2; c++){
	for(int d=0; d<2; d++){
	  dP_i[0] +=
	    pre * ( Itrg(0,c,d,0,c,d,i,y) + Itrg(0,c,d,0,c,d,i,y) );
	  dP_i[1] +=
	    pre * ( Itrg(1,c,d,0,c,d,i,y) + Itrg(0,c,d,1,c,d,i,y) );
	  dP_i[2] +=
	    pre * ( Itrg(1,c,d,1,c,d,i,y) + Itrg(1,c,d,1,c,d,i,y) );
	}//end for d
	
      } //end for c

      //instability in P_11: filtering in pad_power makes this unnecessary
      //if(dP_i[2]/P_i[2] >  10.0) dP_i[2] =  10.0*P_i[2];
      //if(dP_i[2]/P_i[2] < -10.0) dP_i[2] = -10.0*P_i[2];
      
      //non-linear and linear evolution terms for I_{acd,bef}
      for(int j=0; j<nUI; j++){
	dy[2*N_tau*N_mu*NK+(5+j)*NK+i] = 2.0*ee * Aacdbef[JU[j]*NK+i];
	
	int a=aU[j], c=cU[j], d=dU[j], b=bU[j], e=eU[j], f=fU[j]; 
	
	for(int g=0; g<2; g++)
	  dy[2*N_tau*N_mu*NK+(5+j)*NK+i] +=
	    -Xi[i][b][g]*Itrg(a,c,d,g,e,f,i,y)
	    - Xi[i][e][g]*Itrg(a,c,d,b,g,f,i,y)
	    - Xi[i][f][g]*Itrg(a,c,d,b,e,g,i,y);
      }

    } //end if nonlinear

    //power spectrum derivatives
    dy[2*N_tau*N_mu*NK + 2*NK + i] = dP_i[0] / (P_i[0] + 1e-300);
    dy[2*N_tau*N_mu*NK + 3*NK + i] = dP_i[1] / (P_i[1] + 1e-300);
    dy[2*N_tau*N_mu*NK + 4*NK + i] = dP_i[2] / (P_i[2] + 1e-300);
    
  }//end for i (loop over wave numbers)

  return GSL_SUCCESS;
}

///////////////////////////////// EVOLUTION ////////////////////////////////////

//evolve from aeta_in to input redshift
int evolve_to_z(double z, double *w, const double *Ncb_in,
		const struct cosmoparam C){
  
  //initialize perturbations at eta=0
  double aeta_eq = C.Omega_rel_0 / C.Omega_cb_0;
  double dcb_in = aeta_in + (2.0/3.0)*aeta_eq;
  for(int ieq=0; ieq<N_EQ; ieq++) w[ieq] = 0;
  
#pragma omp parallel for schedule(dynamic)
  for(int i=0; i<NK; i++){

    double k = KMIN * exp(DLNK * i);
    
    //CDM+Baryon perturbations
    w[2*N_tau*N_mu*NK + 0*NK + i] = Ncb_in[i] * dcb_in;
    w[2*N_tau*N_mu*NK + 1*NK + i] = Ncb_in[i] * aeta_in;
    w[2*N_tau*N_mu*NK + 2*NK + i] = log( sq(Ncb_in[i] * dcb_in) );
    w[2*N_tau*N_mu*NK + 3*NK + i] = log( sq(Ncb_in[i]) * dcb_in * aeta_in );
    w[2*N_tau*N_mu*NK + 4*NK + i] = log( sq(Ncb_in[i] * aeta_in) );

    //neutrino perturbations: monopoles only
    for(int t=0; t<N_tau; t++){
      double m_t = C.m_hdm_eV/tau_t_eV(t);
      double kfs2 = 1.5 * m_t*m_t * Hc0h2 * C.Omega_m_0 * aeta_in;
      double kfs = sqrt(kfs2), kpkfs = k + kfs, kpkfs2 = kpkfs*kpkfs;
      double Ft = (1.0-C.f_hdm_0) * kfs2 / (kpkfs2 - C.f_hdm_0*kfs2);
      double dlnFt = k*kpkfs / (kpkfs2 - C.f_hdm_0*kfs2);
      w[(2*t+0)*N_mu*NK + 0*NK + i] = Ft * ycb0l(i,w);
      w[(2*t+1)*N_mu*NK + 0*NK + i] = dlnFt*yhdm0(t,0,i,w) + Ft*ycb1l(i,w);
    }

  } //end i for
  
  //initialize GSL ODE integration
  int status = GSL_SUCCESS;
  struct cosmoparam par;
  copy_cosmoparam(C,&par);
  
  double eta0 = 0, aeta1 = 1.0/(1.0+z), eta1 = log(aeta1/aeta_in);
  
  gsl_odeiv2_system sys = {der, NULL, N_EQ, &par};
  
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys,
						       gsl_odeiv2_step_rkf45,
						       PARAM_DETA0,
						       PARAM_EABS,
						       PARAM_EREL);
  gsl_odeiv2_driver_set_hmax(d, 0.1);
  
  //integrate to input redshift, printing results at regular intervals
  double deta = 0.01, eta = eta0, etai = deta;
  
  while(etai < eta1
	&& status == GSL_SUCCESS){
    etai = fmin(eta+deta, eta1);
    status = gsl_odeiv2_driver_apply(d, &eta, etai, w);
  }
  
  //clean up and quit
  gsl_odeiv2_driver_free(d);
  return 0;
}

int evolve_step(double z0, double z1, double *w, const struct cosmoparam C){
  
  //initialize GSL ODE integration
  //int par = nonlin;
  struct cosmoparam par;
  copy_cosmoparam(C,&par);
  
  double aeta0 = 1.0/(1.0+z0), eta = log(aeta0/aeta_in), aeta1 = 1.0/(1.0+z1),
    eta1 = log(aeta1/aeta_in);
  
  gsl_odeiv2_system sys = {der, NULL, N_EQ, &par};
  
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys,
						       gsl_odeiv2_step_rkf45,
						       PARAM_DETA0,
						       PARAM_EABS,
						       PARAM_EREL);
  gsl_odeiv2_driver_set_hmax(d, 0.1);

  //integrate to final redshift and print results
  int status = gsl_odeiv2_driver_apply(d, &eta, eta1, w);

  //clean up and quit
  gsl_odeiv2_driver_free(d);
  return 0;
}

//////////////////////////////////// MAIN //////////////////////////////////////
//various tests of code

int main(int argn, char *args[]){

  //initialize
  struct cosmoparam C;
  initialize_cosmoparam(&C,args[1],N_tau);
  tau_t_eV_init(0,C.file_hdm_distribution,C.T_hdm_0_K);
  Pmat0(1,C);
  double N_cb[NK], *y = (double *)malloc(N_EQ*sizeof(double));
  for(int i=0; i<NK; i++) N_cb[i] = 1;

  //redshift list
  double *zn = (double *)malloc((C.num_z_outputs+1)*sizeof(double));
  zn[0] = C.z_nonlinear_initial;
  for(int iz=1; iz<=C.num_z_outputs; iz++) zn[iz] = C.z_outputs[iz-1];
  
  //linear run: store output y and then normalize at finish
  if(C.switch_nonlinear == 0){
    double *yzn = (double *)malloc(C.num_z_outputs*N_EQ*sizeof(double));
    for(int i=0; i<C.num_z_outputs*N_EQ; i++) yzn[i] = 0;

    evolve_to_z(zn[0],y,N_cb,C);
    for(int iz=0; iz<C.num_z_outputs; iz++){
      evolve_step(zn[iz],zn[iz+1],y,C);
      for(int i=0; i<2*N_tau*N_mu*NK + 19*NK; i++) yzn[iz*N_EQ+i] = y[i];
    }

    if(zn[C.num_z_outputs] > 1e-9) evolve_step(zn[C.num_z_outputs],0,y,C);
    for(int i=0; i<NK; i++){
      double dU = C.f_cb_0*y[2*N_tau*N_mu*NK+i] + C.f_hdm_0*d_hdm_mono(i,0,y);
      N_cb[i] *= sqrt(Pmat0(KMIN*exp(DLNK*i),C)) / dU;
    }
    
    for(int iz=0; iz<C.num_z_outputs; iz++){
      for(int ik=0; ik<NK; ik++){
        for(int j=0; j<2*N_tau*N_mu+2; j++) yzn[iz*N_EQ+j*NK+ik] *= N_cb[ik];
        for(int j=2*N_tau*N_mu+2; j<2*N_tau*N_mu+5; j++)
          yzn[iz*N_EQ+j*NK+ik] += 2.0*log(N_cb[ik]);
      }

      printf("###main: output at z=%g\n",zn[iz+1]);
      print_menu(C.switch_print_linear,zn[iz+1],yzn+iz*N_EQ);
      printf("\n\n");
      fflush(stdout);
    }
      
    free(yzn);
  }  
  else{  //non-linear run: normalize power using linear evol
    int Cnl = C.switch_nonlinear;
    C.switch_nonlinear = 0;
    evolve_to_z(0,y,N_cb,C);
    for(int i=0; i<NK; i++){
      double dU = C.f_cb_0*y[2*N_tau*N_mu*NK+i] + C.f_hdm_0*d_hdm_mono(i,0,y);
      N_cb[i] *= sqrt(Pmat0(KMIN*exp(DLNK*i),C)) / dU;
    }
    
    //turn non-linearity back on and integrate forwards
    C.switch_nonlinear = Cnl;
    evolve_to_z(C.z_nonlinear_initial,y,N_cb,C);
    
    for(int iz=0; iz<C.num_z_outputs; iz++){
      evolve_step(zn[iz],zn[iz+1],y,C);
      printf("###main: output at z=%g\n",zn[iz+1]);
      print_menu(C.switch_print_linear,zn[iz+1],y);
      printf("\n\n");
      fflush(stdout);
    }
  }
  
  free(y);
  free(zn);
  tau_t_eV(FREE_TAU_TABLE);
  return 0;
}
