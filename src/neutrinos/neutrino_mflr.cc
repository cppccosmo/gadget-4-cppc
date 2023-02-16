#include "gadgetconfig.h"

#include <fftw3.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <algorithm>
#include <cassert>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>

#include "pcu.h"
#include "neutrinomflr.h"
#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../domain/domain.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../pm/pm.h"
#include "../pm/pm_periodic.h"
#include "../sort/cxxsort.h"
#include "../system/system.h"
#include "../time_integration/timestep.h"

double nulinear::tau_t_eV(int t) {
    if (N_tau==0) return 0.0;
    static int init = 0;
    static double *tau_table_eV;
    
    if(!init) {
        tau_table_eV = (double *)Mem.mymalloc_movable(&tau_table_eV, "tau_table_eV", N_tau * sizeof(double));
        gsl_interp_accel *spline_accel = gsl_interp_accel_alloc();
        gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen,pcu_N);
        gsl_spline_init(spline,pcu_prob,pcu_tau,pcu_N);
        
        for(int t=0; t<N_tau; t++){
            double prob = (0.5+t) / N_tau;
            tau_table_eV[t] = gsl_spline_eval(spline,prob,spline_accel);
        }
        
        gsl_spline_free(spline);
        gsl_interp_accel_free(spline_accel);
        init = 1;
    }
    
    if(t == -1){
        Mem.myfree(tau_table_eV);
        init = 0;
        return 0;
    }
    return tau_table_eV[t];
}

// speed -tau_t / tau0_t of each neutrino species
double nulinear::v_t_eta(int t, double eta) {
    double t_ma = tau_t_eV(t) / ( m_nu_eV * aeta_in * exp(eta) );
    return (t_ma<1 ? t_ma : 1);
}

double nulinear::v2_t_eta(int t, double eta) {
    double vt = Nulinear.v_t_eta(t,eta);
    return vt*vt;
}

double nulinear::v2_t_eta_REL(int t, double eta) {
    double m_aeta_tau = m_nu_eV * aeta_in*exp(eta) / Nulinear.tau_t_eV(t);
    return 1.0 / (1.0 + m_aeta_tau*m_aeta_tau);
}

double nulinear::v_t_eta_REL(int t, double eta){ return sqrt(Nulinear.v2_t_eta_REL(t,eta)); }

// density ratio rho_t(eta)/rho_t(eta_stop) * aeta^2 and its log deriv 
double nulinear::Ec_t_eta(int t, double eta) { return 1.0 / (aeta_in*exp(eta) ); }

double nulinear::dlnEc_t_eta(int t, double eta){ return -1.0; }

double nulinear::Ec_t_eta_REL(int t, double eta) {
    double vt2 = Nulinear.v2_t_eta_REL(t,eta), aeta = aeta_in*exp(eta);
    if(1-vt2 < 1e-12){
        double ma_tau = m_nu_eV * aeta / Nulinear.tau_t_eV(t);
        return sqrt(1.0 + ma_tau*ma_tau) / (aeta*ma_tau);
    }
    return 1.0 / (aeta * sqrt(1.0 - vt2));
}

double nulinear::dlnEc_t_eta_REL(int t, double eta){ return -1.0 - v2_t_eta_REL(t,eta); }

double nulinear::eta_convert(double a) {
    return log(a/aeta_in);
}

/////// Homogeneous cosmology ////////

double nulinear::Ec_de_eta(double eta) {
    double aeta = aeta_in * exp(eta);
    return pow(aeta, -1.0 -3.0*(w0_eos_de + wa_eos_de)) * exp(3.0*wa_eos_de*(aeta-1.0));
}

double nulinear::dlnEc_de_eta(double eta) {
    double aeta = aeta_in * exp(eta);
    return -1.0 - 3.0*(w0_eos_de + wa_eos_de) + 3.0*wa_eos_de*aeta;
}

// conformal Hubble parameter
double nulinear::Hc2_Hc02_eta(double eta) {
    //scale factor
    double aeta = aeta_in * exp(eta), aeta2 = aeta*aeta, Ec_de = Nulinear.Ec_de_eta(eta);
    
    //sum Omega_{t,0} aeta^2 rho_t(eta)/rho_t_0 over CDM, photons, and DE
    double sum_OEc = All.Omega0/aeta + Omega_rel_0/aeta2 + All.OmegaLambda*Ec_de;
    
    //neutrinos
    if(All.OmegaNuLin != 0 || All.OmegaNuPart != 0) {
      for(int t=0; t<N_tau; t++) sum_OEc += Omega_nu_t_0 * Nulinear.Ec_t_eta_REL(t,eta);
    }

    return sum_OEc;
}

double nulinear::Hubble_Extern (double a) {
    return sqrt(Nulinear.Hc2_Hc02_eta(Nulinear.eta_convert(a)))/a;
}

// density fraction in spatially-flat universe
double nulinear::OF_eta(int F, double eta) {
    double Hc02_Hc2 = 1.0 / Nulinear.Hc2_Hc02_eta(eta), aeta = aeta_in*exp(eta);

    if(F == N_tau) // CDM + particle neutrinos
      return (All.Omega0 + All.OmegaNuPart) * pow(aeta, -1.0-3.0*w_eos_cdm) * Hc02_Hc2;
    else if(F == N_tau+1) // photon + massless nu
      return Omega_rel_0 * pow(aeta, -1.0-3.0*w_eos_gam) * Hc02_Hc2;
    else if(F == N_tau+2) // dark energy, assumed Lambda
      return Omega_de_0 * Nulinear.Ec_de_eta(eta) * Hc02_Hc2;
    else if(F<0 || F>N_tau+2) return 0.0; // no fluids  should have these indices
    return Omega_nu_t_0 * Nulinear.Ec_t_eta(F,eta) * Hc02_Hc2;
}

double nulinear::Hc_eta(double eta) { return Hc0h * sqrt(Nulinear.Hc2_Hc02_eta(eta)); }

// d log(Hc) / d eta
double nulinear::dlnHc_eta(double eta) {
    double aeta = aeta_in*exp(eta), aeta2 = aeta*aeta;
    double pre = 1.0 / (2.0 * Nulinear.Hc2_Hc02_eta(eta) );

    double sum_OdEc = -(1.0 + 3.0*w_eos_cdm) * All.Omega0/aeta //CDM
      - (1.0 + 3.0*w_eos_gam) * Omega_rel_0/aeta2 // photons + massless nu
      + Nulinear.dlnEc_de_eta(eta) * Omega_de_0 * Nulinear.Ec_de_eta(eta); // DE

    for(int t=0; t<N_tau; t++) // neutrino fluids
      sum_OdEc += Nulinear.dlnEc_t_eta_REL(t,eta) * Nulinear.Ec_t_eta_REL(t,eta) * Omega_nu_t_0;

    return pre * sum_OdEc;
}

// Poisson equation for Phi
double nulinear::Poisson(double eta, double k, const double *y) {
    double Hc2 = Hc0h2 * Nulinear.Hc2_Hc02_eta(eta), pre = -1.5 * Hc2 / (k*k);
    double sum_Od = OF_eta(N_tau,eta)*y[2*N_tau*N_mu];

// this is the default setting, where the neutrino flows are converted from the slowest to fastest. The flows are therefore excluded from tau = 0.
// // if however one wishes to exclude specific neutrino flows from the middle of the distribution (for convergence tests for example), then comment 
// out the follow section of the code until the comment-line 'end-of-block'.

    int t_min = 0;
#ifdef ADDITIONAL_GRID
    t_min = All.N_tau_part * All.Nu_part_deg;
#endif
    for(int t=t_min; t<N_tau; t++) sum_Od += Nulinear.OF_eta(t,eta) * y[2*t*N_mu];

// end-of-block

// in the case of the above block of code being commented out, i.e. you wish to exclude specific neutrino flows, uncomment the following block of code. 
// // The two integers t_ex_start and t_ex_finish are the neutrino flow numbers you wish to exclude (i.e. converted to particles). The flows excluded are inclusive of the start and finish numbers. 
// // Note: the flow numbers are the 'real' number labels that start from tau=1 to tau_max, and not the c-index that is offsetted (i.e. start from tau=0).
/*
    int t_ex_start = 3;
    int t_ex_finish = 4;

    for(int t=0; t<t_ex_start-1; t++) {
      sum_Od += Nulinear.OF_eta(t,eta) * y[2*t*N_mu];
    }

    for(int t=t_ex_finish; t<N_tau; t++) {
      sum_Od += Nulinear.OF_eta(t,eta) * y[2*t*N_mu];
    }
*/

    return pre * sum_Od;
}

/////// Utility functions ///////

// minimum, maximum functions
inline double nulinear::fmin(double x, double y){ return (x<y ? x : y); }
inline double nulinear::fmax(double x, double y){ return (x>y ? x : y); }

// neutrino density monopole from perturbation array
double nulinear::d_nu_mono(double z, const double *y) {
    double d_mono = 0, norm = 0, aeta = 1.0/(1.0+z), eta = log(aeta/aeta_in);

// this is the default setting, where the neutrino flows are converted from the slowest to fastest. The flows are therefore excluded from tau = 0.
// if however one wishes to exclude specific neutrino flows from the middle of the distribution (for convergence tests for example), then comment out the follow section of the code until the comment-line 'end-of-block'.

    int t_min = 0;
#ifdef ADDITIONAL_GRID
    t_min = All.N_tau_part * All.Nu_part_deg;
#endif
    for(int t=t_min; t<N_tau; t++){
      double E_m = 1.0;
      d_mono += y[2*t*N_mu] * E_m;
      norm += E_m;
    }

// end-of-block

// in the case of the above block of code being commented out, i.e. you wish to exclude specific neutrino flows, uncomment the following block of code. 
// The two integers t_ex_start and t_ex_finish are the neutrino flow numbers you wish to exclude (i.e. converted to particles). The flows excluded are inclusive of the start and finish numbers. 
// Note: the flow numbers are the 'real' number labels that start from tau=1 to tau_max, and not the c-index that is offsetted (i.e. start from tau=0).
/*
    int t_ex_start = 3;
    int t_ex_finish = 3;

    for(int t=0; t<t_ex_start-1; t++) {
      double E_m = 1.0;
      d_mono += y[2*t*N_mu] * E_m;
      norm += E_m;
    }

    for(int t=t_ex_finish; t<N_tau; t++) {
      double E_m = 1.0;
      d_mono += y[2*t*N_mu] * E_m;
      norm += E_m;
    }
*/

    if(norm > 0) {
      return d_mono / norm;
    } else {
      return 0;
    }
}

// neutrino monopole of individual streams
void nulinear::d_nu_mono_stream(double z, const double *y) {
    double norm = 0, aeta = 1.0/(1.0+z), eta = log(aeta/aeta_in);

    for(int t=0; t<N_tau; t++) {
      double d_mono = y[2*t*N_mu];
      printf("%f,", d_mono);
    }
}

void nulinear::print_delta(const double *w) {
/*
    for(int i=0; i<N_mu; i++) {
      mpi_printf("%f,", w[2*i]);
    }
*/
    printf("%g,", w[0]);
}

void nulinear::print_theta(const double *w) { 
/*
    for(int i=0; i<N_mu; i++) {
      mpi_printf("%f,", w[2*i+1]);
    }
*/
    printf("%g,", w[1]);
}

/////// Derivatives ///////

//neutrino perturbation variables: X_{F,ell} with X of delta or theta,
////F referring to fluid (0 to N_tau-1 for nu; N_tau for CDM; etc.)
////and ell referring to Legendre coefficient (0 to N_mu-1)
////  y[0] = delta_{0,0} 
////  y[1] = theta_{0,0}
////  y[2] = delta_{0,1}
////  y[3] = theta_{0,1}
////      ...
////  y[2*N_mu+0] = delta_{1,0}
////  y[2*N_mu+1] = theta_{1,0}
////             ...
////  y[2*N_tau*N_mu - 2] = delta_{N_tau-1,N_mu-1}
////  y[2*N_tau*N_mu - 1] = theta_{N_tau-1,N_mu-1}
////
////cdm perturbation variables come after the neutrinos
////  y[2*N_tau*N_mu + 0] = delta_{CDM}
////  y[2*N_tau*N_mu + 1] = theta_{CDM}

double nulinear::dn(int alpha, int ell, const double y[]) {
    return (ell<0 ? 0 : y[2*alpha*N_mu + 2*ell]);
}

double nulinear::tn(int alpha, int ell, const double y[]) {
    return (ell<0 ? 0 : y[2*alpha*N_mu + 2*ell + 1]); 
}

double nulinear::gt(int alpha, int ell, const double y[]) {
    return Nulinear.tn(alpha,ell,y)*(double)(ell*(ell-1))/(double)((2*ell-1)*(2*ell+1)); 
}

int der(double eta, const double *y, double *dy, void *par) {
    // initialise
    double *pd = (double *)par, k = pd[0], k_H = k/Nulinear.Hc_eta(eta), k2_H2 = k_H*k_H,
      Phi = Nulinear.Poisson(eta,k,y), aeta = Nulinear.aeta_in*exp(eta), dlnHc = Nulinear.dlnHc_eta(eta);

    // neutrino stream perturbations 
    for(int t=0; t<Nulinear.N_tau; t++) {
      double vt = Nulinear.v_t_eta(t,eta), kv_H = vt*k_H;

      // sum over Legendre moments of fluid equations, except the last two which require approximation
      for(int ell=0; ell<Nulinear.N_mu-2; ell++) {
        dy[2*t*Nulinear.N_mu + 2*ell] 
          = kv_H * ( Nulinear.dn(t,ell-1,y) * ell / (2*ell-1)
                      - Nulinear.dn(t,ell+1,y) * (ell+1) / (2*ell+3) )
          + Nulinear.tn(t,ell,y);

        dy[2*t*Nulinear.N_mu + 2*ell + 1]
          = -(1.0 + dlnHc) * Nulinear.tn(t,ell,y)
          - k2_H2 * (ell==0) * Phi
          + kv_H * ( Nulinear.tn(t,ell-1,y) * ell / (2*ell-1)
                      - Nulinear.tn(t,ell+1,y) * (ell+1) / (2*ell+3) );
      }

      // special case for ell = ell_max-1 = N_mu-2
      int ell = Nulinear.N_mu-2;
      dy[2*t*Nulinear.N_mu + 2*ell]
        = kv_H * ( Nulinear.dn(t,ell-1,y) * (ell) / (2*ell-1)
                      - Nulinear.dn(t,ell+1,y) * (ell+1) / (2*ell+3) )
        + Nulinear.tn(t,ell,y);

      dy[2*t*Nulinear.N_mu + 2*ell +1]
        = -(1.0 + dlnHc) * Nulinear.tn(t,ell,y)
        + kv_H * ( Nulinear.tn(t,ell-1,y) * (ell) / (2*ell-1)
                      - Nulinear.tn(t,ell+1,y) * (ell+1) / (2*ell+3) );

      // special case for ell = ell_max = N_mu-1
      ell = Nulinear.N_mu-1;
      dy[2*t*Nulinear.N_mu + 2*ell]
        = kv_H * ( Nulinear.dn(t,ell-1,y) / (2*ell-1)
                      - 3.0 * (ell) * Nulinear.dn(t,ell,y) / (2*ell+1)
                      + 4.0 * (ell-1) * Nulinear.dn(t,ell-1,y) / (2*ell-1)
                      - (ell-2) * Nulinear.dn(t,ell-2,y) / (2*ell-3) )
        + Nulinear.tn(t,ell,y);

      dy[2*t*Nulinear.N_mu + 2*ell + 1]
        = -(1.0 + dlnHc) * Nulinear.tn(t,ell,y)
        + kv_H * ( Nulinear.tn(t,ell-1,y) / (2*ell-1)
                      - 3.0 * (ell) * Nulinear.tn(t,ell,y) / (2*ell+1)
                      + 4.0 * (ell-1) * Nulinear.tn(t,ell-1,y) / (2*ell+1)
                      - (ell-2) * Nulinear.tn(t,ell-2,y) / (2*ell-3) );
    }

    // CDM perturbations 
    dy[2*Nulinear.N_tau*Nulinear.N_mu + 0] = y[2*Nulinear.N_tau*Nulinear.N_mu + 1];
    dy[2*Nulinear.N_tau*Nulinear.N_mu + 1] = -(1.0 + dlnHc) * y[2*Nulinear.N_tau*Nulinear.N_mu + 1] - k2_H2*Phi;

    return GSL_SUCCESS;
}

/////// Evolution ///////

// evolve from aeta_in to input redshift 
int nulinear::evolve_to_z(double k, double z, double *w) {
    // initialise perturbations at eta=0
    double c_in=1, Omega_m_0 = All.Omega0+All.OmegaNuLin+All.OmegaNuPart, fnu0 = All.OmegaNuLin/Omega_m_0;
    double aeta_eq = Omega_rel_0 / All.Omega0;
    for(int F=0; F<N_EQ; F++) w[F] = 0;

    // CDM + Baryon perturbations 
    w[2*N_tau*N_mu + 0] = c_in * (aeta_in + (2.0/3.0)*aeta_eq);
    w[2*N_tau*N_mu + 1] = c_in * aeta_in;

    // neutrino perturbations: monopoles only
    for(int t=0; t<N_tau; t++) {
      double m_t = m_nu_eV/Nulinear.tau_t_eV(t);
      double kfs2 = 1.5*m_t*m_t * Hc0h2 * Omega_m_0 * aeta_in;
      double kfs = sqrt(kfs2), kpkfs = k + kfs, kpkfs2 = kpkfs * kpkfs;
      double Ft = (1.0-fnu0) * kfs2 / (kpkfs2 - fnu0*kfs2);
      double dlnFt = k*kpkfs / (kpkfs2 - fnu0*kfs2);
      w[2*t*N_mu + 0] = Ft * w[2*N_tau*N_mu + 0];
      w[2*t*N_mu + 1] = dlnFt*w[2*t*N_mu+0] + Ft*w[2*N_tau*N_mu+1];
    }

    // initialise GSL ODE integration
    int status = GSL_SUCCESS;
    double eta0  = 0, aeta1 = 1.0/(1.0+z), eta1 = log(aeta1/aeta_in), par = k;

    gsl_odeiv2_system sys = {der, NULL, N_EQ, &par};

    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, PARAM_DETA0, PARAM_EABS, PARAM_EREL);

    gsl_odeiv2_driver_set_hmax(d, 0.1);

    // integrate to input redshift
    double deta = 0.01, eta = eta0, etai = deta;

    while(etai < eta1 && status == GSL_SUCCESS) {
      etai = fmin(eta+deta, eta1);
      status = gsl_odeiv2_driver_apply(d, &eta, etai, w);
    }

    // clean up and quit
    gsl_odeiv2_driver_free(d);
    return 0;
}

int nulinear::evolve_step(double k, double z0, double z1, double *w) {
    // initialise GSL ODE integration
    double aeta0 = 1.0/(1.0+z0), eta = log(aeta0/aeta_in), aeta1 = 1.0/(1.0+z1),
    eta1 = log(aeta1/aeta_in), par = k;

    gsl_odeiv2_system sys = {der, NULL, N_EQ, &par};
  
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys,
						       gsl_odeiv2_step_rkf45,
						       PARAM_DETA0,
						       PARAM_EABS,
						       PARAM_EREL);
    gsl_odeiv2_driver_set_hmax(d, 0.1);

    // integrate to final redshift and print results
    int status = gsl_odeiv2_driver_apply(d, &eta, eta1, w);

    // clean up and quit
    gsl_odeiv2_driver_free(d);
    return 0;
}

double nulinear::poisson_mod_fac(double k, double a) {
    double fnu = All.OmegaNuLin/(All.Omega0+All.OmegaNuLin+All.OmegaNuPart);
    double kfs = 1.5 * sqrt(All.Time * (All.Omega0 + All.OmegaNuLin + All.OmegaNuPart)) * Nulinear.m_nu_eV_parser();

    return ((k + kfs) * (k + kfs)) / ((k + kfs) * (k + kfs) - kfs * kfs * fnu); 
}

double nulinear::compute_deviation(double k, double z, double *w) {
    double fnu = All.OmegaNuLin / (All.Omega0 + All.OmegaNuPart + All.OmegaNuLin);
    return 1. + fnu / (1.-fnu) * Nulinear.d_nu_mono(z,w)/w[N_nu_tot];
}

double nulinear::m_nu_eV_parser(void) {
    return m_nu_eV;
}

int nulinear::N_EQ_parser(void) {
    return N_EQ;
}

int nulinear::N_nu_tot_parser(void) {
    return N_nu_tot;
}

int nulinear::N_tau_parser(void) {
    return N_tau;
}

int nulinear::N_mu_parser(void) {
    return N_mu;
}
