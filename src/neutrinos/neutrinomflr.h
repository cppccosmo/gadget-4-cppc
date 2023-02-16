#ifndef NEUTRINOMFLR_H
#define NEUTRINOMFLR_H

#include "../data/dtypes.h"
#include "../main/main.h"
#include "gadgetconfig.h"

struct nulinear
{
public:
    const double N_nu_eff = 3.044;
    const double N_nu_massive = 3.044;
    static const int N_tau = 20; // sets the number of neutrino fluids 
    static const int N_mu = 20; // sets the highest legendre moment

#define N_nu_tot (2*N_tau*N_mu)
#define N_EQ (2*N_tau*N_mu+2)

    // homogeneous evolution functions
    double tau_t_eV(int t);
    double v_t_eta(int t, double eta);
    double v2_t_eta(int t, double eta);
    double v2_t_eta_REL(int t, double eta);
    double Ec_t_eta(int t, double eta);
    double dlnEc_t_eta(int t, double eta);
    double v_t_eta_REL(int t, double eta);
    double Ec_t_eta_REL(int t, double eta);
    double dlnEc_t_eta_REL(int t, double eta);
    double eta_convert(double a);
    double Ec_de_eta(double eta);
    double dlnEc_de_eta(double eta);
    double Hc2_Hc02_eta(double eta);
    double Hc_eta(double eta);
    static double Hubble_Extern(double a);
    double dlnHc_eta(double eta);
    double OF_eta(int F, double eta);
    double Poisson(double eta, double k, const double *y);
    double m_nu_eV_parser(void);
    int N_EQ_parser(void);
    int N_nu_tot_parser(void);
    int N_tau_parser(void);
    int N_mu_parser(void);

    inline double fmin(double x, double y);
    inline double fmax(double x, double y);
    double d_nu_mono(double z, const double *y);
    void d_nu_mono_stream(double z, const double *y);
    void print_delta(const double *w);
    void print_theta(const double *w);

    double dn(int alpha, int ell, const double y[]);
    double tn(int alpha, int ell, const double y[]);
    double gt(int alpha, int ell, const double y[]);
    //int der(double eta, const double *y, double *dy, void *par);
    
    int evolve_to_z(double k, double z, double *w);
    int evolve_step(double k, double z0, double z1, double *w);
    double compute_deviation(double k, double z, double *w);

    // SuperEasy linear response functions
    double poisson_mod_fac(double k, double a);

    double y_nu[N_EQ*PMGRID];
    int initialisation_switch = 1;
    double TimeOld;
//private:
    const double Hc0h = 3.33564095198152e-04;
#define Hc0h2 (Hc0h*Hc0h)
    const double aeta_in = 1e-3;
#define eta_stop (-log(aeta_in))
    
    const double T_CMB_0_K = 2.726;
   
#define m_nu_eV (93.259*(All.OmegaNuLin+All.OmegaNuPart)*All.HubbleParam*All.HubbleParam/N_nu_massive)
#define Omega_nu_t_0 ((All.OmegaNuLin+All.OmegaNuPart)/N_tau)
   
#define T_CMB_0_K_4 (T_CMB_0_K*T_CMB_0_K*T_CMB_0_K*T_CMB_0_K)
#define T_NUREL_0_K (0.713765855503608*T_CMB_0_K)
#define m_T_nu (m_nu_eV * 11604.51812 / T_NUREL_0_K)
#define Omega_gam_0 ((4.46911743913795e-07)*T_CMB_0_K_4/(All.HubbleParam*All.HubbleParam))
#define Omega_nurel_0 (0.227107317660239*(N_nu_eff-N_nu_massive)*Omega_gam_0)
#define Omega_rel_0 (Omega_gam_0+Omega_nurel_0)
#define Omega_de_0 (1.0-All.Omega0-All.OmegaNuPart-All.OmegaNuLin-Omega_rel_0)
    
    const double w_eos_cdm = 0;
    const double w_eos_gam = 0.33333333333333333;
    const double w0_eos_de = -1.0;
    const double wa_eos_de = 0.0;

    const double PARAM_DETA0 = 1e-6;
    const double PARAM_EABS = 0;
    const double PARAM_EREL = 1e-6;
    const double C_kfs = 0.759364372216230; //sqrt(ln(2) / zeta(3)) 
   
    //int der(double eta, const double *y, double *dy, void *par);
};

extern nulinear Nulinear;

#endif
