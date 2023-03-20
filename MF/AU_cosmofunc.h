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

////////////////////////////////// CONSTANTS ///////////////////////////////////
//All dimensionful quantities are in units of Mpc/h to the appropriate power,
//unless otherwise noted.  Values declared as const double or const int  may be
//modified by the user, while those in #define statements are derived parameters
//which should not be changed.

//conformal hubble today
const double Hc0h = 3.33564095198152e-04; //(1e2/299792.458)
#define Hc0h2 (Hc0h*Hc0h)

//initial scale factor, and max value of eta=ln(a/a_in)
const double aeta_in = 1e-3; 
#define eta_stop (-log(aeta_in))

//gsl tolerance parameters
const double COSMOFUNC_EABS = 0; //absolute error tolerance
const double COSMOFUNC_EREL = 1e-5; //relative error tolerance

//////////////////////////////////// NEUTRINOS /////////////////////////////////

//neutrino fluid parameters
const int N_tau = 20; //number of neutrino streams; maximum 900 for this pcu.h
const int N_mu = 20; //number of multipoles to track for each stream

//total number of equations:
//  2*N_tau*N_mu*NK (delta, theta for N_tau streams * N_mu moments * NK wave#)
//  + 2*NK (delta and theta for linear CDM+Baryons)
//  + 17*NK (delta, theta, cross, I for non-linear CDM+Baryons)
#define N_EQ (2*N_tau*N_mu*NK + 19*NK)

//homogeneous-universe momentum [eV], used to identify neutrino streams
const int FREE_TAU_TABLE = -4375643; //some negative integer, pick any
const int DEBUG_NU_MOMENTA = 1;

double tau_t_eV_init(int t, const char *f0file, double T_hdm_0_K){

  if(N_tau==0) return 0.0;
  static int init = 0;
  static double *tau_table_eV;

  if(!init){

    //read input file
    FILE *fp;
    if( (fp=fopen(f0file,"r")) == NULL ){
      printf("ERROR: File %s not found.  Quitting.\n",f0file);
      exit(1);
    }
    
    int pcu_N = 0;
    //double pcu_tau[COSMOPARAM_TAU_MAX_Q], pcu_prob[COSMOPARAM_TAU_MAX_Q],
    double *pcu_tau  = malloc(COSMOPARAM_TAU_MAX_Q * sizeof(double));
    double *pcu_prob = malloc(COSMOPARAM_TAU_MAX_Q * sizeof(double));
    double q0=0, q0Last=0, f0=0, f0Last=0, T_hdm_0_eV=T_hdm_0_K/11604.525;
    pcu_prob[0] = 0;
    char buf[1000];

    while( fgets(buf, sizeof buf, fp) ){
      while(*buf=='#' || *buf=='\n'){ fgets(buf, sizeof buf, fp); }
      q0Last = q0;
      f0Last = f0;
      sscanf(buf,"%lg %lg",&q0,&f0);

      pcu_tau[pcu_N] = q0 * T_hdm_0_eV;
      if(pcu_N>0) pcu_prob[pcu_N] = pcu_prob[pcu_N-1]
		 + 2.0*M_PI*(q0-q0Last) * (q0*q0*f0 + q0Last*q0Last*f0Last);
      pcu_N++;
    }

    //normalize and truncate tail for which pcu=1
    double norm_pcu = 1.0 / pcu_prob[pcu_N-1];
    for(int i=0; i<pcu_N; i++) pcu_prob[i] *= norm_pcu;
    for(int i=pcu_N-2; i>0; i--) 
      pcu_N -= (pcu_prob[i]>=1.0) || (pcu_prob[i]>=pcu_prob[i+1]);

    tau_table_eV = malloc(N_tau * sizeof(double));
    gsl_interp_accel *spline_accel = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen,pcu_N);
    gsl_spline_init(spline,pcu_prob,pcu_tau,pcu_N);

    if(DEBUG_NU_MOMENTA) printf("#tau_t_eV: momenta [eV]:");
    
    for(int t=0; t<N_tau; t++){
      double prob = (0.5+t) / N_tau;
      tau_table_eV[t] = gsl_spline_eval(spline,prob,spline_accel);
      if(DEBUG_NU_MOMENTA) printf(" %g",tau_table_eV[t]);
    }

    if(DEBUG_NU_MOMENTA){ printf("\n"); fflush(stdout); }
    gsl_spline_free(spline);
    gsl_interp_accel_free(spline_accel);
    free(pcu_tau);
    free(pcu_prob);
    init = 1;
  }

  if(t == FREE_TAU_TABLE){
    free(tau_table_eV);
    init = 0;
    return 0;
  }
  return tau_table_eV[t];
}

double tau_t_eV(int t){ char dum; return tau_t_eV_init(t,&dum,0); }

//speed -tau_t / tau0_t of each neutrino species
double v_t_eta(int t, double eta, const struct cosmoparam C){
  double t_ma = tau_t_eV(t) / ( C.m_hdm_eV * aeta_in*exp(eta) );
  return (t_ma<1 ? t_ma : 1);
}

double v2_t_eta(int t, double eta, const struct cosmoparam C){
  double vt = v_t_eta(t,eta,C);
  return vt*vt;
}

//density ratio rho_t(eta)/rho_t(eta_stop) * aeta^2 and its log deriv
double Ec_t_eta(int t, double eta){ return 1.0/ ( aeta_in*exp(eta) ); }

double dlnEc_t_eta(int t, double eta){ return -1.0; }

//relativistic versions of the above, for Hubble rate calculation
double v2_t_eta_REL(int t, double eta, const struct cosmoparam C){
  double m_aeta_tau = C.m_hdm_eV * aeta_in*exp(eta) / tau_t_eV(t);
  return 1.0 / (1.0 + m_aeta_tau*m_aeta_tau);
}

double v_t_eta_REL(int t, double eta, const struct cosmoparam C){
  return sqrt(v2_t_eta_REL(t,eta,C));
}

double Ec_t_eta_REL(int t, double eta, const struct cosmoparam C){
  double vt2 = v2_t_eta_REL(t,eta,C), aeta = aeta_in*exp(eta);
  if(1-vt2 < 1e-12){
    double ma_tau = C.m_hdm_eV * aeta / tau_t_eV(t);
    return sqrt(1.0 + ma_tau*ma_tau) / (aeta*ma_tau); 
  }
  return 1.0 / (aeta * sqrt(1.0 - vt2));
}

double dlnEc_t_eta_REL(int t, double eta, const struct cosmoparam C){
  return -1.0 - v2_t_eta_REL(t,eta,C);
}

//////////////////////////// HOMOGENEOUS COSMOLOGY /////////////////////////////

//a(eta)^2 * rho_de(eta) / rho_de_0 and its derivative
double Ec_de_eta(double eta, const struct cosmoparam C){
  double aeta = aeta_in * exp(eta);
  return pow(aeta,-1.0 - 3.0*(C.w0_eos_de + C.wa_eos_de)) *
    exp(3.0*C.wa_eos_de*(aeta-1.0));
}

double dlnEc_de_eta(double eta, const struct cosmoparam C){
  double aeta = aeta_in * exp(eta);
  return -1.0 - 3.0*(C.w0_eos_de + C.wa_eos_de) + 3.0*C.wa_eos_de*aeta;
}

//conformal hubble parameter
double Hc2_Hc02_eta(double eta, const struct cosmoparam C){

  //scale factor
  double aeta = aeta_in * exp(eta), aeta2 = aeta*aeta, Ec_de = Ec_de_eta(eta,C);

  //sum Omega_{t,0} aeta^2 rho_t(eta)/rho_t_0 over CDM, photons, and DE
  double sum_OEc = C.Omega_cb_0/aeta + C.Omega_rel_0/aeta2 + C.Omega_de_0*Ec_de;
  
  //hot dark matter, using relativistic Omega_hdm(eta)
  for(int t=0; t<N_tau; t++) sum_OEc += C.Omega_hdm_t_0 * Ec_t_eta_REL(t,eta,C);

  return sum_OEc;
}

double Hc_eta(double eta, const struct cosmoparam C){
  return Hc0h * sqrt(Hc2_Hc02_eta(eta,C));
}

//d log(Hc) / d eta
double dlnHc_eta(double eta, const struct cosmoparam C){
  
  double aeta = aeta_in*exp(eta), aeta2 = aeta*aeta;
  double pre = 1.0 / ( 2.0 * Hc2_Hc02_eta(eta,C) );
  
  double sum_OdEc = -(1.0 + 3.0*C.w_eos_cdm) *  C.Omega_cb_0/aeta //CDM
    - (1.0 + 3.0*C.w_eos_gam) * C.Omega_rel_0/aeta2 //photons + massless nu
    + dlnEc_de_eta(eta,C) * C.Omega_de_0 * Ec_de_eta(eta,C); //DE
  
  for(int t=0; t<N_tau; t++)//neutrino fluids
    sum_OdEc +=  dlnEc_t_eta_REL(t,eta,C)*Ec_t_eta_REL(t,eta,C)*C.Omega_hdm_t_0;
  
  return pre * sum_OdEc;
}

//density fraction in spatially-flat universe
double OF_eta(int F, double eta, const struct cosmoparam C){
  
  double Hc02_Hc2 = 1.0/Hc2_Hc02_eta(eta,C), aeta = aeta_in*exp(eta);

  if(F == N_tau) //CDM
    return C.Omega_cb_0 * pow(aeta,-1.0-3.0*C.w_eos_cdm) * Hc02_Hc2;
  else if(F == N_tau+1) //photons + massless nu
    return C.Omega_rel_0 * pow(aeta,-1.0-3.0*C.w_eos_gam) * Hc02_Hc2;
  else if(F == N_tau+2) //dark energy, assumed Lambda
    return C.Omega_de_0 * Ec_de_eta(eta,C) * Hc02_Hc2;
  else if(F<0 || F>N_tau+2) return 0.0; //no fluids should have these indices
  return C.Omega_hdm_t_0 * Ec_t_eta(F,eta) * Hc02_Hc2;
}

/////////////////////////// INHOMOGENEOUS COSMOLOGY ////////////////////////////

//Poisson equation for Phi; use linear or nonlinear delta_cb
double Poisson_lin(double eta, int ik, const double *y,
                   const struct cosmoparam C){
  double k = KMIN * exp(DLNK * ik);
  double Hc2 = Hc0h2 * Hc2_Hc02_eta(eta,C), pre = -1.5 * Hc2 / (k*k);
  double sum_Od = OF_eta(N_tau,eta,C) * y[2*N_tau*N_mu*NK + ik];
  for(int t=0; t<N_tau; t++) sum_Od += OF_eta(t,eta,C) * y[2*t*N_mu*NK + ik];
  return pre * sum_Od;
}

double Poisson_nonlin(double eta, int ik, const double *y,
                      const struct cosmoparam C){
  double k = KMIN * exp(DLNK * ik), ee = exp(eta);
  double Hc2 = Hc0h2 * Hc2_Hc02_eta(eta,C), pre = -1.5 * Hc2 / (k*k);
  double dcb = ee * sqrt(  exp( y[2*N_tau*N_mu*NK + 2*NK + ik] ) );
  double sum_Od = OF_eta(N_tau,eta,C) * dcb;
  for(int t=0; t<N_tau; t++) sum_Od += OF_eta(t,eta,C) * y[2*t*N_mu*NK + ik];
  return pre * sum_Od;
}

//Eisenstein-Hu no-wiggle transfer function
double T_EH(double k, const struct cosmoparam C){
  double G_eff = C.Omega_m_0 * C.h *
    ( C.alpha_G + (1.0-C.alpha_G)/(1.0 + pow(0.43*k*C.sound_horiz,4)) );
  double q_EH = k * C.Theta_CMB_27_Sq / G_eff;
  double L_EH = log(2.0*M_E + 1.8*q_EH);
  double C_EH = 14.2 + 731.0/(1.0+62.5*q_EH);
  return L_EH / (L_EH + q_EH*q_EH*C_EH);
}

//transfer function and power spectrum interpolation
//  Strictly speaking, each time a function is called, I should compare
//  C to the stored value before returning the result from the existing
//  spline.  In order to save time, skip that for now, since I only intend
//  to use the code for a single cosmological model at a time.

#define NMAX_TRANSFER_INTERP (10000)

//interpolate total matter power spectrum from CAMB file
double Tmat0(double k, const struct cosmoparam C){

  static int init = 0;
  static double kTmin, kTmax;
  static gsl_interp_accel *acc;
  static gsl_spline *spl_T_Teh_lnk; //spline T/T_EH vs log(k)

  if(!init){

    int nt = 0;
    double lnkT[NMAX_TRANSFER_INTERP], Ttot_Teh[NMAX_TRANSFER_INTERP], t[13], T0;

    FILE *fp;
    if( (fp=fopen(C.file_transfer_function,"r")) == NULL ){
      printf("ERROR: Tmat0: Could not read transfer file %s. Quitting.\n",
             C.file_transfer_function);
      exit(1);
    }

    char line[10000];

    while( fgets(line,sizeof line, fp) && !feof(fp)){ 
      if(*line != '#'){
        sscanf(line,"%lg %lg %lg %lg %lg %lg %lg",t,t+1,t+2,t+3,t+4,t+5,t+6);
        //sscanf(line,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
        //       t,t+1,t+2,t+3,t+4,t+5,t+6,t+7,t+8,t+9,t+10,t+11,t+12);
        double Teh = T_EH(t[0],C);
        if(nt==0) T0 = t[6];
        lnkT[nt] = log(t[0]);
        Ttot_Teh[nt] = t[6] / (Teh * T0);
        nt++;
      }
    }

    kTmin = exp(lnkT[0])    * 1.000000001;
    kTmax = exp(lnkT[nt-1]) * 0.999999999;

    acc = gsl_interp_accel_alloc();
    spl_T_Teh_lnk = gsl_spline_alloc(gsl_interp_cspline,nt);
    gsl_spline_init(spl_T_Teh_lnk, lnkT, Ttot_Teh, nt);
    
    init = 1;
  }

  double kt = k;
  if(kt<kTmin) kt = kTmin;
  if(kt>kTmax) kt = kTmax; 
  
  return T_EH(k,C) * gsl_spline_eval(spl_T_Teh_lnk, log(kt), acc);;
}

//total matter power spectrum at z=0
double integrand_Pmat0(double lnkR, void *input){
  struct cosmoparam *C = (struct cosmoparam *)input;
  double R = 8.0, kR = exp(lnkR), k = kR/R, T = Tmat0(k,*C), W = 1.0-0.1*kR*kR;
  if(kR > 1e-2){
    double kR2 = kR*kR, kR3 = kR*kR2;
    W = 3.0 * (sin(kR)/kR3 - cos(kR)/kR2);
  }
  return pow(k,3.0+C->n_s) * T*T * W*W / (2.0 * M_PI*M_PI);
}

double Pmat0(double k, const struct cosmoparam C){

  static int init = 0;
  static double norm = 0;

  if(!init){
    struct cosmoparam par;
    copy_cosmoparam(C,&par);
    double s82U, err, x0=-15, x1=15;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function = &integrand_Pmat0;
    F.params = &par;
    gsl_integration_qag(&F, x0, x1, COSMOFUNC_EABS, COSMOFUNC_EREL,
			1000, 6, w, &s82U, &err);
    gsl_integration_workspace_free(w);
    norm = C.sigma_8 * C.sigma_8 / s82U;
    init = 1;

    //TESTING!!
    printf("#Pmat0: found norm = %g\n",norm);
    fflush(stdout);
    
  }

  double T = Tmat0(k,C);
  return norm * pow(k,C.n_s) * T*T;
}

