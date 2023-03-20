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

#define COSMOPARAM_MAX_REDSHIFTS (1000)
#define COSMOPARAM_MAX_CHAR_LEN (1000)
#define COSMOPARAM_NU_EFF (3.044)
#define COSMOPARAM_NU_MASSIVE (0.0)
#define COSMOPARAM_TAU_MAX_Q (100000)

const int COSMOPARAM_DEBUG_INIT = 1;

struct cosmoparam{

  //user-defined parameters
  double n_s;
  double sigma_8;
  double h;
  double Omega_m_0;
  double Omega_b_0;
  double Omega_hdm_0;
  double T_CMB_0_K;
  double w0_eos_de;
  double wa_eos_de;

  //code switches: 0 or 1
  int switch_nonlinear;
  int switch_1loop;
  int switch_print_linear;
  int switch_print_rsd;

  //inputs and outputs
  double z_nonlinear_initial;
  int num_z_outputs;
  double z_outputs[COSMOPARAM_MAX_REDSHIFTS];
  char file_transfer_function[COSMOPARAM_MAX_CHAR_LEN];
  int num_hdm_approx;
  double T_hdm_0_K;
  double m_hdm_eV;
  double g_hdm;
  char file_hdm_distribution[COSMOPARAM_MAX_CHAR_LEN];
  char file_nu_transfer_root[COSMOPARAM_MAX_CHAR_LEN];
  int num_interp_redshifts;
  double z_interp_redshifts[COSMOPARAM_MAX_REDSHIFTS];

  //fixed or derived parameters
  double Omega_cb_0;
  double Omega_hdm_t_0;
  double Omega_gam_0;
  double Omega_hdmrel_0;
  double Omega_hdmgam_0;
  double Omega_rel_0;
  double Omega_de_0;
  double Omega_m_h2_0;
  double Omega_b_h2_0;
  
  int N_tau;
  double N_nu_eff;
  double N_nu_massive;
  double f_cb_0;
  double f_hdm_0;
  
  double w_eos_cdm;
  double w_eos_gam;

  double alpha_G;
  double sound_horiz;
  double Theta_CMB_27_Sq;
};

int print_cosmoparam(const struct cosmoparam C, int verbosity){
  
  if(verbosity > 0){
    printf("#cosmoparam: n_s=%g, sigma_8=%g, h=%g, Omega_m_0=%g, Omega_b_0=%g, T_CMB_0_K=%g, w0_eos_de=%g, wa_eos_de=%g, T_hdm_0_K=%g, m_hdm_eV=%g, g_hdm=%g\n",
	   C.n_s, C.sigma_8, C.h, C.Omega_m_0, C.Omega_b_0, 
	   C.T_CMB_0_K, C.w0_eos_de, C.wa_eos_de,
	   C.T_hdm_0_K, C.m_hdm_eV, C.g_hdm);
    fflush(stdout);
  }

  if(verbosity > 1){
    printf("#cosmoparam: switch_nonlinear=%i, switch_1loop=%i, switch_print_linear=%i, switch_print_rsd=%i\n",
      C.switch_nonlinear, C.switch_1loop, C.switch_print_linear, 
      C.switch_print_rsd);
    fflush(stdout);
  }

  if(verbosity > 2){
    printf("#cosmoparam: z_nonlinear_initial=%g\n",C.z_nonlinear_initial);
    printf("#cosmoparam: z_outputs[%i]:", C.num_z_outputs);
    for(int i=0; i<C.num_z_outputs; i++) printf(" %g",C.z_outputs[i]);
    printf("\n");
    fflush(stdout);
  }

  return 0;
}

double integrate_hdm_distribution_function(struct cosmoparam *C){
  FILE *fp;
  if( (fp=fopen(C->file_hdm_distribution,"r")) == NULL ){
    printf("ERROR: Distribution function file %s not found.  Quitting.\n",
	   C->file_hdm_distribution);
    exit(1);
  }

  int N = 0;
  char buf[1000];
  double *q_hdm = malloc(COSMOPARAM_TAU_MAX_Q * sizeof(double));
  double *fqq_hdm = malloc(COSMOPARAM_TAU_MAX_Q * sizeof(double));

  while(fgets(buf, sizeof buf, fp) && N<COSMOPARAM_TAU_MAX_Q){
    double f_hdm;
    while(*buf=='#' || *buf=='\n'){ fgets(buf, sizeof buf, fp); }
    sscanf(buf,"%lg %lg", q_hdm+N, &f_hdm);
    fqq_hdm[N] = f_hdm*q_hdm[N]*q_hdm[N];
    N++;
  }

  if(N >= COSMOPARAM_TAU_MAX_Q){
    printf("#integrate_hdm_distribution_function: WARNING! Maximum number \n");
    printf("#  of points %i reached, stopping at q=%g.\n",
	   COSMOPARAM_TAU_MAX_Q, q_hdm[N-1]);
    fflush(stdout);
  }
  
  double result = ncint_composite_trapezoid(N,q_hdm,fqq_hdm);
  free(q_hdm);
  free(fqq_hdm);
  return result;
} 

int initialize_cosmoparam(struct cosmoparam *C, const char *params, int N_tau){
  FILE *fp;
  if( (fp=fopen(params,"r")) == NULL ){
    printf("ERROR: File %s not found.  Quitting.\n",params);
    exit(1);
  }

  char buf[1000], buf2[1000], *pbuf = buf;
  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%lg",&C->n_s);

  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%lg",&C->sigma_8);
  
  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%lg",&C->h);
  
  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%lg",&C->Omega_m_0);
  
  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%lg",&C->Omega_b_0);
  
  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%lg",&C->Omega_hdm_0);
  
  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%lg",&C->T_CMB_0_K);
  
  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%lg",&C->w0_eos_de);
  
  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%lg",&C->wa_eos_de);

  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%i",&C->switch_nonlinear);

  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%i",&C->switch_1loop);

  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%i",&C->switch_print_linear);

  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%i",&C->switch_print_rsd);

  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%lg",&C->z_nonlinear_initial);

  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%i",&C->num_z_outputs);

  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  for(int i=0; i<C->num_z_outputs; i++){
    sscanf(pbuf,"%s",buf2);
    sscanf(buf2,"%lg",&C->z_outputs[i]);
    pbuf += strlen(buf2)+1;
  }
    
  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%s",C->file_transfer_function);
  
  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%i",&C->num_hdm_approx);
  if(C->num_hdm_approx != 1){
    printf("ERROR: num_hdm_approx != 1.  Only mflr supported.\n");
    exit(1);
  }

  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%lg",&C->T_hdm_0_K);
  
  //do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  //sscanf(buf,"%lg",&C->m_hdm_eV);

  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%lg",&C->g_hdm);

  do{ fgets(buf, sizeof buf, fp); } while(*buf=='#' || *buf=='\n');
  sscanf(buf,"%s",C->file_hdm_distribution);

  //set fixed/derived parameters
  C->N_tau = N_tau;
  C->N_nu_eff = COSMOPARAM_NU_EFF;
  C->N_nu_massive = COSMOPARAM_NU_MASSIVE;

  C->Omega_cb_0 = C->Omega_m_0 - C->Omega_hdm_0;
  double nu_hdm_K3 = C->g_hdm * C->T_hdm_0_K*C->T_hdm_0_K*C->T_hdm_0_K
    * integrate_hdm_distribution_function(C);
  //C->Omega_hdm_0 = 3.9931e-4 * nu_hdm_K3 * C->m_hdm_eV / (C->h*C->h); 
  C->m_hdm_eV = C->Omega_hdm_0*C->h*C->h / (0.09905 * nu_hdm_K3);

  //TESTING!!
  printf("#initialize_cosmoparam: nu_hdm_K^3=%1.14g, Omega_hdm_0*h^2=%1.14g\n",
	 nu_hdm_K3, C->Omega_hdm_0 * C->h * C->h);
  fflush(stdout);
  

  C->Omega_hdm_t_0 = C->Omega_hdm_0 / N_tau;
  C->Omega_gam_0 = 4.46911743913795e-07 * pow(C->T_CMB_0_K,4) / (C->h*C->h);
  C->Omega_hdmrel_0 = 0.227107317660239
    * (C->N_nu_eff-C->N_nu_massive) * C->Omega_gam_0;
  C->Omega_hdmgam_0 = (1.0+0.227107317660239*C->N_nu_eff)*C->Omega_gam_0;
  C->Omega_rel_0 = C->Omega_gam_0 + C->Omega_hdmrel_0;
  C->Omega_de_0 = 1.0 - C->Omega_cb_0 - C->Omega_hdm_0 - C->Omega_rel_0;
  C->Omega_m_h2_0 = C->Omega_m_0 * C->h	* C->h;
  C->Omega_b_h2_0 = C->Omega_b_0 * C->h * C->h;
  
  C->w_eos_cdm = 0;
  C->w_eos_gam = 0.333333333333333333;
  C->f_cb_0 = C->Omega_cb_0 / C->Omega_m_0;
  C->f_hdm_0 = C->Omega_hdm_0 / C->Omega_m_0;

  C->sound_horiz = 44.5*C->h*log(9.83/C->Omega_m_h2_0)
    / sqrt(1.0 + 10.0 * pow(C->Omega_b_h2_0,0.75)); /*[Mpc/h]*/ 
  double rbm = C->Omega_b_0 / C->Omega_m_0;
  C->alpha_G = 1.0 - 0.328*log(431.0*C->Omega_m_h2_0) * rbm
    + 0.38 * log(22.3*C->Omega_m_h2_0) * rbm*rbm;
  C->Theta_CMB_27_Sq = pow(C->T_CMB_0_K/2.7,2);
  
  if(COSMOPARAM_DEBUG_INIT) print_cosmoparam(*C,COSMOPARAM_DEBUG_INIT); 

  fclose(fp);
  return 0;
}
  
int copy_cosmoparam(const struct cosmoparam B, struct cosmoparam *C){

  C->n_s = B.n_s;
  C->sigma_8 = B.sigma_8;
  C->h = B.h;
  C->Omega_m_0 = B.Omega_m_0;
  C->Omega_b_0 = B.Omega_b_0;
  C->T_CMB_0_K = B.T_CMB_0_K;
  C->w0_eos_de = B.w0_eos_de;
  C->wa_eos_de = B.wa_eos_de;
  C->T_hdm_0_K = B.T_hdm_0_K;
  C->m_hdm_eV  = B.m_hdm_eV;
  C->g_hdm     = B.g_hdm;

  C->switch_nonlinear = B.switch_nonlinear;
  C->switch_1loop = B.switch_1loop;
  C->switch_print_linear = B.switch_print_linear;
  C->switch_print_rsd = B.switch_print_rsd;

  C->z_nonlinear_initial = B.z_nonlinear_initial;
  C->num_z_outputs = B.num_z_outputs;
  for(int i=0; i<B.num_z_outputs; i++) C->z_outputs[i] = B.z_outputs[i];
  strcpy(C->file_transfer_function,B.file_transfer_function);
  C->num_hdm_approx = B.num_hdm_approx;
  strcpy(C->file_hdm_distribution, B.file_hdm_distribution);
  strcpy(C->file_nu_transfer_root, B.file_nu_transfer_root);
  C->num_interp_redshifts = B.num_interp_redshifts;
  for(int i=0; i<B.num_interp_redshifts; i++)
    C->z_interp_redshifts[i] = B.z_interp_redshifts[i];
  
  //set fixed/derived parameters
  C->N_tau = B.N_tau;
  C->N_nu_eff = B.N_nu_eff;
  C->N_nu_massive = B.N_nu_massive;
  C->Omega_cb_0 = B.Omega_cb_0;
  C->Omega_hdm_0 = B.Omega_hdm_0;
  C->Omega_hdm_t_0 = B.Omega_hdm_t_0;
  C->Omega_gam_0 = B.Omega_gam_0;
  C->Omega_hdmrel_0 = B.Omega_hdmrel_0;
  C->Omega_hdmgam_0 = B.Omega_hdmgam_0;
  C->Omega_rel_0 = B.Omega_rel_0;
  C->Omega_de_0 = B.Omega_de_0;
  C->Omega_m_h2_0 = B.Omega_m_h2_0;
  C->Omega_b_h2_0 = B.Omega_b_h2_0;
  
  C->w_eos_cdm = 0;
  C->w_eos_gam = 0.333333333333333333;
  C->f_cb_0 = B.f_cb_0;
  C->f_hdm_0 = B.f_hdm_0;

  C->sound_horiz = B.sound_horiz;
  C->alpha_G = B.alpha_G;
  C->Theta_CMB_27_Sq = B.Theta_CMB_27_Sq;

  return 0;
}

double cosmoparam_fdiff(double x, double y){
  return 2.0 * fabs(x-y) / (fabs(x) + fabs(y) + 1e-100);
}

int isequal_cosmoparam(const struct cosmoparam B, const struct cosmoparam C){
  int equal = (B.switch_nonlinear == C.switch_nonlinear);
  equal = equal && ( B.switch_1loop == C.switch_1loop );
  equal	= equal	&& ( B.switch_print_linear == C.switch_print_linear );
  equal = equal && ( B.switch_print_rsd == C.switch_print_rsd );
  equal = equal && ( B.num_z_outputs == C.num_z_outputs );
  equal = equal && ( B.num_hdm_approx == C.num_hdm_approx );
  equal = equal && ( B.N_tau == C.N_tau );
  if(!equal) return 0;
  if(strcmp(B.file_transfer_function, C.file_transfer_function) != 0) return 0;
  if(strcmp(B.file_hdm_distribution, C.file_hdm_distribution) != 0) return 0;

  double fdmax = cosmoparam_fdiff(B.n_s,C.n_s);
  fdmax = fmax( fdmax, cosmoparam_fdiff(B.sigma_8,     C.sigma_8) );
  fdmax	= fmax( fdmax, cosmoparam_fdiff(B.h,           C.h) );
  fdmax = fmax( fdmax, cosmoparam_fdiff(B.Omega_m_0,   C.Omega_m_0) );
  fdmax = fmax( fdmax, cosmoparam_fdiff(B.Omega_b_0,   C.Omega_b_0) );
  fdmax = fmax( fdmax, cosmoparam_fdiff(B.Omega_hdm_0, C.Omega_hdm_0) );
  fdmax = fmax( fdmax, cosmoparam_fdiff(B.T_CMB_0_K,   C.T_CMB_0_K) );
  fdmax = fmax( fdmax, cosmoparam_fdiff(B.w0_eos_de,   C.w0_eos_de) );
  fdmax = fmax( fdmax, cosmoparam_fdiff(B.wa_eos_de,   C.wa_eos_de) );
  fdmax = fmax( fdmax, cosmoparam_fdiff(B.T_hdm_0_K,   C.T_hdm_0_K) );
  fdmax = fmax( fdmax, cosmoparam_fdiff(B.m_hdm_eV,    C.m_hdm_eV) );
  fdmax = fmax( fdmax, cosmoparam_fdiff(B.g_hdm,       C.g_hdm) );

  for(int i=0; i<B.num_z_outputs; i++)
    fdmax = fmax( fdmax, cosmoparam_fdiff(B.z_outputs[i],C.z_outputs[i]) );

  return (fdmax < 1e-6);
}
