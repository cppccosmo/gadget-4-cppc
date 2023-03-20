////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    Copyright 2022 Amol Upadhye                                             //
//                                                                            //
//    This file is part of MuFLR.                                             //
//                                                                            //
//    MuFLR is free software: you can redistribute it and/or modify           //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation, either version 3 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    MuFLR is distributed in the hope that it will be useful,                //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with MuFLR.  If not, see <http://www.gnu.org/licenses/>.          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

//useful functions: square, cube, quadrature add
double sq(double x){ return x*x; }
double cu(double x){ return x*x*x; }
inline double qu(double x){ return sq(sq(x)); }
inline double qadd(double x, double y){ return sqrt(sq(x)+sq(y)); }

//tilt applied to power spectrum before fft
#define FASTPT_TILT_NU (-2)
const int nu_int = FASTPT_TILT_NU;
const double nu = FASTPT_TILT_NU;

//power spectrum and Fourier coefficient window functions
inline double W_edge(double x){ return x - sin(2.0*M_PI*x)/(2.0*M_PI); }

double WP(double lnk){
  if(lnk <= LNK_PAD_WINLO) return 0;
  else if(lnk < LNK_PAD_WINLI) 
    return W_edge((lnk-LNK_PAD_WINLO)/(LNK_PAD_WINLI-LNK_PAD_WINLO));
  else if(lnk < LNK_PAD_WINRI) return 1;
  else if(lnk < LNK_PAD_WINRO) 
    return W_edge((LNK_PAD_WINRO-lnk)/(LNK_PAD_WINRO-LNK_PAD_WINRI));
  return 0;
}

double WC(int n){ //coeffs in GSL halfcomplex notation
  const int nl = 1*NKP/8, nc = NKP/2, nr = 7*NKP/8, Dn = 3*NKP/8;
  if(n<=nl || n>=nr) return 1;
  else if(n<nc) return W_edge((double)(nc-n)/Dn);
  else if(n<nr) return W_edge((double)(n-nc)/Dn);
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
// fast-pt computations of A, R, P_T 
// (see McEwen, Fang, Hirata, Blazek 1603.04826)

////////// functions derived from Gamma

int g_MFHB_dddd(double mu, double reKappa, double imKappa, double *g_out){

  //compute Gamma function
  gsl_sf_result lnGtop, lnGbot, argTop, argBot;
  gsl_sf_lngamma_complex_e(0.5*(mu+reKappa+1), 0.5*imKappa,  &lnGtop, &argTop);
  gsl_sf_lngamma_complex_e(0.5*(mu-reKappa+1), -0.5*imKappa, &lnGbot, &argBot);
  
  //g_out[] = { |g|, arg(g) }
  g_out[0] = exp(lnGtop.val - lnGbot.val);
  g_out[1] = argTop.val - argBot.val;
  return 0;
}

int f_MFHB_ddd(double reRho, double imRho, double *f_out){

  double g[2], pre = 0.5*sqrt(M_PI) * pow(2.0,reRho);
  g_MFHB_dddd(0.5,reRho-0.5,imRho,g);
  f_out[0] = pre * g[0];
  f_out[1] = imRho*M_LN2 + g[1];
  return 0;
}

//frontends for the above functions
int f_MFHB(int alpha, int beta, int h, double *f){
  int n = (h<=NKP ? h : h-2*NKP);
  return f_MFHB_ddd(-4.0-2.0*nu-(double)(alpha+beta),-2.0*M_PI*n/(DLNK*NKP),f);
}

int g_reg_MFHB(int m, double *g_out){ //regularized version for ell=0, alpha=-2
  int n = (m<=NKP/2 ? m : m-NKP);
  return f_MFHB_ddd(nu, 2.0*M_PI*n/(DLNK*NKP), g_out); 
}

int g_MFHB(int ell, int alpha, int m, double *g){
  if(m==0 && alpha==ell-nu_int){ g[0]=0; g[1]=0; return 0; }
  if(alpha==-2 && ell==0) return g_reg_MFHB(m,g);
  int n = (m<=NKP/2 ? m : m-NKP);
  return g_MFHB_dddd(0.5+(double)ell, 1.5+nu+(double)alpha,
		     2.0*M_PI*n/(DLNK*NKP), g);
}

////////// frontends for gsl fft routines

//forward fft, output replaces input array
int fft(double *x, int N){ return gsl_fft_real_radix2_transform(x,1,N); }

//inverse fft
int ifft(double *x, int N){ return gsl_fft_halfcomplex_radix2_inverse(x,1,N); }

//backward fft, identical to ifft except for lack of normalization factor
int bfft(double *x, int N){ return gsl_fft_halfcomplex_radix2_backward(x,1,N); }

//convolution fo real functions; assume arrays of equal length
int convolve(int N, double *in0, double *in1, double *out){
  fft(in0,N);
  fft(in1,N);

  //out is now in halfcomplex format
  out[0] = in0[0]*in1[0];
  out[N/2] = in0[N/2]*in1[N/2];

  for(int i=1; i<N/2; i++){
    out[i] = in0[i]*in1[i] - in0[N-i]*in1[N-i];
    out[N-i] = in0[i]*in1[N-i] + in0[N-i]*in1[i];
  }

  ifft(out,N);
  return 0;
}

//convolution for halfcomplex functions; assume arrays already of equal length 
int iconvolve(int N, double *in0, double *in1, double *out){
  ifft(in0,N);
  ifft(in1,N);
  for(int i=0; i<N; i++) out[i] = in0[i] * in1[i] * N;
  fft(out,N);
  return 0;
}

//convolution for complex arrays, with even elements respresenting
//real values and odd elements imaginary values
//For linear convolution, final half of input arrays should be zero-padded.
int cconvolve(int N, double *in0, double *in1, double *out){
  gsl_fft_complex_radix2_forward(in0,1,N);
  gsl_fft_complex_radix2_forward(in1,1,N);
  for(int i=0; i<N; i++){
    out[2*i] = in0[2*i]*in1[2*i] - in0[2*i+1]*in1[2*i+1];
    out[2*i+1] = in0[2*i+1]*in1[2*i] + in0[2*i]*in1[2*i+1];
  }
  gsl_fft_complex_radix2_inverse(out,1,N);
  return 0;
}

int convolve_bruteforce(int N, double *in0, double *in1, double *out){
  for(int i=0; i<N; i++) out[i] = 0;

  for(int n=0; n<N; n++){
    for(int m=0; m<=n; m++) out[n] += in0[m] * in1[n-m];
    for(int m=n+1; m<N; m++) out[n] += in0[m] * in1[N+m-n];
  }

  return 0;
}

//regularized J_{2,-2,0} from McEwen++ 1603.04826
int Jreg_MFHB(const double *Palpha, const double *Pbeta, double *Ji){
  const int alpha=2, beta=-2, ell=0;

  //fft, after multiplying by power law 
  double ca[NKP], cb[NKP];
  for(int i=0; i<NKP; i++){
    double lnk = LNK_PAD_MIN + DLNK*i, k_nnu = exp(-nu*lnk);
    ca[i] = Palpha[i] * k_nnu;
    cb[i] = Pbeta[i] * k_nnu;
  }
  fft(ca,NKP);
  fft(cb,NKP);
  for(int i=0; i<NKP; i++){
    double win = WC(i);
    ca[i] *= win;
    cb[i] *= win;
  }

  double cga[4*NKP], cgb[4*NKP], ga[2], gb[2];
  for(int i=0; i<4*NKP; i++){ cga[i]=0; cgb[i]=0; }

  g_MFHB(ell,alpha,0,ga);
  g_MFHB(ell,beta,0,gb);
  ga[0] *= pow(2.0,1.5+nu+alpha);
  ga[1] += 2.0*M_PI     * 0     /(DLNK*NKP) * M_LN2;
  cga[0] = ca[0]*ga[0]*cos(ga[1]);
  cga[1] = ca[0]*ga[0]*sin(ga[1]);
  cgb[0] = cb[0]*gb[0]*cos(gb[1]);
  cgb[1] = cb[0]*gb[0]*sin(gb[1]);

  //halfcomplex convolution
  for(int i=1; i<NKP/2; i++){
    g_MFHB(ell,alpha,i,ga);
    g_MFHB(ell,beta,i,gb);
    ga[0] *= pow(2.0,1.5+nu+alpha);
    ga[1] += 2.0*M_PI*i/(DLNK*NKP)* M_LN2;

    //fullcomplex convolution for halfcomplex ca and cb
    cga[2*i] = ca[i]*ga[0]*cos(ga[1]) - ca[NKP-i]*ga[0]*sin(ga[1]);
    cga[2*i+1] = ca[i]*ga[0]*sin(ga[1]) + ca[NKP-i]*ga[0]*cos(ga[1]);
    cga[4*NKP-2*i] = ca[i]*ga[0]*cos(ga[1]) - ca[NKP-i]*ga[0]*sin(ga[1]);
    cga[4*NKP-2*i+1] = -ca[i]*ga[0]*sin(ga[1]) - ca[NKP-i]*ga[0]*cos(ga[1]);
    cgb[2*i] = cb[i]*gb[0]*cos(gb[1]) - cb[NKP-i]*gb[0]*sin(gb[1]);
    cgb[2*i+1] = cb[i]*gb[0]*sin(gb[1]) + cb[NKP-i]*gb[0]*cos(gb[1]);
    cgb[4*NKP-2*i] = cb[i]*gb[0]*cos(gb[1]) - cb[NKP-i]*gb[0]*sin(gb[1]);
    cgb[4*NKP-2*i+1] = -cb[i]*gb[0]*sin(gb[1]) - cb[NKP-i]*gb[0]*cos(gb[1]);
  }

  //convolve to get C_h, then IFFT
  double C_h_cmplx[4*NKP], C_h[2*NKP], Cf_h[2*NKP], f[2];
  cconvolve(2*NKP,cga,cgb,C_h_cmplx);
  
  //recover halfcomplex C_h
  C_h[0] = C_h_cmplx[0];
  C_h[NKP] = C_h_cmplx[2*NKP];
  for(int i=0; i<NKP; i++){
    C_h[i] = C_h_cmplx[2*i];
    C_h[2*NKP-i] = C_h_cmplx[2*i+1];
  }

  f_MFHB(alpha,beta,0,f);
  Cf_h[0] = C_h[0] * f[0] * cos(f[1]);

  for(int i=1; i<=NKP; i++){
    double MC = qadd(C_h[i],C_h[2*NKP-i]), AC = atan2(C_h[2*NKP-i],C_h[i]);
    if(i==NKP){ MC = C_h[NKP]; AC = 0; }
    f_MFHB(alpha,beta,i,f);

    double MCf = MC * f[0], ACf = AC + f[1];
    if(i==NKP) ACf = 0;
    Cf_h[2*NKP-i] = MCf * sin(ACf);
    Cf_h[i] = MCf * cos(ACf);
  }

  bfft(Cf_h,2*NKP); //now it's real

  //assemble J array
  double sl = (ell%2==0 ? 1 : -1);
  double pre = sl / (2.0*M_PI*M_PI*NKP*NKP) * sqrt(2.0/M_PI);
  for(int i=0; i<NKP; i++){
    double ki=exp(LNK_PAD_MIN+DLNK*i), k_npm2=pow(ki,3.0+2.0*nu+alpha+beta);
    Ji[i] = pre * k_npm2 * Cf_h[2*i];
  }

  return 0;
}

//J_{alpha,beta,ell} from McEwen++ 1603.04826
int J_MFHB(int alpha, int beta, int ell, 
           const double *Palpha, const double *Pbeta, double *Ji){

  //do we want the regularized version?
  if(ell==0 && alpha==2 && beta==-2) return Jreg_MFHB(Palpha,Pbeta,Ji);
  if(ell==0 && alpha==-2 && beta==2) return Jreg_MFHB(Pbeta,Palpha,Ji);

  //fft, after multiplying by power law
  double ca[NKP], cb[NKP];
  for(int i=0; i<NKP; i++){
    double lnk = LNK_PAD_MIN + DLNK*i, k_nnu = exp(-nu*lnk);
    ca[i] = Palpha[i] * k_nnu;
    cb[i] = Pbeta[i] * k_nnu;
  }
  fft(ca,NKP);
  fft(cb,NKP);
  for(int i=0; i<NKP; i++){ 
    double win = WC(i);
    ca[i] *= win;
    cb[i] *= win;
  }

  //combine c with g, into 2N element array
  double cga[2*NKP], cgb[2*NKP], ga[2], gb[2];
  for(int i=0; i<2*NKP; i++){ cga[i]=0; cgb[i]=0; }

  g_MFHB(ell,alpha,0,ga);
  g_MFHB(ell,beta,0,gb);
  cga[0] = ca[0]*ga[0];
  cgb[0] = cb[0]*gb[0];

  for(int i=1; i<NKP/2; i++){
    g_MFHB(ell,alpha,i,ga);
    g_MFHB(ell,beta,i,gb);

    //use for fft convolution
    cga[i] = ca[i]*ga[0]*cos(ga[1]) - ca[NKP-i]*ga[0]*sin(ga[1]);
    cga[2*NKP-i] = ca[i]*ga[0]*sin(ga[1]) + ca[NKP-i]*ga[0]*cos(ga[1]);
    cgb[i] = cb[i]*gb[0]*cos(gb[1]) - cb[NKP-i]*gb[0]*sin(gb[1]);
    cgb[2*NKP-i] = cb[i]*gb[0]*sin(gb[1]) + cb[NKP-i]*gb[0]*cos(gb[1]);
  }

  //convolve to get C_h, then IFFT
  double C_h[2*NKP], Cftau_h[2*NKP], f[2];
  iconvolve(2*NKP,cga,cgb,C_h);

  f_MFHB(alpha,beta,0,f);
  Cftau_h[0] = C_h[0] * f[0] * cos(f[1]);
  
  for(int i=1; i<=NKP; i++){
    double MC = qadd(C_h[i],C_h[2*NKP-i]), AC = atan2(C_h[2*NKP-i],C_h[i]);
    if(i==NKP){ MC = C_h[NKP]; AC = 0; }
    f_MFHB(alpha,beta,i,f);

    double tau = 2.0*M_PI*i/(DLNK*NKP);
    double MCftau = MC * f[0], ACftau = AC + f[1] + M_LN2*tau;
    Cftau_h[2*NKP-i] = MCftau * sin(ACftau);
    Cftau_h[i] = MCftau * cos(ACftau);
  }

  bfft(Cftau_h,2*NKP); //now it's real

  //assemble J array
  double sl = (ell%2==0 ? 1 : -1);
  double pre = sl / (2.0*M_PI*M_PI*NKP*NKP);
  for(int i=0; i<NKP; i++){
    double ki=exp(LNK_PAD_MIN+DLNK*i), k2_npm2=pow(ki*2,3.0+2.0*nu+alpha+beta);
    Ji[i] = pre * k2_npm2 * Cftau_h[2*i];
  }

  return 0;
}

const int tZ = 10; //number of Taylor expansion terms to keep
const double epsZ = 1e-2; //switch to expansion for r<epsZ or r>1/epsZ
 
double Zreg_n(int n, double r){

  if(n<0) return Zreg_n(-n,1.0/r);
  
  double Z = 0, lnkq = log(fabs((1.0+r)/(1.0-r)));

  switch(n){
  case 0:
    Z = 1.0;
    break;
  case 1:
    if(r < epsZ){
      for(int m=0; m<tZ; m++) Z += 2.0*pow(r,2.0*m+1.0)*(1.0-r) / (2.0*m+1.0); }
    else if(r > 1.0/epsZ){
      for(int m=0; m<tZ; m++) Z += 2.0*pow(r,-2.0*m-1.0)*(1.0-r)/(2.0*m+1.0); }
    else if(r == 1){ Z = 0.0; }
    else Z = (1.0 - r) * lnkq;
    break;
  case 2:
    if( r < epsZ){
      Z = 2.0*r;
      for(int m=0; m<tZ; m++) 
        Z += 2.0*pow(r,2.0*m+3.0) / ((2.0*m+1.0)*(2.0*m+3.0));
    }
    else if(r > 1.0/epsZ){
      for(int m=0; m<tZ; m++)
        Z += 2.0*pow(r,-2.0*m-1.0) / ((2.0*m+1.0)*(2.0*m+3.0));
    }
    else if(r == 1){ Z = 1.0; }
    else Z = r + 0.5*(1.0-r*r)*lnkq;
    break;
  case 3:
    if( r < epsZ){
      Z = r*r;
      for(int m=0; m<tZ; m++)
        Z += (1.0-cu(r))*pow(r,2*m+1) / (2.0*m+1.0);
    }
    else if(r > 1.0/epsZ){
      for(int m=0; m<tZ; m++)
        Z += pow(r,-2*m) *((2.0*m+3.0)/r-2.0*m-1.0) / ((2.0*m+1.0)*(2.0*m+3.0));
    }
    else if(r == 1){ Z = 1.0; }
    else Z = sq(r) + 0.5*(1.0-cu(r))*lnkq;
    break;
  case 4:
    if( r < epsZ){
      Z = (4.0/3.0) * (r + cu(r));
      for(int m=0; m<tZ; m++)
        Z += -4.0*pow(r,2*m+5) / ((2.0*m+1.0)*(2.0*m+5.0));
    }
    else if(r > 1.0/epsZ){
      for(int m=0; m<tZ; m++)
        Z += 4.0 * pow(r,-2*m-1) / ((2.0*m+1.0)*(2.0*m+5.0));
    }
    else if(r == 1){ Z = 4.0/3.0; }
    else Z = cu(r) + r/3.0 + 0.5*(1.0-sq(sq(r)))*lnkq;
    break;
  case 5:
    if( r < epsZ){
      Z = sq(sq(r)) + sq(r)/3.0;
      for(int m=0; m<tZ; m++)
        Z += (1.0-cu(r)*sq(r))*pow(r,2*m+1) / (2.0*m + 1.0);
    }
    else if(r > 1.0/epsZ){
      for(int m=0; m<tZ; m++)
        Z += pow(r,-2*m) *((2.0*m+5.0)/r-2.0*m-1.0) / ((2.0*m+1.0)*(2.0*m+5.0));
    }
    else if(r == 1){ Z = 4.0/3.0; }
    else Z = sq(sq(r)) + sq(r)/3.0 + 0.5*(1.0-cu(r)*sq(r))*lnkq;
    break;
  default:
    printf("#ERROR in Zreg_n for n=%i: ",n);
    printf("kernel terms only defined for |n|<=5.  Quitting.");
    exit(1);
    break;
  }

  return Z;
}

int PZ_reg(int n, const double *Pq, const double *Pk, double *PZn){

  //for s[m]=log(q_m): Fs[m] = Pq(q_m) and Gs[m] = q_m^{-3}*Z(1/q_m)
  double Fs[4*NKP], Gs[4*NKP], FGconv[4*NKP];
  for(int i=0; i<4*NKP; i++){ Fs[i]=(i<NKP ? Pq[i] : 0);   Gs[i] = 0; }
  
  for(int i=0; i<NKP; i++){

    //r>1
    double si=DLNK*(i-NKP), r=exp(-si), r2=sq(r), r3=r*r2, r4=sq(r2), r5=r*r4;
    double Zi = Zreg_n(n,r);
    Gs[i] = Zi * r3;
  }

  Gs[NKP] = Zreg_n(n,1.0); //r=1   //12.0 + 10.0 + 100.0 - 42.0;

  for(int i=NKP+1; i<2*NKP; i++){
    //r<1
    double si=DLNK*(i-NKP), r=exp(-si), r2=sq(r), r3=r*r2, r4=sq(r2), r5=r*r4;
    double Zi = Zreg_n(n,r);
    Gs[i] = Zi * r3;
  }

  convolve_bruteforce(4*NKP,Fs,Gs,FGconv);
  //convolve(4*NKP,Fs,Gs,FGconv);

  //double pre = DLNK / (252.0 * sq(2.0*M_PI));
  double pre = DLNK / (2.0 * sq(M_PI));
  for(int i=0; i<NKP; i++){
    double lnk = LNK_PAD_MIN + DLNK*i, k = exp(lnk), k3=k*k*k;
    double kfac = k3 * Pk[i];
    PZn[i] = pre * kfac * FGconv[i+NKP];
  }

  return 0;
}

const int nJ = 63, nJn0 = 63;

const int ell_n[7]   = {0, 0, 1, 2, 2, 3, 4};
const int alpha_n[7] = {0, 2, 1, 0, 2, 1, 0};

const int elln0_n[7]   = {0, 2, 4, 0, 2, 4, 6};
const int alphan0_n[7] = {0, 0, 0, 2, 2, 2, 2};
const int betan0_n[7]  = {2, 2, 2, 2, 2, 2, 2};

const int Z_n[7] = {0, 1, -1, 3, -3, 5, -5};

int compute_Aacdbef(double eta, const double *Ppad, double *Aacdbef){
  
  //initialize to zero; only compute nonzero components
  for(int i=0; i<64*NK; i++) Aacdbef[i] = 0;
  
  //extrapolate power spectra.  P ~ k^{n_s} at low k, Eisenstein-Hu at high k
  double P[3*NKP];
  for(int i=0; i<NKP; i++){
    double k = exp(LNK_PAD_MIN + DLNK*i), Win = WP(LNK_PAD_MIN + DLNK*i);
    P[i] = Ppad[i] * Win;
    P[i+NKP] = Ppad[i+NKP] * Win;
    P[i+2*NKP] = Ppad[i+2*NKP] * Win;
  }

  //full TRG calculation with input P; do improved 1-loop later.
  double J[nJ][NKP], Jn0[nJn0][NKP], PZ[nJ][NKP];

#pragma omp parallel for schedule(dynamic)
  for(int iJ=0; iJ<nJ; iJ++){
    int n = iJ/9, iabcd = iJ%9, iab = iabcd/3, icd = iabcd%3;
    J_MFHB(alpha_n[n],-alpha_n[n],ell_n[n], &P[iab*NKP],&P[icd*NKP], J[iJ]);
    //PZ_reg(Z_n[n], &P[iab*NKP],&P[icd*NKP], PZ[iJ]);
  }

#pragma omp parallel for schedule(dynamic)
  for(int iJ=0; iJ<nJ; iJ+=3){
    int n = iJ/9, iabcd = iJ%9, iab = iabcd/3, icd = 0;
    PZ_reg(Z_n[n], &P[iab*NKP],&P[icd*NKP], PZ[iJ]);

    for(int i=0; i<NKP; i++){
      PZ[iJ+1][i] = PZ[iJ][i] * P[1*NKP+i] / (P[0*NKP+i] + 1e-100);
      PZ[iJ+2][i] = PZ[iJ][i] * P[2*NKP+i] / (P[0*NKP+i] + 1e-100);
    }
  }

#pragma omp parallel for schedule(dynamic)
  for(int iJ=0; iJ<nJn0; iJ++){
    int n = iJ/9, iabcd = iJ%9, iab = iabcd/3, icd = iabcd%3;
    J_MFHB(alphan0_n[n], betan0_n[n], elln0_n[n], 
	   &P[iab*NKP], &P[icd*NKP], Jn0[iJ]);
  }

#pragma omp parallel for schedule(dynamic)
  for(int i=0; i<NK; i++){
    double k = exp(LNKMIN + DLNK*i), k2 = k*k, 
      pre_A = k / (4.0*M_PI), pre_R = 1.0/(2.0*M_PI*k);
    double Jterms = 0, PZterms = 0;

    //A_{acd,bef}
    Jterms = J[4*9+1][NSHIFT+i]/6 
      +J[2*9+1][NSHIFT+i]/2 
      +J[0*9+1][NSHIFT+i]/4 
      +J[1*9+1][NSHIFT+i]/12 
      + J[3*9+3][NSHIFT+i]/6 
      + J[2*9+3][NSHIFT+i]/4 
      + J[2*9+1][NSHIFT+i]/4 
      + J[0*9+3][NSHIFT+i]/3;
    PZterms = -PZ[0*9+1][NSHIFT+i]/12.0 
      + (PZ[4*9+3][NSHIFT+i] 
	 - PZ[2*9+3][NSHIFT+i] 
	 + PZ[0*9+3][NSHIFT+i] 
	 + PZ[1*9+3][NSHIFT+i]/2
         - PZ[3*9+1][NSHIFT+i] 
	 + PZ[1*9+1][NSHIFT+i] 
	 + PZ[0*9+1][NSHIFT+i]*3 
	 - PZ[2*9+1][NSHIFT+i]/2) / 16;
    Aacdbef[8*NK+i] = pre_A * (Jterms + PZterms);

    Jterms = J[4*9+2][NSHIFT+i]/6 
      + J[2*9+2][NSHIFT+i]/2 
      +J[0*9+2][NSHIFT+i]/4 
      +J[1*9+2][NSHIFT+i]/12
      + J[3*9+4][NSHIFT+i]/6 
      + J[2*9+4][NSHIFT+i]/4 
      + J[2*9+4][NSHIFT+i]/4 
      + J[0*9+4][NSHIFT+i]/3;
    PZterms = 0;
    Aacdbef[9*NK+i] = pre_A * (Jterms + PZterms);

    Jterms = J[4*9+4][NSHIFT+i]/6 
      +J[2*9+4][NSHIFT+i]/2 
      +J[0*9+4][NSHIFT+i]/4 
      +J[1*9+4][NSHIFT+i]/12
      + J[3*9+6][NSHIFT+i]/6 
      + J[2*9+6][NSHIFT+i]/4 
      + J[2*9+2][NSHIFT+i]/4 
      + J[0*9+6][NSHIFT+i]/3;
    PZterms = -PZ[0*9+4][NSHIFT+i]/12.0
      + (PZ[4*9+6][NSHIFT+i] 
	 - PZ[2*9+6][NSHIFT+i] 
	 + PZ[0*9+6][NSHIFT+i] 
	 + PZ[1*9+6][NSHIFT+i]/2
         - PZ[3*9+4][NSHIFT+i] 
	 + PZ[1*9+4][NSHIFT+i] 
	 + PZ[0*9+4][NSHIFT+i]*3 
	 - PZ[2*9+4][NSHIFT+i]/2) / 16;
    Aacdbef[10*NK+i] = pre_A * (Jterms + PZterms);

    Jterms = J[4*9+5][NSHIFT+i]/6 
      +J[2*9+5][NSHIFT+i]/2 
      +J[0*9+5][NSHIFT+i]/4 
      +J[1*9+5][NSHIFT+i]/12
      + J[3*9+7][NSHIFT+i]/6 
      + J[2*9+7][NSHIFT+i]/4 
      + J[2*9+5][NSHIFT+i]/4 
      + J[0*9+7][NSHIFT+i]/3;
    PZterms = 0;
    Aacdbef[11*NK+i] = pre_A * (Jterms + PZterms);

    Jterms = J[5*9+4][NSHIFT+i]/5 
      + J[3*9+4][NSHIFT+i]/2 
      + J[4*9+4][NSHIFT+i]/6 
      + 0.55*J[2*9+4][NSHIFT+i]
      + J[2*9+4][NSHIFT+i]/4 
      +J[0*9+4][NSHIFT+i]/4 
      + J[1*9+4][NSHIFT+i]/12;
    PZterms = -PZ[0*9+2][NSHIFT+i]/12.0
      + (PZ[4*9+4][NSHIFT+i] 
	 - PZ[2*9+4][NSHIFT+i] 
	 + PZ[0*9+4][NSHIFT+i] 
	 + PZ[1*9+4][NSHIFT+i]/2
         - PZ[3*9+2][NSHIFT+i] 
	 + PZ[1*9+2][NSHIFT+i] 
	 + PZ[0*9+2][NSHIFT+i]*3 
	 - PZ[2*9+2][NSHIFT+i]/2) / 16;
    Aacdbef[12*NK+i] = pre_A * (Jterms + PZterms);

    Jterms = J[5*9+5][NSHIFT+i]/5 
      + J[3*9+5][NSHIFT+i]/2 
      + J[4*9+5][NSHIFT+i]/6 
      + 0.55*J[2*9+5][NSHIFT+i]
      + J[2*9+7][NSHIFT+i]/4 
      +J[0*9+5][NSHIFT+i]/4 
      + J[1*9+5][NSHIFT+i]/12;
    PZterms = 0;
    Aacdbef[13*NK+i] = pre_A * (Jterms + PZterms);

    Jterms = J[5*9+7][NSHIFT+i]/5 
      + J[3*9+7][NSHIFT+i]/2 
      + J[4*9+7][NSHIFT+i]/6 
      + 0.55*J[2*9+7][NSHIFT+i]
      + J[2*9+5][NSHIFT+i]/4 
      +J[0*9+7][NSHIFT+i]/4 
      + J[1*9+7][NSHIFT+i]/12;
    PZterms = -PZ[0*9+5][NSHIFT+i]/12.0
      + (PZ[4*9+7][NSHIFT+i] 
	 - PZ[2*9+7][NSHIFT+i] 
	 + PZ[0*9+7][NSHIFT+i] 
	 + PZ[1*9+7][NSHIFT+i]/2
         - PZ[3*9+5][NSHIFT+i] 
	 + PZ[1*9+5][NSHIFT+i] 
	 + PZ[0*9+5][NSHIFT+i]*3 
	 - PZ[2*9+5][NSHIFT+i]/2) / 16;
    Aacdbef[14*NK+i] = pre_A * (Jterms + PZterms);

    Jterms = J[5*9+8][NSHIFT+i]/5 
      + J[3*9+8][NSHIFT+i]/2 
      + J[4*9+8][NSHIFT+i]/6 
      + 0.55*J[2*9+8][NSHIFT+i]
      + J[2*9+8][NSHIFT+i]/4 
      +J[0*9+8][NSHIFT+i]/4 
      + J[1*9+8][NSHIFT+i]/12;
    PZterms = 0;
    Aacdbef[15*NK+i] = pre_A * (Jterms + PZterms);

    Jterms = (J[5*9+1][NSHIFT+i]/5 
	      + J[3*9+1][NSHIFT+i]/2 
	      + J[4*9+1][NSHIFT+i]/6 
	      + 0.55*J[2*9+1][NSHIFT+i]
              + J[2*9+3][NSHIFT+i]/4 
	      + J[0*9+1][NSHIFT+i]/4 
	      + J[1*9+1][NSHIFT+i]/12) * 2.0;
    PZterms = (-PZ[4*9+1][NSHIFT+i]*2 
	       + PZ[2*9+1][NSHIFT+i]*2 
	       - PZ[0*9+1][NSHIFT+i]*2 
	       - PZ[1*9+1][NSHIFT+i]
               + PZ[6*9+3][NSHIFT+i]*2 
	       - PZ[4*9+3][NSHIFT+i]*4 
	       + PZ[2*9+3][NSHIFT+i]) / 16.0;
    Aacdbef[56*NK+i] = pre_A * (Jterms + PZterms);

    Jterms = J[5*9+2][NSHIFT+i]/5 
      + J[3*9+2][NSHIFT+i]/2 
      + J[4*9+2][NSHIFT+i]/6 
      + 0.55*J[2*9+2][NSHIFT+i]
      + J[2*9+6][NSHIFT+i]/4 
      + J[0*9+2][NSHIFT+i]/4 
      + J[1*9+2][NSHIFT+i]/12
      + J[5*9+4][NSHIFT+i]/5 
      + J[3*9+4][NSHIFT+i]/2 
      + J[4*9+4][NSHIFT+i]/6 
      + 0.55*J[2*9+4][NSHIFT+i]
      + J[2*9+4][NSHIFT+i]/4 
      + J[0*9+4][NSHIFT+i]/4 
      + J[1*9+4][NSHIFT+i]/12;
    PZterms = (-PZ[4*9+4][NSHIFT+i] 
	       + PZ[2*9+4][NSHIFT+i] 
	       - PZ[0*9+4][NSHIFT+i] 
	       - PZ[1*9+4][NSHIFT+i]/2
               + PZ[6*9+6][NSHIFT+i] 
	       - PZ[4*9+6][NSHIFT+i]*2 
	       + PZ[2*9+6][NSHIFT+i]/2) / 16.0;
    Aacdbef[57*NK+i] = pre_A * (Jterms + PZterms);

    Jterms = (J[5*9+5][NSHIFT+i]/5 
	      + J[3*9+5][NSHIFT+i]/2 
	      + J[4*9+5][NSHIFT+i]/6 
	      + 0.55*J[2*9+5][NSHIFT+i]
              + J[2*9+7][NSHIFT+i]/4 
	      + J[0*9+5][NSHIFT+i]/4 
	      + J[1*9+5][NSHIFT+i]/12) * 2.0;
    PZterms = 0;
    Aacdbef[59*NK+i] = pre_A * (Jterms + PZterms);

    Jterms = J[6*9+4][NSHIFT+i]*8/35 
      + 0.4*J[5*9+4][NSHIFT+i] 
      + 0.4*J[5*9+4][NSHIFT+i]
      + J[3*9+4][NSHIFT+i]*19/21 
      + J[4*9+4][NSHIFT+i]/6 
      + J[4*9+4][NSHIFT+i]/6 
      + 0.6*J[2*9+4][NSHIFT+i]
      + 0.6*J[2*9+4][NSHIFT+i] 
      + J[0*9+4][NSHIFT+i]*11/30 
      + J[1*9+4][NSHIFT+i]/12 
      + J[1*9+4][NSHIFT+i]/12;
    PZterms = (-PZ[4*9+2][NSHIFT+i]*2 
	       + PZ[2*9+2][NSHIFT+i]*2 
	       - PZ[0*9+2][NSHIFT+i]*2 
	       - PZ[1*9+2][NSHIFT+i]
               + PZ[6*9+4][NSHIFT+i]*2 
	       - PZ[4*9+4][NSHIFT+i]*4 
	       + PZ[2*9+4][NSHIFT+i]) / 16.0;
    Aacdbef[60*NK+i] = pre_A * (Jterms + PZterms);

    Jterms = J[6*9+5][NSHIFT+i]*8/35 
      + 0.4*J[5*9+5][NSHIFT+i] 
      + 0.4*J[5*9+7][NSHIFT+i]
      + J[3*9+5][NSHIFT+i]*19/21 
      + J[4*9+5][NSHIFT+i]/6 
      + J[4*9+7][NSHIFT+i]/6 
      + 0.6*J[2*9+5][NSHIFT+i]
      + 0.6*J[2*9+7][NSHIFT+i] 
      + J[0*9+5][NSHIFT+i]*11/30 
      + J[1*9+5][NSHIFT+i]/12 
      + J[1*9+7][NSHIFT+i]/12;
    PZterms = (-PZ[4*9+5][NSHIFT+i] 
	       + PZ[2*9+5][NSHIFT+i] 
	       - PZ[0*9+5][NSHIFT+i] 
	       - PZ[1*9+5][NSHIFT+i]/2
               + PZ[6*9+7][NSHIFT+i] 
	       - PZ[4*9+7][NSHIFT+i]*2 
	       + PZ[2*9+7][NSHIFT+i]/2) / 16.0;
    Aacdbef[61*NK+i] = pre_A * (Jterms + PZterms);
    
    Jterms = J[6*9+8][NSHIFT+i]*8/35 
      + 0.4*J[5*9+8][NSHIFT+i] 
      + 0.4*J[5*9+8][NSHIFT+i]
      + J[3*9+8][NSHIFT+i]*19/21 
      + J[4*9+8][NSHIFT+i]/6 
      + J[4*9+8][NSHIFT+i]/6 
      + 0.6*J[2*9+8][NSHIFT+i]
      + 0.6*J[2*9+8][NSHIFT+i] 
      + J[0*9+8][NSHIFT+i]*11/30 
      + J[1*9+8][NSHIFT+i]/12 
      + J[1*9+8][NSHIFT+i]/12;
    PZterms = 0;
    Aacdbef[63*NK+i] = pre_A * (Jterms + PZterms);

    //symmetries: A_{acd,bef} = A_{adc,bfe}
    Aacdbef[16*NK+i] = Aacdbef[8*NK+i];
    Aacdbef[18*NK+i] = Aacdbef[9*NK+i];
    Aacdbef[17*NK+i] = Aacdbef[10*NK+i];
    Aacdbef[19*NK+i] = Aacdbef[11*NK+i];
    Aacdbef[20*NK+i] = Aacdbef[12*NK+i];
    Aacdbef[22*NK+i] = Aacdbef[13*NK+i];
    Aacdbef[21*NK+i] = Aacdbef[14*NK+i];
    Aacdbef[23*NK+i] = Aacdbef[15*NK+i];
    Aacdbef[58*NK+i] = Aacdbef[57*NK+i];
    Aacdbef[62*NK+i] = Aacdbef[61*NK+i];
  }

  return 0;
}

