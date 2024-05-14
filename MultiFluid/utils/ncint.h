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

/*******************************************************************************
composite Newton-Cotes integration
Use closed 4th-degree method in steps until fewer than 4 points remain,
then switch to a lower-degree method.

 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double ncint_boole_coeffs[5] =
  {0.311111111111111, 1.42222222222222, 0.533333333333333, 1.42222222222222,
   0.311111111111111};

double ncint_boole_step(double dx, const double *f){
  double res = 0;
  for(int i=0; i<5; i++) res += ncint_boole_coeffs[i] * f[i];
  return res * dx;
}

const double ncint_simpson38_coeffs[4] = {0.375, 1.125, 1.125, 0.375};

double ncint_simpson38_step(double dx, const double *f){
  double res = 0;
  for(int i=0; i<4; i++) res += ncint_simpson38_coeffs[i] * f[i];
  return res * dx;
}

const double ncint_simpson_coeffs[3] =
  {0.333333333333333, 1.33333333333333, 0.333333333333333};

double ncint_simpson_step(double dx, const double *f){
  double res = 0;
  for(int i=0; i<3; i++) res += ncint_simpson_coeffs[i] * f[i];
  return res * dx;
}

const double ncint_trapezoid_coeffs[2] = {0.5, 0.5};

double ncint_trapezoid_step(double dx, const double *f){
  double res = 0;
  for(int i=0; i<2; i++) res += ncint_trapezoid_coeffs[i] * f[i];
  return res * dx;
}

//composite fixed-step-size Newton-Cotes integrator
double ncint_cf(int N, double dx, const double *f){
  if(N<2) return 0;
  int Nint = N-1; //number of intervals
  double result = 0;
  
  while(Nint >= 6){
    Nint -= 4;
    result += ncint_boole_step(dx, &f[Nint]);
  }

  switch(Nint){
  case 5:
    result += ncint_simpson38_step(dx,&f[2]) + ncint_simpson_step(dx,f);
    break;
  case 4:
    result += ncint_boole_step(dx,f);
    break;
  case 3:
    result += ncint_simpson38_step(dx,f);
    break;
  case 2:
    result += ncint_simpson_step(dx,f);
    break;
  case 1: //should only get here if N=2, otherwise don't use trapezoid rule
    result += ncint_trapezoid_step(dx,f);
    break;
    //case 0: //shouldn't get here
    //break;
  default:
    printf("ERROR in ncint: Nint=%i invalid. Quitting.\n",Nint);
    fflush(stdout);
    abort();
    break;
  }

  return result;   
}

//composite trapezoid-rule integrator for any step sizes
double ncint_composite_trapezoid(int N, const double *x, const double *f){
  double res = 0;
  for(int i=0; i<N-1; i++) res += ncint_trapezoid_step(x[i+1]-x[i],f+i);
  return res;
}

/*****************************************************************************
//variable step size integration: needs GSL to solve linear systems,
//make sure to include it in main program

//single step with variable step sizes 
double ncint_variable_step(int N, const double *x, const double *f){

  //find weights by solving NxN linear system
  gsl_matrix *Xji = gsl_matrix_alloc(N,N), *XLU = gsl_matrix_alloc(N,N);
  double x0jp1[N], xNm1jp1[N]; //x_0^{j+1} and x_{N-1}^{j+1}
  for(int i=0; i<N; i++){
    double xj = 1; //x_i^j
    for(int j=0; j<N; j++){
      gsl_matrix_set(Xji,j,i,xj);
      gsl_matrix_set(XLU,j,i,xj);
      xj *= x[i];
      if(i==0) x0jp1[j] = xj;
      else if(i==N-1) xNm1jp1[j] = xj;
    }
  }

  gsl_vector *Yj = gsl_vector_alloc(N), *wi = gsl_vector_alloc(N),
    *work = gsl_vector_alloc(N);
  for(int j=0; j<N; j++) gsl_vector_set(Yj,j,(xNm1jp1[j]-x0jp1[j])/(1.0+j));

  int s;
  gsl_permutation *p = gsl_permutation_alloc(N);
  gsl_linalg_LU_decomp(XLU,p,&s);
  gsl_linalg_LU_solve(XLU,p,Yj,wi);
  gsl_linalg_LU_refine(Xji,XLU,p,Yj,wi,work);
  
  //compute result, clean up, and quit
  double integral = 0;
  for(int i=0; i<N; i++) integral += gsl_vector_get(wi,i) * f[i];
  
  gsl_matrix_free(Xji);
  gsl_matrix_free(XLU);
  gsl_vector_free(Yj);
  gsl_vector_free(wi);
  gsl_vector_free(work);
  gsl_permutation_free(p);

  return integral;
}

//composite integrator
//Strategy: take steps of size N_composite_maxstep until we approach the end
//of the array, then divide up steps to avoid using the lowest-order method.

const int N_composite_maxstep = 5;

double ncint_vf(int N, const double *x, const double *f){

  //initialize and take care of simplest cases
  if(N<2) return 0;
  else if(N==2) return ncint_trapezoid_step(x[1]-x[0], f);
  else if(N<=N_composite_maxstep) return ncint_variable_step(N, x, f);
  
  int Nint = N-1; //number of intervals
  int DNmax = N_composite_maxstep - 1; 
  double result = 0;

  //integrate in chunks of DNmax
  while(Nint >= DNmax+2){
    Nint -= DNmax;
    result += ncint_variable_step(N_composite_maxstep,&x[Nint], &f[Nint]);
  }

  //split the remainder into two smaller steps if needed
  if(Nint <= DNmax) result += ncint_variable_step(Nint+1, x, f);
  else{
    int NL = Nint/2, NR = Nint-NL;
    result += ncint_variable_step(NL+1, x, f)
      + ncint_variable_step(NR+1, &x[NL], &f[NL]);
  }
    
  return result;
}
*****************************************************************************/
