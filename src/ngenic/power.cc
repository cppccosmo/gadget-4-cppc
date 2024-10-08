/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  power.cc
 *
 *  \brief auxiliary routines for computing the linear power spectrum for the ICs
 */

#include "gadgetconfig.h"

#ifdef NGENIC

#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/mymalloc.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../ngenic/ngenic.h"
#include "../pm/pm_mpi_fft.h"
#include "../system/system.h"
#include "../neutrinos/neutrino.h"

double ngenic::ngenic_power_spec(double k)
{
  double power = 0;

#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
  if(Type == 1)
#endif
    switch(All.PowerSpectrumType)
      {
        case 1:
          power = ngenic_powerspec_eh(k);
          break;

        case 2:
          power = ngenic_powerspec_tabulated(k);
          break;

        default:
          power = ngenic_powerspec_efstathiou(k);
          break;
      }

#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
  if(Type == 2)
    {
      power = PowerSpec_DM_2ndSpecies(k);
    }
#endif

  return power;
}

#ifdef CB_PHASE
double ngenic::ngenic_CB_power_spec(double k) 
{
  double power = 0;

  power = ngenic_CB_powerspec_tabulated(k);

  return power;
}
#endif

double ngenic::ngenic_growth_rate(double k)
{
  double f_growth = 0;

  f_growth = ngenic_f_growth_tabulated(k);

  return f_growth;
}

void ngenic::free_power_table(void) { Mem.myfree(PowerTable); }

#ifdef CB_PHASE
void ngenic::free_CB_power_table(void) { Mem.myfree(CBPowerTable); }
#endif

void ngenic::free_f_growth_table(void) { Mem.myfree(GrowthTable); }

void ngenic::read_power_table(void)
{
  FILE *fd;
  char buf[MAXLEN_PATH];
  double k, p;

  sprintf(buf, All.PowerSpectrumFile);

  if(!(fd = fopen(buf, "r")))
    {
      Terminate("can't read input spectrum in file '%s' on task %d\n", buf, ThisTask);
    }

  NPowerTable = 0;
  do
    {
      if(fscanf(fd, " %lg %lg ", &k, &p) == 2)
        NPowerTable++;
      else
        break;
    }
  while(1);

  fclose(fd);

  mpi_printf("NGENIC_POWER: found %d rows in input spectrum table\n", NPowerTable);

  PowerTable = (pow_table *)Mem.mymalloc("PowerTable", NPowerTable * sizeof(pow_table));

  sprintf(buf, All.PowerSpectrumFile);

  if(!(fd = fopen(buf, "r")))
    {
      Terminate("can't read input spectrum in file '%s' on task %d\n", buf, ThisTask);
    }

  NPowerTable = 0;
  do
    {
      double p;

      if(fscanf(fd, " %lg %lg ", &k, &p) == 2)
        {
          PowerTable[NPowerTable].k = k;
          PowerTable[NPowerTable].P = p;
          NPowerTable++;
        }
      else
        break;
    }
  while(1);

  fclose(fd);

  std::sort(PowerTable, PowerTable + NPowerTable);
}

#ifdef CB_PHASE
void ngenic::read_CB_power_table(void)
{
  FILE *fdcb;
  char buf[MAXLEN_PATH];
  double k, p;

  sprintf(buf, All.CBPowerSpectrumFile);

  if(!(fdcb = fopen(buf, "r")))
    {
      Terminate("can't read input CB spectrum in file '%s' on task %d\n", buf, ThisTask);
    }

  NCBPowerTable = 0;
  do
    {
      if(fscanf(fdcb, " %lg %lg ", &k, &p) == 2)
        NCBPowerTable++;
      else
        break;
    }
  while(1);

  fclose(fdcb);

  mpi_printf("NGENIC_POWER: found %d rows in input CB spectrum table\n", NCBPowerTable);

  CBPowerTable = (pow_table *)Mem.mymalloc("CBPowerTable", NCBPowerTable * sizeof(pow_table));

  sprintf(buf, All.CBPowerSpectrumFile);

  if(!(fdcb = fopen(buf, "r")))
    {
      Terminate("can't read input CB spectrum in file '%s' on task %d\n", buf, ThisTask);
    }

  NCBPowerTable = 0;
  do
    {
      double p;

      if(fscanf(fdcb, " %lg %lg ", &k, &p) == 2)
        {
          CBPowerTable[NCBPowerTable].k = k;
          CBPowerTable[NCBPowerTable].P = p;
          NCBPowerTable++;
        }
      else
        break;
    }
  while(1);

  fclose(fdcb);
  
  std::sort(CBPowerTable, CBPowerTable + NCBPowerTable);   
}
#endif

void ngenic::read_f_growth_table(void)
{
  FILE *fd;
  char buf[MAXLEN_PATH];
  double k, f;

  sprintf(buf, All.GrowthRateFile);

  if(!(fd = fopen(buf, "r")))
    {
      Terminate("can't read input growth rate in file '%s' on task %d\n", buf, ThisTask);
    }

  NGrowthTable = 0;
  do 
    {
      if(fscanf(fd, " %lg %lg ", &k, &f) == 2)
        NGrowthTable++;
      else
        break;
    }
  while(1);

  fclose(fd);

  mpi_printf("NGENIC_POWER: found %d rows in growth rate table\n", NGrowthTable);

  GrowthTable = (growth_table *)Mem.mymalloc("GrowthTable", NGrowthTable * sizeof(growth_table));
  
  sprintf(buf, All.GrowthRateFile);

  if(!(fd = fopen(buf, "r")))
    {
      Terminate("can't read input growth rate in file '%s' on task %d\n", buf, ThisTask);
    }

  NGrowthTable = 0;
  do
    {
      double f;
      
      if(fscanf(fd, " %lg %lg ", &k, &f) == 2)
        {
          GrowthTable[NGrowthTable].k = k;
          GrowthTable[NGrowthTable].f = f;
          NGrowthTable++;
        }
      else
        break;
    }
  while(1);
  
  fclose(fd);

  std::sort(GrowthTable, GrowthTable + NGrowthTable);
}

void ngenic::ngenic_initialize_powerspectrum(void)
{
  double res;

  AA = 6.4 / All.ShapeGamma * (3.085678e24 / All.UnitLength_in_cm);
  BB = 3.0 / All.ShapeGamma * (3.085678e24 / All.UnitLength_in_cm);
  CC = 1.7 / All.ShapeGamma * (3.085678e24 / All.UnitLength_in_cm);
  nu = 1.13;

  R8 = 8 * (3.085678e24 / All.UnitLength_in_cm); /* 8 Mpc/h */

  if(All.PowerSpectrumType == 2)
    read_power_table();

  if(All.ReNormalizeInputSpectrum == 0 && All.PowerSpectrumType == 2)
    {
      Norm = 1.0;
      /* tabulated file is already at the initial redshift */
      Dplus = 1.0;
    }
  else
    {
#ifdef DIFFERENT_TRANSFER_FUNC
      Type = 1;
#endif
      Norm = 1.0;
      res  = ngenic_tophat_sigma2(R8);

      if(ThisTask == 0 && All.PowerSpectrumType == 2)
        printf("\nNormalization of spectrum in file:  Sigma8 = %g\n", sqrt(res));

      Norm = All.Sigma8 * All.Sigma8 / res;

      if(ThisTask == 0 && All.PowerSpectrumType == 2)
        printf("Normalization adjusted to  Sigma8=%g   (Normfac=%g)\n\n", All.Sigma8, Norm);

      Dplus = ngenic_growth_factor(All.cf_atime, 1.0);
    }
  mpi_printf("NGENIC: Dplus=%g\n", Dplus);
}

#ifdef CB_PHASE
void ngenic::ngenic_initialize_CB_powerspectrum(void)
{
  read_CB_power_table();
}
#endif

void ngenic::ngenic_initialize_growth_rate(void)
{
  read_f_growth_table();
}

double ngenic::ngenic_powerspec_tabulated(double k)
{
  double kold = k;

  k *= (All.InputSpectrum_UnitLength_in_cm / All.UnitLength_in_cm);  // convert to h/Mpc

  if(k < PowerTable[0].k || k > PowerTable[NPowerTable - 1].k)
    return 0;

  int binlow  = 0;
  int binhigh = NPowerTable - 1;

  while(binhigh - binlow > 1)
    {
      int binmid = (binhigh + binlow) / 2;
      if(k < PowerTable[binmid].k)
        binhigh = binmid;
      else
        binlow = binmid;
    }

  double dk = PowerTable[binhigh].k - PowerTable[binlow].k;

  if(dk == 0)
    Terminate("dk == 0");

  double u = (k - PowerTable[binlow].k) / dk;

  double P = (1 - u) * PowerTable[binlow].P + u * PowerTable[binhigh].P;

  return P;
}

#ifdef CB_PHASE
double ngenic::ngenic_CB_powerspec_tabulated(double k)
{
  double kold = k;

  k *= (All.InputSpectrum_UnitLength_in_cm / All.UnitLength_in_cm);  // convert to h/Mpc

  if(k < CBPowerTable[0].k || k > CBPowerTable[NCBPowerTable - 1].k)
    return 0;

  int binlow = 0;
  int binhigh = NCBPowerTable - 1;

  while(binhigh - binlow > 1)
    {
      int binmid = (binhigh + binlow) / 2;
      if(k < CBPowerTable[binmid].k)
        binhigh = binmid;
      else
        binlow = binmid;
    }

  double dk = CBPowerTable[binhigh].k - CBPowerTable[binlow].k;

  if(dk == 0)
    Terminate("dk == 0");

  double u = (k - CBPowerTable[binlow].k) / dk;

  double P = (1 - u) * CBPowerTable[binlow].P + u * CBPowerTable[binhigh].P;

  return P;
}
#endif

double ngenic::ngenic_f_growth_tabulated(double k)
{
  k *= (All.InputSpectrum_UnitLength_in_cm / All.UnitLength_in_cm);  // convert to h/Mpc

  if(k < GrowthTable[0].k || k > GrowthTable[NGrowthTable - 1].k)
    return 0;

  int binlow = 0;
  int binhigh = NGrowthTable - 1;

  while(binhigh - binlow > 1)
    {
      int binmid = (binhigh + binlow) / 2;
      if(k < GrowthTable[binmid].k)
        binhigh = binmid;
      else
        binlow = binmid;
    }

  double dk = GrowthTable[binhigh].k - GrowthTable[binlow].k;

  if(dk == 0)
    Terminate("dk == 0");

  double u = (k - GrowthTable[binlow].k)/dk;
  
  double f = (1 - u) * GrowthTable[binlow].f + u * GrowthTable[binhigh].f;

  return f;
}

double ngenic::ngenic_powerspec_efstathiou(double k)
{
  return Norm * k / pow(1 + pow(AA * k + pow(BB * k, 1.5) + CC * CC * k * k, nu), 2 / nu);
}

double ngenic::ngenic_powerspec_eh(double k) /* Eisenstein & Hu */ { return Norm * k * pow(ngenic_tk_eh(k), 2); }

double ngenic::ngenic_tk_eh(double k) /* from Martin White */
{
  double q, theta, ommh2, a, s, gamma, L0, C0;
  double tmp;
  double omegam, ombh2;

  /* other input parameters */

  omegam = All.Omega0;
  ombh2  = All.OmegaBaryon * All.HubbleParam * All.HubbleParam;

  if(All.OmegaBaryon == 0)
    ombh2 = 0.04 * All.HubbleParam * All.HubbleParam;

  k *= (3.085678e24 / All.UnitLength_in_cm); /* convert to h/Mpc */

  theta = 2.728 / 2.7;
  ommh2 = omegam * All.HubbleParam * All.HubbleParam;
  s     = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * All.HubbleParam;
  a     = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2 + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
  gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
  gamma *= omegam * All.HubbleParam;
  q   = k * theta * theta / gamma;
  L0  = log(2. * exp(1.) + 1.8 * q);
  C0  = 14.2 + 731. / (1. + 62.5 * q);
  tmp = L0 / (L0 + C0 * q * q);
  return (tmp);
}

double ngenic::ngenic_tophat_sigma2(double R)
{
  const int worksize = 1000000;

  double result, abserr, kmin, kmax;
  gsl_function F;

  myparams par = {R, this};

  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(worksize);
  F.function                           = &sigma2_int;
  F.params                             = &par;

  if(All.PowerSpectrumType == 2)
    {
      kmin = pow(10.0, PowerTable[0].k) * (All.UnitLength_in_cm / All.InputSpectrum_UnitLength_in_cm);
      kmax = pow(10.0, PowerTable[NPowerTable - 1].k) * (All.UnitLength_in_cm / All.InputSpectrum_UnitLength_in_cm);
    }
  else
    {
      kmin = 1.0e-15 / R;
      kmax = 1.0e3 / R;
    }

  if(All.PowerSpectrumType == 2)
    {
      /* because of the oscillatory behaviour of the integrand, the gsl_integration_qag() has trouble with its error estimates
       * when the function is piece-wise interpolated. That's why we integrate the tabulated function segment by segment.
       */

      /* first get a rough result with up to 10% relative error */
      gsl_integration_qag(&F, log(kmin), log(kmax), 0, 0.1, worksize, GSL_INTEG_GAUSS15, workspace, &result, &abserr);

      /* now set a low absolute error bound for each segment */
      double errbound = 1.0e-8 / NPowerTable * result;
      result          = 0.0;

      for(int i = 0; i < NPowerTable - 2; i++)
        {
          double k0 = pow(10.0, PowerTable[i].k) * (All.UnitLength_in_cm / All.InputSpectrum_UnitLength_in_cm);
          double k1 = pow(10.0, PowerTable[i + 1].k) * (All.UnitLength_in_cm / All.InputSpectrum_UnitLength_in_cm);
          double x;

          gsl_integration_qag(&F, log(k0), log(k1), errbound, 0, worksize, GSL_INTEG_GAUSS15, workspace, &x, &abserr);

          result += x;
        }
    }
  else
    {
      /* for the smooth analytic function, we integrate directly with a relative error estimate */
      gsl_integration_qag(&F, log(kmin), log(kmax), 0, 1.0e-8, worksize, GSL_INTEG_GAUSS15, workspace, &result, &abserr);
    }

  gsl_integration_workspace_free(workspace);

  return result;
}

double ngenic::ngenic_growth_factor(double astart, double aend) { return ngenic_growth(aend) / ngenic_growth(astart); }

double ngenic::ngenic_growth(double a)
{
  double hubble_a;

  hubble_a = sqrt(All.Omega0 / (a * a * a) + (1 - All.Omega0 - All.OmegaLambda) / (a * a) + All.OmegaLambda);

  const int worksize = 100000;

  double result, abserr;
  gsl_function F;

  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(worksize);
  F.function                           = &ngenic_growth_int;

  gsl_integration_qag(&F, 0, a, 0, 1.0e-8, worksize, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return hubble_a * result;
}

double ngenic::ngenic_f1_omega(double a)
{
  double omega_a;

  omega_a = All.Omega0 / (All.Omega0 + a * (1 - All.Omega0 - All.OmegaLambda) + a * a * a * All.OmegaLambda);

  return pow(omega_a, 0.6);
}

double ngenic::ngenic_f2_omega(double a)
{
  double omega_a;

  omega_a = All.Omega0 / (All.Omega0 + a * (1 - All.Omega0 - All.OmegaLambda) + a * a * a * All.OmegaLambda);

  return 2 * pow(omega_a, 4.0 / 7);
}

#endif
