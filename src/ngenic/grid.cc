/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  grid.cc
 *
 *  \brief routines for setting up an unperturbed particle load
 */

#include "gadgetconfig.h"

#if defined (CREATE_GRID) || defined (ADDITIONAL_GRID) || defined (CREATE_HDM_GRID)

#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../logs/timer.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../ngenic/ngenic.h"
#include "../system/system.h"

void ngenic::create_grid(void)
{
  long long gridSize    = All.GridSize;
  long long partTotal   = gridSize * gridSize * gridSize;
  long long partPerTask = partTotal / NTask;

  Sp->RegionLen     = All.BoxSize;
  Sp->FacCoordToInt = pow(2.0, BITS_FOR_POSITIONS) / Sp->RegionLen;
  Sp->FacIntToCoord = Sp->RegionLen / pow(2.0, BITS_FOR_POSITIONS);

  All.Time = All.TimeBegin;

  double masstot = All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G) * All.BoxSize * All.BoxSize * All.BoxSize;

  double m = masstot / (partTotal);

  for(int i = 0; i < NTYPES; i++)
    All.MassTable[i] = 0.;

  All.MassTable[1] = m;

  Sp->NumGas  = 0;
  Sp->NumPart = partPerTask;

  if(ThisTask == NTask - 1)
    {
      Sp->NumPart = partTotal - Sp->NumPart * (NTask - 1);
    }

  int max_load, max_sphload;
  MPI_Allreduce(&Sp->NumPart, &max_load, 1, MPI_INT, MPI_MAX, Communicator);
  MPI_Allreduce(&Sp->NumGas, &max_sphload, 1, MPI_INT, MPI_MAX, Communicator);

#ifdef GENERATE_GAS_IN_ICS
  Sp->TotNumGas  = partTotal;
  Sp->TotNumPart = 2 * partTotal;
  max_sphload    = max_load;
  max_load *= 2;
#else
  Sp->TotNumPart = partTotal;
  Sp->TotNumGas  = 0;
#endif

  Sp->MaxPart    = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
  Sp->MaxPartSph = max_sphload / (1.0 - 2 * ALLOC_TOLERANCE);

  Sp->allocate_memory();

  for(int i = 0; i < Sp->NumPart; i++)
    {
      long long ipcell = ThisTask * partPerTask + i;
      int x            = ipcell / (All.GridSize * All.GridSize);
      int xr           = ipcell % (All.GridSize * All.GridSize);
      int y            = xr / All.GridSize;
      int z            = xr % All.GridSize;

      double xyz[3];
      xyz[0] = x * All.BoxSize / All.GridSize;
      xyz[1] = y * All.BoxSize / All.GridSize;
      xyz[2] = z * All.BoxSize / All.GridSize;

      Sp->pos_to_intpos(xyz, Sp->P[i].IntPos);

      Sp->P[i].Vel[0] = 0.;
      Sp->P[i].Vel[1] = 0.;
      Sp->P[i].Vel[2] = 0.;

      Sp->P[i].ID.set(ipcell + 1);

      Sp->P[i].setType(1);
    }

  mpi_printf("NGENIC: generated grid of size %d\n", All.GridSize);
}

#ifdef ADDITIONAL_GRID
void ngenic::additional_grid(void)
{
  mpi_printf("NGENIC: additional grid called\n");

  long long gridSize    = All.GridSize;
  long long partTotal   = gridSize * gridSize * gridSize;
  long long partPerTask = partTotal / NTask;

  Sp->RegionLen     = All.BoxSize;
  Sp->FacCoordToInt = pow(2.0, BITS_FOR_POSITIONS) / Sp->RegionLen;
  Sp->FacIntToCoord = Sp->RegionLen / pow(2.0, BITS_FOR_POSITIONS);

  All.Time = All.TimeBegin;

  //double masstot = (All.OmegaNuPart / All.N_tau_part) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G) * All.BoxSize * All.BoxSize * All.BoxSize;

  //double m = masstot / (partTotal);

  All.MassTable[All.N_tau_part+1] = 0; // this needs to be set to 0, otherwise Gadget won't write the individual particle masses

  Sp->NumICPart = Sp->NumPart;

  Sp->NumGas   = 0;

  mpi_printf("MaxPart from IC = %d\n", Sp->MaxPart);

  int NewPart = 0;

  if(ThisTask != NTask - 1)
    {
      NewPart = partPerTask;
    }

  if(ThisTask == NTask - 1)
    {
      NewPart = partTotal - partPerTask * (NTask - 1);
    }

  printf("NewPart on Task %d is %d\n", ThisTask, NewPart);

  int max_load, max_sphload;
  MPI_Allreduce(&Sp->NumPart, &max_load, 1, MPI_INT, MPI_MAX, Communicator);
  MPI_Allreduce(&Sp->NumGas, &max_sphload, 1, MPI_INT, MPI_MAX, Communicator);

  int max_load_ic;
  MPI_Allreduce(&Sp->NumICPart, &max_load_ic, 1, MPI_INT, MPI_MAX, Communicator);

  Sp->TotNumPart += partTotal;
  Sp->TotNumGas  = 0;

  Sp->MaxPart    = max_load_ic / (1.0 - 2 * ALLOC_TOLERANCE);
  Sp->MaxPartSph = max_sphload / (1.0 - 2 * ALLOC_TOLERANCE);

  Sp->NumPart += NewPart;

  int max_load_new;
  MPI_Allreduce(&Sp->NumPart, &max_load_new, 1, MPI_INT, MPI_MAX, Communicator);

  int MaxPartNew = max_load_new / (1.0 - 2 * ALLOC_TOLERANCE);

  mpi_printf("MaxPart with new particles = %d\n", MaxPartNew);

  mpi_printf("TotNumPart = %lld, TotNumICPart = %lld, NumPart = %d, NumICPart = %d\n", Sp->TotNumPart, Sp->TotNumICPart, Sp->NumPart, Sp->NumICPart);

/*
  for(int i=0; i<10; i++)
    {
      Sp->print_particle_info(i);
    }
*/
  Sp->reallocate_memory_maxpart(MaxPartNew);
/*
  for(int i=0; i<10; i++)
    {
      Sp->print_particle_info(i);
    }
*/

  long long IDMaxCount = 0;
  int PartTag = 0;

  for(int i = 0; i < Sp->NumICPart; i++)
    {
      if(Sp->P[i].ID.get() > IDMaxCount) 
        {
          IDMaxCount = Sp->P[i].ID.get();
          PartTag = i;
        }
    }

  MPI_Barrier(Communicator);

  long long IDMaxGlobal = 0;
  MPI_Allreduce(&IDMaxCount, &IDMaxGlobal, 1, MPI_LONG_LONG_INT, MPI_MAX, Communicator);

  mpi_printf("Max ID from IC = %lld\n", IDMaxGlobal);

  for(int i = Sp->NumICPart; i < Sp->NumPart; i++)
    {
      long long ipcell = ThisTask * partPerTask + i - Sp->NumICPart;
      int x            = ipcell / (All.GridSize * All.GridSize);
      int xr           = ipcell % (All.GridSize * All.GridSize);
      int y            = xr / All.GridSize;
      int z            = xr % All.GridSize;

      double xyz[3];
      xyz[0] = x * All.BoxSize / All.GridSize;
      xyz[1] = y * All.BoxSize / All.GridSize;
      xyz[2] = z * All.BoxSize / All.GridSize;

      Sp->pos_to_intpos(xyz, Sp->P[i].IntPos);

      Sp->P[i].Vel[0] = 0.;
      Sp->P[i].Vel[1] = 0.;
      Sp->P[i].Vel[2] = 0.;

      Sp->P[i].ID.set(ipcell + 1 + IDMaxGlobal); // the IDMaxGlobal offsets the ID of new particles to avoid doubling with the IC particles existing ID.

      Sp->P[i].setType(All.N_tau_part + 1);
    }
}
#endif

#ifdef CREATE_HDM_GRID
void ngenic::create_hdm_grid(void)
{
  mpi_printf("NGENIC: create hdm grid called\n");
  
  long long gridSize    = All.GridSize;
  long long partTotal   = gridSize * gridSize * gridSize;
  long long partPerTask = partTotal / NTask;

  Sp->RegionLen     = All.BoxSize;
  Sp->FacCoordToInt = pow(2.0, BITS_FOR_POSITIONS) / Sp->RegionLen;
  Sp->FacIntToCoord = Sp->RegionLen / pow(2.0, BITS_FOR_POSITIONS);

  All.Time = All.TimeBegin;

  //double masstot = (All.OmegaNuPart / All.N_tau_part) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G) * All.BoxSize * All.BoxSize * All.BoxSize;
  //double m = masstot / (partTotal);

  /* Set the mass array to zero to let gagdet write particle masses */
  for (int ii=0;ii<All.N_hdm_types;ii++) 
	  All.MassTable[ii+2] = 0.;

  Sp->NumICPart = Sp->NumPart;

  Sp->NumGas   = 0;

  mpi_printf("MaxPart from IC = %d\n", Sp->MaxPart);

  int NewPart = 0;

  if(ThisTask != NTask - 1)
    {









}
#endif



#endif
