/*******************************************************************************
 * \copyright   This file is part of the GADGET4 N-body/SPH code developed
 * \copyright   by Volker Springel. Copyright (C) 2014-2020 by Volker Springel
 * \copyright   (vspringel@mpa-garching.mpg.de) and all contributing authors.
 *******************************************************************************/

/*! \file  ngenic.cc
 *
 *  \brief sets up cosmological initial conditions
 */

#include "gadgetconfig.h"

#ifdef NGENIC

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <algorithm>

#include "../data/allvars.h"
#include "../data/dtypes.h"
#include "../data/intposconvert.h"
#include "../data/mymalloc.h"
#include "../logs/logs.h"
#include "../main/simulation.h"
#include "../mpi_utils/mpi_utils.h"
#include "../mpi_utils/shared_mem_handler.h"
#include "../ngenic/ngenic.h"
#include "../pm/pm_mpi_fft.h"
#include "../system/system.h"

#ifdef GRIDX
#undef GRIDX
#undef GRIDY
#undef GRIDZ
#undef INTCELL
#endif

#define GRIDX (NGENIC)
#define GRIDY (NGENIC)
#define GRIDZ (NGENIC)

#define INTCELL ((~((MyIntPosType)0)) / GRIDX + 1)

#define GRIDz (GRIDZ / 2 + 1)
#define GRID2 (2 * GRIDz)

#define FI(x, y, z) (((large_array_offset)GRID2) * (GRIDY * (x) + (y)) + (z))
#define FC(c, z) (((large_array_offset)GRID2) * ((c)-myplan.firstcol_XY) + (z))

#if(GRIDZ > 1024)
typedef long long large_array_offset; /* use a larger data type in this case so that we can always address all cells of the 3D grid
                                         with a single index */
#else
typedef unsigned int large_array_offset;
#endif

#ifdef NUMPART_PER_TASK_LARGE
typedef long long large_numpart_type; /* if there is a risk that the local particle number times 8 overflows a 32-bit integer, this
                                         data type should be used */
#else
typedef int large_numpart_type;
#endif

void ngenic::ngenic_displace_additional_particles(void)
{
}

void ngenic::ngenic_displace_particles(void)
{
  TIMER_START(CPU_NGENIC);

  mpi_printf("NGENIC: computing displacement fields...\n");

  All.set_cosmo_factors_for_current_time();

  double vel_prefac1 = All.cf_atime * All.cf_hubble_a * ngenic_f1_omega(All.cf_atime);
  double vel_prefac2 = All.cf_atime * All.cf_hubble_a * ngenic_f2_omega(All.cf_atime);

  mpi_printf("f_of_k = %g\n", ngenic_f1_omega(All.cf_atime));

  if(All.NLR > 0) {
    mpi_printf("NGENIC: scale-dependent growth used.\n");
    vel_prefac1 /= ngenic_f1_omega(All.cf_atime);
  } else {
    vel_prefac1 /= ngenic_f1_omega(All.cf_atime);
  }

  vel_prefac1 /= sqrt(All.cf_atime); /* converts to Gadget velocity */
  vel_prefac2 /= sqrt(All.cf_atime); /* converts to Gadget velocity */

  mpi_printf("NGENIC: vel_prefac1= %g  hubble_a=%g   fom1=%g\n", vel_prefac1, All.cf_hubble_a, ngenic_f1_omega(All.cf_atime));
  mpi_printf("NGENIC: vel_prefac2= %g  hubble_a=%g   fom2=%g\n", vel_prefac2, All.cf_hubble_a, ngenic_f2_omega(All.cf_atime));

  rnd_generator_conjugate = gsl_rng_alloc(gsl_rng_ranlxd1);
  rnd_generator           = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rnd_generator, All.NgenicSeed);

  ngenic_initialize_powerspectrum();
#ifdef CB_PHASE
  ngenic_initialize_CB_powerspectrum();
#endif
  ngenic_initialize_growth_rate();

  ngenic_initialize_ffts();

  if(!(seedtable = (unsigned int *)Mem.mymalloc("seedtable", NGENIC * NGENIC * sizeof(unsigned int))))
    Terminate("could not allocate seed table");

  for(int i = 0; i < NGENIC / 2; i++)
    {
      for(int j = 0; j < i; j++)
        seedtable[i * NGENIC + j] = 0x7fffffff * gsl_rng_uniform(rnd_generator);

      for(int j = 0; j < i + 1; j++)
        seedtable[j * NGENIC + i] = 0x7fffffff * gsl_rng_uniform(rnd_generator);

      for(int j = 0; j < i; j++)
        seedtable[(NGENIC - 1 - i) * NGENIC + j] = 0x7fffffff * gsl_rng_uniform(rnd_generator);

      for(int j = 0; j < i + 1; j++)
        seedtable[(NGENIC - 1 - j) * NGENIC + i] = 0x7fffffff * gsl_rng_uniform(rnd_generator);

      for(int j = 0; j < i; j++)
        seedtable[i * NGENIC + (NGENIC - 1 - j)] = 0x7fffffff * gsl_rng_uniform(rnd_generator);

      for(int j = 0; j < i + 1; j++)
        seedtable[j * NGENIC + (NGENIC - 1 - i)] = 0x7fffffff * gsl_rng_uniform(rnd_generator);

      for(int j = 0; j < i; j++)
        seedtable[(NGENIC - 1 - i) * NGENIC + (NGENIC - 1 - j)] = 0x7fffffff * gsl_rng_uniform(rnd_generator);

      for(int j = 0; j < i + 1; j++)
        seedtable[(NGENIC - 1 - j) * NGENIC + (NGENIC - 1 - i)] = 0x7fffffff * gsl_rng_uniform(rnd_generator);
    }

  if(Shmem.Island_NTask != Shmem.World_NTask)
    {
      // We actually have multiple shared memory nodes in which we set aside one MPI rank for shared memory communictiona.
      // In this casem move the seedtable to the communication rank in order to consume this memory only once on the node

      if(Shmem.Island_ThisTask == 0)
        {
          size_t tab_len = NGENIC * NGENIC * sizeof(unsigned int);

          MPI_Send(&tab_len, sizeof(tab_len), MPI_BYTE, Shmem.MyShmRankInGlobal, TAG_TABLE_ALLOC, MPI_COMM_WORLD);
          MPI_Send(seedtable, tab_len, MPI_BYTE, Shmem.MyShmRankInGlobal, TAG_DMOM, MPI_COMM_WORLD);
        }

      Mem.myfree(seedtable);

      ptrdiff_t off;
      MPI_Bcast(&off, sizeof(ptrdiff_t), MPI_BYTE, Shmem.Island_NTask - 1, Shmem.SharedMemComm);

      seedtable = (unsigned int *)((char *)Shmem.SharedMemBaseAddr[Shmem.Island_NTask - 1] + off);
    }

  ngenic_distribute_particles();

#ifdef CB_PHASE
  /* get the phase information from the snapshot */
  FILE *fdgrid;
  char bufgrid[MAXLEN_PATH];

  mpi_printf("NGENIC: reading in the 3D delta file from the restart package\n");

  sprintf(bufgrid, All.deltagridFILE);

  if(!(fdgrid = fopen(bufgrid, "r"))) 
    {
      Terminate("can't read file\n");
    }

  /* read from the deltagrid file the first long int element, it should tell you how much data is stored in the file */
  long int maxfftglobalsize;
  fread(&maxfftglobalsize, sizeof(long int), 1, fdgrid);

  maxfftglobalsize += 2;
  /* now try figure out how much data each process needs to read from the file */
  /* the idea is that nslab_y / PMGRID tells you how big of a slab this particular process is responsible for */
  /* each process is only reading the phase info of its local slab */
  long int local_size_to_read = maxfftglobalsize * myplan.nslab_y / PMGRID;

  /* now initialise a real fftw array to hold the phase info, it will be fourier transformed into k-space */
  fft_real *deltagrid_local;
  deltagrid_local = (fft_real *)Mem.mymalloc_movable_clear(&deltagrid_local, "deltagrid_local", local_size_to_read * sizeof(fft_real));

  //printf("local_size_to_read = %ld on Task %d\n", local_size_to_read, ThisTask);
  /* set the location in the file to start reading */
  /* this indexes from the slabstart_y to ensure the correct spot is found even when the slab division is not equal */
  long int local_start_point = maxfftglobalsize * myplan.slabstart_y / PMGRID; 

  //mpi_printf("maxfftglobalsize = %ld\n", maxfftglobalsize);

  /* fseek sets the file pointer to the spot in the file where the array will start reading from */
  /* important to offset by one long int that is for storing the maxfftglobalsize */
  fseek(fdgrid, local_start_point * sizeof(fft_real) + sizeof(long int), SEEK_SET);

  /* fread fetches the local size amount of info from the file and stores it into the deltagrid_local array */
  fread(deltagrid_local, sizeof(fft_real), local_size_to_read, fdgrid);
  
  fclose(fdgrid);

  MPI_Barrier(Communicator);
#endif

  /* allocate displacement vectors */
#ifndef ADDITIONAL_GRID
  Pdisp = (disp_data *)Mem.mymalloc_clear("disp_data", Sp->NumPart * sizeof(disp_data));
#else
  int PdispSize = Sp->NumPart - Sp->NumICPart;
  Pdisp = (disp_data *)Mem.mymalloc_clear("disp_data", PdispSize * sizeof(disp_data));
#endif

  /* allocate velocity vectors */
#ifndef ADDITIONAL_GRID
  Pvel = (disp_data *)Mem.mymalloc_clear("vel_data", Sp->NumPart * sizeof(disp_data));
#else
  int PvelSize = Sp->NumPart - Sp->NumICPart;
  Pvel = (disp_data *)Mem.mymalloc_clear("vel_data", PvelSize * sizeof(disp_data));
#endif

  /* allocate thermal velocity vector */
#ifdef THERMAL_VEL_IC
  mpi_printf("NGENIC: generating random thermal velocity for IC\n");

  int vel_vec_size = Sp->NumPart - Sp->NumICPart;
  thermal_vel_vec = (double *)Mem.mymalloc_clear("thermal_vel_vec", vel_vec_size * 3 * sizeof(thermal_vel_vec)); 

  double stream_vel_phys = All.stream_vel / All.cf_atime; // the physical velocity is defined as the tau_i / (a * m_nu), the input stream_vel is only tau_i / m_nu 
  
  mpi_printf("NGENIC: thermal velocity of the stream at a = %g is %g km / s\n", All.cf_atime, stream_vel_phys);

  // populate the vector with randomly drawn 3D velocities
  for(int i = 0; i<vel_vec_size; i++)
    {
      double rand_x = gsl_ran_gaussian(rnd_generator, 1.);
      double rand_y = gsl_ran_gaussian(rnd_generator, 1.);
      double rand_z = gsl_ran_gaussian(rnd_generator, 1.);

      double norm_fac = sqrt(rand_x * rand_x + rand_y * rand_y + rand_z * rand_z);

      rand_x /= norm_fac;
      rand_y /= norm_fac;
      rand_z /= norm_fac;

      thermal_vel_vec[3*i+0] = stream_vel_phys * rand_x / sqrt(All.cf_atime);
      thermal_vel_vec[3*i+1] = stream_vel_phys * rand_y / sqrt(All.cf_atime);
      thermal_vel_vec[3*i+2] = stream_vel_phys * rand_z / sqrt(All.cf_atime); 
    }
#endif

#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
  for(Type = MinType; Type <= MaxType; Type++)
#endif
    {
#ifdef NGENIC_2LPT

      /* allocate temporary buffers for second derivatives */
      fft_real *d2phi1[3];
      for(int axes = 0; axes < 3; axes++)
        d2phi1[axes] = (fft_real *)Mem.mymalloc_clear("d2Phi1", maxfftsize * sizeof(fft_real));

      for(int axes = 0; axes < 3; axes++)
        {
          mpi_printf("NGENIC_2LPT: Computing secondary source term, derivatices %d %d\n", axes, axes);

          fft_real *disp = (fft_real *)Mem.mymalloc("disp", maxfftsize * sizeof(fft_real));

          ngenic_setup_modes_in_kspace((fft_complex *)disp, NULL, NULL);
          ngenic_get_derivate_from_fourier_field(axes, axes, (fft_complex *)disp);

          memcpy(d2phi1[axes], disp, maxfftsize * sizeof(fft_real));

          Mem.myfree(disp);
        }

      /* allocate second source potential */
      fft_real *Phi2 = (fft_real *)Mem.mymalloc_movable(&Phi2, "Phi2", maxfftsize * sizeof(fft_real));

      for(size_t n = 0; n < maxfftsize; n++)
        Phi2[n] = d2phi1[0][n] * d2phi1[1][n] + d2phi1[0][n] * d2phi1[2][n] + d2phi1[1][n] * d2phi1[2][n];

      for(int axes = 2; axes >= 0; axes--)
        Mem.myfree_movable(d2phi1[axes]);

      for(int i = 0; i < 3; i++)
        for(int j = i + 1; j < 3; j++)
          {
            mpi_printf("NGENIC_2LPT: Computing secondary source term, derivatices %d %d\n", i, j);

            fft_real *disp = (fft_real *)Mem.mymalloc("disp", maxfftsize * sizeof(fft_real));

            ngenic_setup_modes_in_kspace((fft_complex *)disp, NULL, NULL);
            ngenic_get_derivate_from_fourier_field(i, j, (fft_complex *)disp);

            for(size_t n = 0; n < maxfftsize; n++)
              Phi2[n] -= disp[n] * disp[n];

            Mem.myfree(disp);
          }

      mpi_printf("NGENIC_2LPT: Secondary source term computed in real space\n");

      /* Do a forward inplace-FFT to get Phi2 in Fourier space */
      ngenic_compute_transform_of_source_potential(Phi2);

      mpi_printf("NGENIC_2LPT: Done transforming it to k-space\n");

      for(int axes = 0; axes < 3; axes++)
        {
          mpi_printf("NGENIC_2LPT: Obtaining second order displacements for axes=%d\n", axes);

          fft_real *disp = (fft_real *)Mem.mymalloc("disp", maxfftsize * sizeof(fft_real));

          memcpy(disp, Phi2, maxfftsize * sizeof(fft_real));

          ngenic_get_derivate_from_fourier_field(axes, -1, (fft_complex *)disp);

          ngenic_readout_disp(disp, axes, 3.0 / 7, 3.0 / 7 * vel_prefac2, -1);

          Mem.myfree(disp);
        }

      Mem.myfree(Phi2);
#endif

      /* the mass variation IC does not require repeating through the vector directions indexed by axes */
#ifdef ADDITIONAL_GRID
      fft_real *deltaField = (fft_real *)Mem.mymalloc("deltaField", maxfftsize * sizeof(fft_real));

#ifndef CB_PHASE
      ngenic_setup_modes_in_kspace((fft_complex *)deltaField, NULL, NULL);  // set up the modes with the random seed 

      ngenic_get_delta_from_fourier_field((fft_complex *)deltaField);  // simply inverse transform the deltaField to get d(x)

      ngenic_readout_mass(deltaField);  // essentially the same function as readout_disp but only need to do it once since mass field is scalar
#else
      // if request CB_PHASE then read in the full 3D rhogrid from the restart snapshot and use the phase there instead
      // deltagrid holds all the info read in from the input file 
    
      ngenic_setup_modes_in_kspace((fft_complex *)deltaField, NULL, (fft_complex *)deltagrid_local);
 
      ngenic_get_delta_from_fourier_field((fft_complex *)deltaField);

      ngenic_readout_mass(deltaField); 
#endif

      MPI_Barrier(Communicator);
      
#endif

      /* now carry out Zeldovich approximation, yielding first order displacements */
      for(int axes = 0; axes < 3; axes++)
        {
          mpi_printf("NGENIC_2LPT: Obtaining Zeldovich displacements for axes=%d\n", axes);

          fft_real *disp = (fft_real *)Mem.mymalloc("disp", maxfftsize * sizeof(fft_real));

          fft_real *vel;
          vel = (fft_real *)Mem.mymalloc("vel", maxfftsize * sizeof(fft_real));

#ifdef CB_PHASE
          ngenic_setup_modes_in_kspace((fft_complex *)disp, (fft_complex *)vel, (fft_complex *)deltagrid_local);
#else
          ngenic_setup_modes_in_kspace((fft_complex *)disp, (fft_complex *)vel, NULL);
#endif 

          ngenic_get_derivate_from_fourier_field(axes, -1, (fft_complex *)disp);

          ngenic_get_derivate_from_fourier_field(axes, -1, (fft_complex *)vel);

          ngenic_readout_disp(disp, axes, 1.0, vel_prefac1, 0);

          ngenic_readout_disp(vel, axes, 1.0, vel_prefac1, 1);

          Mem.myfree(vel);

          Mem.myfree(disp);
        }

#ifdef ADDITIONAL_GRID
          Mem.myfree(deltaField);
#endif
    }

  /* now add displacement to Lagrangian coordinates  */
  double maxdisp = 0;
  double maxvel  = 0;

  int n_start = 0;

#ifdef ADDITIONAL_GRID 
  n_start = Sp->NumICPart;
  mpi_printf("skipped %d IC particles\n", n_start);
#endif

  int max_vel_part = 0;

  for(int n = n_start; n < Sp->NumPart; n++)
    {
      double posdiff[3] = {Pdisp[n-n_start].deltapos[0], Pdisp[n-n_start].deltapos[1], Pdisp[n-n_start].deltapos[2]};

      MyIntPosType delta[3];
      Sp->pos_to_signedintpos(posdiff, (MySignedIntPosType *)delta);

      for(int axes = 0; axes < 3; axes++)
        {
          Sp->P[n].IntPos[axes] += delta[axes];

          if(Pdisp[n-n_start].deltapos[axes] > maxdisp)
            maxdisp = Pdisp[n-n_start].deltapos[axes];

          //if(Sp->P[n].Vel[axes] > maxvel)
            //maxvel = Sp->P[n].Vel[axes];
          if(Sp->P[n].Vel[axes] > maxvel)
            {
              maxvel = Sp->P[n].Vel[axes];
              max_vel_part = n;
            }
        }

      //mpi_printf("mass[%d] = %f\n", n, Sp->P[n].getMass());
    }

  //printf("maxvel = %g on task %d\n", maxvel, ThisTask);

  double max_disp_global, maxvel_global;
  MPI_Reduce(&maxdisp, &max_disp_global, 1, MPI_DOUBLE, MPI_MAX, 0, Communicator);
  MPI_Reduce(&maxvel, &maxvel_global, 1, MPI_DOUBLE, MPI_MAX, 0, Communicator);

  double max_3D_vel = sqrt(Sp->P[max_vel_part].Vel[0] * Sp->P[max_vel_part].Vel[0] + Sp->P[max_vel_part].Vel[1] * Sp->P[max_vel_part].Vel[1] + Sp->P[max_vel_part].Vel[2] * Sp->P[max_vel_part].Vel[2]);
  printf("on Task = %d, max_3D_vel = %g, maxvel = %g\n", ThisTask, max_3D_vel, maxvel);

  mpi_printf("\nNGENIC: Maximum displacement: %g, in units of the part-spacing= %g\n\n", max_disp_global,
             max_disp_global / (All.BoxSize / NGENIC));
  mpi_printf("\nNGENIC: Maximum velocity component: %g\n\n", maxvel_global);

#ifdef THERMAL_VEL_IC
  Mem.myfree(thermal_vel_vec);
#endif
  Mem.myfree(Pvel);
  Mem.myfree(Pdisp);

#ifdef CB_PHASE
  Mem.myfree(deltagrid_local);
#endif

  Mem.myfree(partin);
  Mem.myfree(Rcvpm_offset);
  Mem.myfree(Rcvpm_count);
  Mem.myfree(Sndpm_offset);
  Mem.myfree(Sndpm_count);

  if(Shmem.Island_NTask != Shmem.World_NTask)
    {
      if(Shmem.Island_ThisTask == 0)
        {
          // need to send this flag to the correct processor rank (our shared memory handler) so that the table is freed there
          size_t tab_len = NGENIC * NGENIC * sizeof(unsigned int);
          MPI_Send(&tab_len, sizeof(tab_len), MPI_BYTE, Shmem.MyShmRankInGlobal, TAG_TABLE_FREE, MPI_COMM_WORLD);
        }
    }
  else
    {
      Mem.myfree(seedtable);
    }

#ifndef FFT_COLUMN_BASED
  my_slab_based_fft_free(&myplan);
#else
  my_column_based_fft_free(&myplan);
#endif

  FFTW(destroy_plan)(myplan.forward_plan_zdir);
  FFTW(destroy_plan)(myplan.forward_plan_ydir);
  FFTW(destroy_plan)(myplan.forward_plan_xdir);

  FFTW(destroy_plan)(myplan.backward_plan_zdir);
  FFTW(destroy_plan)(myplan.backward_plan_ydir);
  FFTW(destroy_plan)(myplan.backward_plan_xdir);

  if(All.PowerSpectrumType == 2)
    {
      free_f_growth_table();
#ifdef CB_PHASE
      free_CB_power_table();
#endif
      free_power_table();
    }

  gsl_rng_free(rnd_generator);
  gsl_rng_free(rnd_generator_conjugate);

  print_spec();

  TIMER_STOP(CPU_NGENIC);
}

void ngenic::ngenic_distribute_particles(void)
{
  Sndpm_count  = (size_t *)Mem.mymalloc("Sndpm_count", NTask * sizeof(size_t));
  Sndpm_offset = (size_t *)Mem.mymalloc("Sndpm_offset", NTask * sizeof(size_t));
  Rcvpm_count  = (size_t *)Mem.mymalloc("Rcvpm_count", NTask * sizeof(size_t));
  Rcvpm_offset = (size_t *)Mem.mymalloc("Rcvpm_offset", NTask * sizeof(size_t));

#ifdef FFT_COLUMN_BASED
  int columns         = GRIDX * GRIDY;
  int avg             = (columns - 1) / NTask + 1;
  int exc             = NTask * avg - columns;
  int tasklastsection = NTask - exc;
  int pivotcol        = tasklastsection * avg;
#endif

  /* determine the slabs/columns each particles accesses */
  {
    size_t *send_count = Sndpm_count;

    for(int j = 0; j < NTask; j++)
      send_count[j] = 0;

    int i_start = 0;

#ifdef ADDITIONAL_GRID
    i_start = Sp->NumICPart;
    mpi_printf("skipped %d IC particles\n", i_start);
#endif

    for(int i = i_start; i < Sp->NumPart; i++)
      {
        int slab_x  = Sp->P[i].IntPos[0] / INTCELL;
        int slab_xx = slab_x + 1;

        if(slab_xx >= GRIDX)
          slab_xx = 0;

#ifndef FFT_COLUMN_BASED
        int task0 = myplan.slab_to_task[slab_x];
        int task1 = myplan.slab_to_task[slab_xx];

        send_count[task0]++;
        if(task0 != task1)
          send_count[task1]++;
#else
        int slab_y  = Sp->P[i].IntPos[1] / INTCELL;
        int slab_yy = slab_y + 1;

        if(slab_yy >= GRIDY)
          slab_yy = 0;

        int column0 = slab_x * GRIDY + slab_y;
        int column1 = slab_x * GRIDY + slab_yy;
        int column2 = slab_xx * GRIDY + slab_y;
        int column3 = slab_xx * GRIDY + slab_yy;

        int task0, task1, task2, task3;

        if(column0 < pivotcol)
          task0 = column0 / avg;
        else
          task0 = (column0 - pivotcol) / (avg - 1) + tasklastsection;

        if(column1 < pivotcol)
          task1 = column1 / avg;
        else
          task1 = (column1 - pivotcol) / (avg - 1) + tasklastsection;

        if(column2 < pivotcol)
          task2 = column2 / avg;
        else
          task2 = (column2 - pivotcol) / (avg - 1) + tasklastsection;

        if(column3 < pivotcol)
          task3 = column3 / avg;
        else
          task3 = (column3 - pivotcol) / (avg - 1) + tasklastsection;

        send_count[task0]++;
        if(task1 != task0)
          send_count[task1]++;
        if(task2 != task1 && task2 != task0)
          send_count[task2]++;
        if(task3 != task0 && task3 != task1 && task3 != task2)
          send_count[task3]++;
#endif
      }
  }

  /* collect thread-specific offset table and collect the results from the other threads */
  Sndpm_offset[0] = 0;
  for(int i = 1; i < NTask; i++)
    {
      int ind      = i;
      int ind_prev = i - 1;

      Sndpm_offset[ind] = Sndpm_offset[ind_prev] + Sndpm_count[ind_prev];
    }

  MPI_Alltoall(Sndpm_count, sizeof(size_t), MPI_BYTE, Rcvpm_count, sizeof(size_t), MPI_BYTE, Communicator);

  nimport = 0, nexport = 0, Rcvpm_offset[0] = 0, Sndpm_offset[0] = 0;
  for(int j = 0; j < NTask; j++)
    {
      nexport += Sndpm_count[j];
      nimport += Rcvpm_count[j];

      if(j > 0)
        {
          Sndpm_offset[j] = Sndpm_offset[j - 1] + Sndpm_count[j - 1];
          Rcvpm_offset[j] = Rcvpm_offset[j - 1] + Rcvpm_count[j - 1];
        }
    }

  /* allocate import and export buffer */
  partin  = (partbuf *)Mem.mymalloc("partin", nimport * sizeof(partbuf));
  partout = (partbuf *)Mem.mymalloc("partout", nexport * sizeof(partbuf));

  {
    size_t *send_count  = Sndpm_count;
    size_t *send_offset = Sndpm_offset;

    for(int j = 0; j < NTask; j++)
      send_count[j] = 0;

    /* fill export buffer */

    int i_start = 0;

#ifdef ADDITIONAL_GRID 
    i_start = Sp->NumICPart;
    mpi_printf("skipped %d IC particles\n", i_start);
#endif

    for(int i = i_start; i < Sp->NumPart; i++)
      {
        int slab_x  = Sp->P[i].IntPos[0] / INTCELL;
        int slab_xx = slab_x + 1;

        if(slab_xx >= GRIDX)
          slab_xx = 0;

#ifndef FFT_COLUMN_BASED
        int task0 = myplan.slab_to_task[slab_x];
        int task1 = myplan.slab_to_task[slab_xx];

        size_t ind0 = send_offset[task0] + send_count[task0]++;
        for(int j = 0; j < 3; j++)
          partout[ind0].IntPos[j] = Sp->P[i].IntPos[j];

        if(task0 != task1)
          {
            size_t ind1 = send_offset[task1] + send_count[task1]++;
            for(int j = 0; j < 3; j++)
              partout[ind1].IntPos[j] = Sp->P[i].IntPos[j];
          }
#else
        int slab_y  = Sp->P[i].IntPos[1] / INTCELL;
        int slab_yy = slab_y + 1;

        if(slab_yy >= GRIDY)
          slab_yy = 0;

        int column0 = slab_x * GRIDY + slab_y;
        int column1 = slab_x * GRIDY + slab_yy;
        int column2 = slab_xx * GRIDY + slab_y;
        int column3 = slab_xx * GRIDY + slab_yy;

        int task0, task1, task2, task3;

        if(column0 < pivotcol)
          task0 = column0 / avg;
        else
          task0 = (column0 - pivotcol) / (avg - 1) + tasklastsection;

        if(column1 < pivotcol)
          task1 = column1 / avg;
        else
          task1 = (column1 - pivotcol) / (avg - 1) + tasklastsection;

        if(column2 < pivotcol)
          task2 = column2 / avg;
        else
          task2 = (column2 - pivotcol) / (avg - 1) + tasklastsection;

        if(column3 < pivotcol)
          task3 = column3 / avg;
        else
          task3 = (column3 - pivotcol) / (avg - 1) + tasklastsection;

        size_t ind0 = send_offset[task0] + send_count[task0]++;
        for(int j = 0; j < 3; j++)
          partout[ind0].IntPos[j] = Sp->P[i].IntPos[j];

        if(task1 != task0)
          {
            size_t ind1 = send_offset[task1] + send_count[task1]++;

            for(int j = 0; j < 3; j++)
              partout[ind1].IntPos[j] = Sp->P[i].IntPos[j];
          }
        if(task2 != task1 && task2 != task0)
          {
            size_t ind2 = send_offset[task2] + send_count[task2]++;

            for(int j = 0; j < 3; j++)
              partout[ind2].IntPos[j] = Sp->P[i].IntPos[j];
          }
        if(task3 != task0 && task3 != task1 && task3 != task2)
          {
            size_t ind3 = send_offset[task3] + send_count[task3]++;

            for(int j = 0; j < 3; j++)
              partout[ind3].IntPos[j] = Sp->P[i].IntPos[j];
          }
#endif
      }
  }

  int flag_big = 0, flag_big_all;
  for(int i = 0; i < NTask; i++)
    if(Sndpm_count[i] * sizeof(partbuf) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, Communicator);

  /* exchange particle data */
  myMPI_Alltoallv(partout, Sndpm_count, Sndpm_offset, partin, Rcvpm_count, Rcvpm_offset, sizeof(partbuf), flag_big_all, Communicator);

  Mem.myfree(partout);
}

void ngenic::ngenic_compute_transform_of_source_potential(fft_real *pot)
{
  fft_real *workspace = (fft_real *)Mem.mymalloc("workspace", maxfftsize * sizeof(fft_real));

#ifndef FFT_COLUMN_BASED
  my_slab_based_fft(&myplan, &pot[0], &workspace[0], +1);
#else
  my_column_based_fft(&myplan, pot, workspace, +1);  // result is in workspace, not in Phi2
  memcpy(pot, workspace, maxfftsize * sizeof(fft_real));
#endif

  Mem.myfree(workspace);

  double normfac = 1 / (((double)GRIDX) * GRIDY * GRIDZ);

  for(size_t n = 0; n < maxfftsize; n++)
    pot[n] *= normfac;
}

void ngenic::ngenic_get_delta_from_fourier_field(fft_complex *fft_of_grid)
{
#ifdef FFT_COLUMN_BASED
  if(myplan.second_transposed_firstcol == 0)
    fft_of_grid[0][0] = fft_of_grid[0][1] = 0.0;
#else
  if(myplan.slabstart_y == 0)
    fft_of_grid[0][0] = fft_of_grid[0][1] = 0.0;
#endif

  /* Do the inverse FFT to get the displacement field */
  fft_real *workspace = (fft_real *)Mem.mymalloc("workspace", maxfftsize * sizeof(fft_real));

#ifndef FFT_COLUMN_BASED
  my_slab_based_fft(&myplan, &fft_of_grid[0], &workspace[0], -1);
#else
  my_column_based_fft(&myplan, fft_of_grid, workspace, -1);  // result is in workspace
  memcpy(fft_of_grid, workspace, maxfftsize * sizeof(fft_real));
#endif

  Mem.myfree(workspace);
}

/* this function returns the component 'axes' (0, 1, or 2) of the gradient of a field phi,
 * which is the solution of nabla^2 phi = grid.
 * We input the Fourier transform of grid to the function, and this field is overwritten with
 * the gradient.
 */
void ngenic::ngenic_get_derivate_from_fourier_field(int axes1, int axes2, fft_complex *fft_of_grid)
{
  double kfacx = 2.0 * M_PI / All.BoxSize;
  double kfacy = 2.0 * M_PI / All.BoxSize;
  double kfacz = 2.0 * M_PI / All.BoxSize;

#ifdef FFT_COLUMN_BASED
  for(large_array_offset ip = 0; ip < myplan.second_transposed_ncells; ip += GRIDX)
    {
      large_array_offset ipcell = ip + ((large_array_offset)myplan.second_transposed_firstcol) * GRIDX;
      int y                     = ipcell / (GRIDX * GRIDz);
      int yr                    = ipcell % (GRIDX * GRIDz);
      int z                     = yr / GRIDX;
      if(yr % GRIDX != 0)  // Note: check that x-columns are really complete
        Terminate("x-column seems incomplete. This is not expected");
#else
  for(int y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
    for(int z = 0; z < GRIDz; z++)
      {
#endif

      for(int x = 0; x < GRIDX; x++)
        {
          int xx, yy, zz;

          if(x >= (GRIDX / 2))
            xx = x - GRIDX;
          else
            xx = x;
          if(y >= (GRIDY / 2))
            yy = y - GRIDY;
          else
            yy = y;
          if(z >= (GRIDZ / 2))
            zz = z - GRIDZ;
          else
            zz = z;

          double kvec[3];
          kvec[0] = kfacx * xx;
          kvec[1] = kfacy * yy;
          kvec[2] = kfacz * zz;

          double kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];

          double smth = 1;

#ifdef CORRECT_CIC
          if(axes2 >= 0)
            {
              /* do deconvolution of CIC interpolation */
              double fx = 1, fy = 1, fz = 1;
              if(kvec[0] != 0)
                {
                  fx = (kvec[0] * All.BoxSize / 2) / NGENIC;
                  fx = sin(fx) / fx;
                }
              if(kvec[1] != 0)
                {
                  fy = (kvec[1] * All.BoxSize / 2) / NGENIC;
                  fy = sin(fy) / fy;
                }
              if(kvec[2] != 0)
                {
                  fz = (kvec[2] * All.BoxSize / 2) / NGENIC;
                  fz = sin(fz) / fz;
                }
              double ff = 1 / (fx * fy * fz);
              smth      = ff * ff;
              /* end deconvolution */
            }
#endif

#ifndef FFT_COLUMN_BASED
          large_array_offset elem = ((large_array_offset)GRIDz) * (GRIDX * (y - myplan.slabstart_y) + x) + z;
#else
            large_array_offset elem = ip + x;
#endif

          fft_real re = smth * fft_of_grid[elem][0];
          fft_real im = smth * fft_of_grid[elem][1];

          if(axes2 < 0)
            {
              /* first derivative */
              fft_of_grid[elem][0] = (kmag2 > 0.0 ? -kvec[axes1] / kmag2 * im : 0.0);
              fft_of_grid[elem][1] = (kmag2 > 0.0 ? kvec[axes1] / kmag2 * re : 0.0);
            }
          else
            {
              /* second derivative */
              fft_of_grid[elem][0] = (kmag2 > 0.0 ? kvec[axes1] * kvec[axes2] / kmag2 * re : 0.0);
              fft_of_grid[elem][1] = (kmag2 > 0.0 ? kvec[axes1] * kvec[axes2] / kmag2 * im : 0.0);
            }
        }
    }

#ifdef FFT_COLUMN_BASED
  if(myplan.second_transposed_firstcol == 0)
    fft_of_grid[0][0] = fft_of_grid[0][1] = 0.0;
#else
  if(myplan.slabstart_y == 0)
    fft_of_grid[0][0] = fft_of_grid[0][1] = 0.0;
#endif

  /* Do the inverse FFT to get the displacement field */
  fft_real *workspace = (fft_real *)Mem.mymalloc("workspace", maxfftsize * sizeof(fft_real));

#ifndef FFT_COLUMN_BASED
  my_slab_based_fft(&myplan, &fft_of_grid[0], &workspace[0], -1);
#else
  my_column_based_fft(&myplan, fft_of_grid, workspace, -1);  // result is in workspace
  memcpy(fft_of_grid, workspace, maxfftsize * sizeof(fft_real));
#endif

  Mem.myfree(workspace);
}

void ngenic::ngenic_setup_modes_in_kspace(fft_complex *fft_of_grid, fft_complex *fft_of_vel, fft_complex *fft_of_delta)
{
  double fac = pow(1. / All.BoxSize, 1.5);  /* removed the factor of (2*Pi)^1.5 from the stock version, as the power spectrum input is now in units of (Mpc/h)^3. */

  /* clear local FFT-mesh */
  memset(fft_of_grid, 0, maxfftsize * sizeof(fft_real));

  if(fft_of_vel != NULL) 
    memset(fft_of_vel, 0, maxfftsize * sizeof(fft_real)); 

  mpi_printf("NGENIC: setting up modes in kspace...\n");

  double kfacx = 2.0 * M_PI / All.BoxSize;
  double kfacy = 2.0 * M_PI / All.BoxSize;
  double kfacz = 2.0 * M_PI / All.BoxSize;

//  printf("Task = %d, slabstart_y = %d, nslab_y = %d\n", ThisTask, myplan.slabstart_y, myplan.nslab_y);

#ifdef FFT_COLUMN_BASED
  for(large_array_offset ip = 0; ip < myplan.second_transposed_ncells; ip += GRIDX)
    {
      large_array_offset ipcell = ip + ((large_array_offset)myplan.second_transposed_firstcol) * GRIDX;
      int y                     = ipcell / (GRIDX * GRIDz);
      int yr                    = ipcell % (GRIDX * GRIDz);
      int z                     = yr / GRIDX;
      if(yr % GRIDX != 0)  // Note: check that x-columns are really complete
        Terminate("x-column seems incomplete. This is not expected");
#else
  for(int y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
    for(int z = 0; z < GRIDz; z++)
      {
#endif

      // let's use the y and z plane here, because the x-column is available in full for both FFT schemes
      gsl_rng_set(rnd_generator, seedtable[y * NGENIC + z]);

      // we also create the modes for the conjugate column so that we can fulfill the reality constraint
      // by using the conjugate of the corresponding mode if needed
      int y_conj, z_conj;
      if(y > 0)
        y_conj = GRIDY - y;
      else
        y_conj = 0;

      if(z > 0)
        z_conj = GRIDZ - z;
      else
        z_conj = 0;

      gsl_rng_set(rnd_generator_conjugate, seedtable[y_conj * NGENIC + z_conj]);

#ifndef NGENIC_FIX_MODE_AMPLITUDES
      double mode_ampl[GRIDX], mode_ampl_conj[GRIDX];
#endif
      double mode_phase[GRIDX], mode_phase_conj[GRIDX];

      // in this loop we precompute the modes for both columns, from low-k to high-k,
      // so that after an increase of resolution, one gets the same modes plus new ones
      for(int xoff = 0; xoff < GRIDX / 2; xoff++)
        for(int side = 0; side < 2; side++)
          {
            int x;
            if(side == 0)
              x = xoff;
            else
              x = GRIDX - 1 - xoff;

            double phase      = gsl_rng_uniform(rnd_generator) * 2 * M_PI;
            double phase_conj = gsl_rng_uniform(rnd_generator_conjugate) * 2 * M_PI;

#ifdef NGENIC_MIRROR_PHASES
            phase += M_PI;
            if(phase >= 2 * M_PI)
              phase -= 2 * M_PI;

            phase_conj += M_PI;
            if(phase_conj >= 2 * M_PI)
              phase_conj -= 2 * M_PI;
#endif
            mode_phase[x]      = phase;
            mode_phase_conj[x] = phase_conj;

#ifndef NGENIC_FIX_MODE_AMPLITUDES
            double ampl;
            do
              {
                ampl = gsl_rng_uniform(rnd_generator);
              }
            while(ampl == 0);

            double ampl_conj;
            do
              {
                ampl_conj = gsl_rng_uniform(rnd_generator_conjugate);
              }
            while(ampl_conj == 0);

            mode_ampl[x] = ampl;

            mode_ampl_conj[x] = ampl_conj;
#endif
          }

      // now let's populate the full x-column of modes
      for(int x = 0; x < GRIDX; x++)
        {
          int xx, yy, zz;

          if(x >= (GRIDX / 2))
            xx = x - GRIDX;
          else
            xx = x;
          if(y >= (GRIDY / 2))
            yy = y - GRIDY;
          else
            yy = y;
          if(z >= (GRIDZ / 2)) 
            zz = z - GRIDZ;
          else
            zz = z;

          double kvec[3];
          kvec[0] = kfacx * xx;
          kvec[1] = kfacy * yy;
          kvec[2] = kfacz * zz;

          double kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
          double kmag  = sqrt(kmag2);

          if(All.SphereMode == 1)
            {
              if(kmag * All.BoxSize / (2 * M_PI) > All.NSample / 2) /* select a sphere in k-space */
                continue;
            }
          else
            {
              if(fabs(kvec[0]) * All.BoxSize / (2 * M_PI) > All.NSample / 2)
                continue;
              if(fabs(kvec[1]) * All.BoxSize / (2 * M_PI) > All.NSample / 2)
                continue;
              if(fabs(kvec[2]) * All.BoxSize / (2 * M_PI) > All.NSample / 2)
                continue;
            }

          double p_of_k = ngenic_power_spec(kmag);

#ifdef CB_PHASE
          double pcb_of_k = ngenic_CB_power_spec(kmag);
#endif

          double f_of_k = 0.;
          if(fft_of_vel != NULL)
            {
              f_of_k = ngenic_growth_rate(kmag);
            }

          /* Note: kmag and p_of_k are unaffected by whether or not we use the conjugate mode */

          int conjugate_flag = 0;

          if(z == 0 || z == GRIDZ / 2)
            {
              if(x > GRIDX / 2 && x < GRIDX)
                conjugate_flag = 1;
              else if(x == 0 || x == GRIDX / 2)
                {
                  if(y > GRIDY / 2 && y < GRIDY)
                    conjugate_flag = 1;
                  else if(y == 0 || y == GRIDX / 2)
                    {
                      continue;
                    }
                }
            }

          // determine location of conjugate mode in x column
          int x_conj;

          if(x > 0)
            x_conj = GRIDX - x;
          else
            x_conj = 0;

#ifndef NGENIC_FIX_MODE_AMPLITUDES
          if(conjugate_flag)
            p_of_k *= -log(mode_ampl_conj[x_conj]);
          else
            p_of_k *= -log(mode_ampl[x]);

/*
#ifdef CB_PHASE
	  if(conjugate_flag)
	    pcb_of_k *= -log(mode_ampl_conj[x_conj]);
          else
            pcb_of_k *= -log(mode_ampl[x]);
#endif
*/
#endif

#ifndef CB_PHASE
          double delta = fac * sqrt(p_of_k) / Dplus; /* scale back to starting redshift */
#else
          /* if CB_PHASE is requested, ignore the gaussian random field, the phases are imported from the restart package */
          /* the -log(mode_ampl) factor is dropped here since that is part of the Box-Muller transform */
          double delta = fac * sqrt(ngenic_power_spec(kmag)); //fac * sqrt(ngenic_power_spec(kmag)); // delta here is the d_nu(k) 
          double delta_cb = fac * sqrt(ngenic_CB_power_spec(kmag));

          //mpi_printf("delta_cb = %g, pcb_of_k = %g\n", delta_cb, pcb_of_k);

          /* fft_of_delta[elem] / (fac * sqrt(ngenic_CB_power_spec(kmag))) is the gaussian random field with mean 0 and std dev 1 */
#endif

#ifndef FFT_COLUMN_BASED
          large_array_offset elem = ((large_array_offset)GRIDz) * (GRIDX * (y - myplan.slabstart_y) + x) + z;
#else
            large_array_offset elem = ip + x;
#endif

          if(conjugate_flag)
            {
#ifndef CB_PHASE
              fft_of_grid[elem][0] = delta * cos(mode_phase_conj[x_conj]);
              fft_of_grid[elem][1] = -delta * sin(mode_phase_conj[x_conj]);
#else
	      if(fft_of_delta != NULL)
                {
                  if(delta_cb == 0)
                    {
                      fft_of_grid[elem][0] = 0.0;
                      fft_of_grid[elem][1] = 0.0;
                    } else {
                      fft_of_grid[elem][0] = delta * fft_of_delta[elem][0] / delta_cb;
                      fft_of_grid[elem][1] = delta * fft_of_delta[elem][1] / delta_cb;
                    }
                } else
                {
                  // this part is really only here to make sure the disp field doesn't call a random patch of memory
                  // it's not doing anything, I'm just too lazy to restructure the whole thing 
                  // if all goees right, this part of code should never be called
                  fft_of_grid[elem][0] = delta * cos(mode_phase_conj[x_conj]);
                  fft_of_grid[elem][1] = -delta * sin(mode_phase_conj[x_conj]);
                }
#endif

              if(fft_of_vel != NULL)
                {
#ifndef CB_PHASE
                  fft_of_vel[elem][0] = delta * f_of_k * cos(mode_phase_conj[x_conj]);
                  fft_of_vel[elem][1] = -delta * f_of_k * sin(mode_phase_conj[x_conj]);
#else
                  if(fft_of_delta != NULL)
                    {
		      if(delta_cb == 0) 
                        {
                          fft_of_vel[elem][0] = 0.0;
                          fft_of_vel[elem][1] = 0.0;
                        } else {
                          fft_of_vel[elem][0] = delta * f_of_k * fft_of_delta[elem][0] / delta_cb;
                          fft_of_vel[elem][1] = delta * f_of_k * fft_of_delta[elem][1] / delta_cb;
		        }
                    } else
                    {
		      printf("wrong place\n");
                    }
#endif
                }
            }
          else
            {
#ifndef CB_PHASE
              fft_of_grid[elem][0] = delta * cos(mode_phase[x]);
              fft_of_grid[elem][1] = delta * sin(mode_phase[x]);
#else
              if(fft_of_delta != NULL)
                {
                  if(delta_cb == 0)
                    {
                      fft_of_grid[elem][0] = 0.0;
                      fft_of_grid[elem][1] = 0.0;
                    } else {
                      fft_of_grid[elem][0] = delta * fft_of_delta[elem][0] / delta_cb;
                      fft_of_grid[elem][1] = delta * fft_of_delta[elem][1] / delta_cb;
		    }
                } else 
                {
                  fft_of_grid[elem][0] = delta * cos(mode_phase[x]);
                  fft_of_grid[elem][1] = delta * sin(mode_phase[x]);
                }
#endif

              if(fft_of_vel != NULL)
                {
#ifndef CB_PHASE
                  fft_of_vel[elem][0] = delta * f_of_k * cos(mode_phase[x]);
                  fft_of_vel[elem][1] = delta * f_of_k * sin(mode_phase[x]);
#else
                  if(fft_of_delta != NULL)
                    {
		      if(delta_cb == 0) 
 		        {
			  fft_of_vel[elem][0] = 0.0;
                          fft_of_vel[elem][1] = 0.0;
    			} else {
                          fft_of_vel[elem][0] = delta * f_of_k * fft_of_delta[elem][0] / delta_cb;
                          fft_of_vel[elem][1] = delta * f_of_k * fft_of_delta[elem][1] / delta_cb;
                        }
                    } else 
                    {
		      printf("wrong place\n");
                    }
#endif
                }
            }
        }
    }
}

#ifdef ADDITIONAL_GRID
void ngenic::ngenic_readout_mass(fft_real *grid)
{
#ifdef FFT_COLUMN_BASED
  int columns         = GRIDX * GRIDY;
  int avg             = (columns - 1) / NTask + 1;
  int exc             = NTask * avg - columns;
  int tasklastsection = NTask - exc;
  int pivotcol        = tasklastsection * avg;
#endif

  long long gridSize    = All.GridSize;
  long long partTotal   = gridSize * gridSize * gridSize;
  long long partPerTask = partTotal / NTask;

  double masstot = 0;

  if(All.N_tau_part > 0)
    {
      masstot = (All.OmegaNuPart / All.N_tau_part) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G) * All.BoxSize * All.BoxSize * All.BoxSize;
    } else {
      mpi_printf("no particle neutrinos are present\n");
      masstot = 0;
    }

  double m = masstot / (partTotal);

  mpi_printf("mass = %g\n", m);

  double *flistin  = (double *)Mem.mymalloc("flistin", nimport * sizeof(double));
  double *flistout = (double *)Mem.mymalloc("flistout", nexport * sizeof(double));

  for(size_t i = 0; i < nimport; i++)
    {
      flistin[i] = 0;

      int slab_x         = partin[i].IntPos[0] / INTCELL;
      MyIntPosType rmd_x = partin[i].IntPos[0] % INTCELL;

      int slab_y         = partin[i].IntPos[1] / INTCELL;
      MyIntPosType rmd_y = partin[i].IntPos[1] % INTCELL;

      int slab_z         = partin[i].IntPos[2] / INTCELL;
      MyIntPosType rmd_z = partin[i].IntPos[2] % INTCELL;

      double dx = rmd_x * (1.0 / INTCELL);
      double dy = rmd_y * (1.0 / INTCELL);
      double dz = rmd_z * (1.0 / INTCELL);

      int slab_xx = slab_x + 1;
      int slab_yy = slab_y + 1;
      int slab_zz = slab_z + 1;

      if(slab_xx >= GRIDX)
        slab_xx = 0;
      if(slab_yy >= GRIDY)
        slab_yy = 0;
      if(slab_zz >= GRIDZ)
        slab_zz = 0;

#ifndef FFT_COLUMN_BASED
      if(myplan.slab_to_task[slab_x] == ThisTask)
        {
          slab_x -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += grid[FI(slab_x, slab_y, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FI(slab_x, slab_y, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz) +
                        grid[FI(slab_x, slab_yy, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) +
                        grid[FI(slab_x, slab_yy, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(myplan.slab_to_task[slab_xx] == ThisTask)
        {
          slab_xx -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += grid[FI(slab_xx, slab_y, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FI(slab_xx, slab_y, slab_zz)] * (dx) * (1.0 - dy) * (dz) +
                        grid[FI(slab_xx, slab_yy, slab_z)] * (dx) * (dy) * (1.0 - dz) +
                        grid[FI(slab_xx, slab_yy, slab_zz)] * (dx) * (dy) * (dz);
        }
#else
      int column0 = slab_x * GRIDY + slab_y;
      int column1 = slab_x * GRIDY + slab_yy;
      int column2 = slab_xx * GRIDY + slab_y;
      int column3 = slab_xx * GRIDY + slab_yy;

      if(column0 >= myplan.firstcol_XY && column0 <= myplan.lastcol_XY)
        {
          flistin[i] += grid[FC(column0, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FC(column0, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz);
        }
      if(column1 >= myplan.firstcol_XY && column1 <= myplan.lastcol_XY)
        {
          flistin[i] +=
              grid[FC(column1, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) + grid[FC(column1, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(column2 >= myplan.firstcol_XY && column2 <= myplan.lastcol_XY)
        {
          flistin[i] +=
              grid[FC(column2, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) + grid[FC(column2, slab_zz)] * (dx) * (1.0 - dy) * (dz);
        }

      if(column3 >= myplan.firstcol_XY && column3 <= myplan.lastcol_XY)
        {
          flistin[i] += grid[FC(column3, slab_z)] * (dx) * (dy) * (1.0 - dz) + grid[FC(column3, slab_zz)] * (dx) * (dy) * (dz);
        }
#endif
    }

  /* exchange the potential component data */
  int flag_big = 0, flag_big_all;
  for(int i = 0; i < NTask; i++)
    if(Sndpm_count[i] * sizeof(double) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
 *    * transfer the data in chunks.
 *       */
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, Communicator);

  /* exchange  data */
  myMPI_Alltoallv(flistin, Rcvpm_count, Rcvpm_offset, flistout, Sndpm_count, Sndpm_offset, sizeof(double), flag_big_all, Communicator);

  {
    size_t *send_count  = Sndpm_count;
    size_t *send_offset = Sndpm_offset;

    for(int j = 0; j < NTask; j++)
      send_count[j] = 0;

    int i_start = 0;

    i_start = Sp->NumICPart;
    mpi_printf("skipped %d IC particles\n", i_start);

    for(int i = i_start; i < Sp->NumPart; i++)
      {
        int slab_x  = Sp->P[i].IntPos[0] / INTCELL;
        int slab_xx = slab_x + 1;

        if(slab_xx >= GRIDX)
          slab_xx = 0;

#ifndef FFT_COLUMN_BASED
        int task0 = myplan.slab_to_task[slab_x];
        int task1 = myplan.slab_to_task[slab_xx];

        double value = flistout[send_offset[task0] + send_count[task0]++];

        if(task0 != task1)
          value += flistout[send_offset[task1] + send_count[task1]++];

#else
        int slab_y  = Sp->P[i].IntPos[1] / INTCELL;
        int slab_yy = slab_y + 1;

        if(slab_yy >= GRIDY)
          slab_yy = 0;

        int column0 = slab_x * GRIDY + slab_y;
        int column1 = slab_x * GRIDY + slab_yy;
        int column2 = slab_xx * GRIDY + slab_y;
        int column3 = slab_xx * GRIDY + slab_yy;

        int task0, task1, task2, task3;

        if(column0 < pivotcol)
          task0 = column0 / avg;
        else
          task0 = (column0 - pivotcol) / (avg - 1) + tasklastsection;

        if(column1 < pivotcol)
          task1 = column1 / avg;
        else
          task1 = (column1 - pivotcol) / (avg - 1) + tasklastsection;

        if(column2 < pivotcol)
          task2 = column2 / avg;
        else
          task2 = (column2 - pivotcol) / (avg - 1) + tasklastsection;

        if(column3 < pivotcol)
          task3 = column3 / avg;
        else
          task3 = (column3 - pivotcol) / (avg - 1) + tasklastsection;

        double value = flistout[send_offset[task0] + send_count[task0]++];

        if(task1 != task0)
          value += flistout[send_offset[task1] + send_count[task1]++];

        if(task2 != task1 && task2 != task0)
          value += flistout[send_offset[task2] + send_count[task2]++];

        if(task3 != task0 && task3 != task1 && task3 != task2)
          value += flistout[send_offset[task3] + send_count[task3]++];
#endif

	Sp->P[i].setMass((1.+value)*m);
      }
  }

  Mem.myfree(flistout);
  Mem.myfree(flistin);

}
#endif

void ngenic::ngenic_readout_disp(fft_real *grid, int axis, double pfac, double vfac, int vel_option)
{
#ifdef FFT_COLUMN_BASED
  int columns         = GRIDX * GRIDY;
  int avg             = (columns - 1) / NTask + 1;
  int exc             = NTask * avg - columns;
  int tasklastsection = NTask - exc;
  int pivotcol        = tasklastsection * avg;
#endif

  double *flistin  = (double *)Mem.mymalloc("flistin", nimport * sizeof(double));
  double *flistout = (double *)Mem.mymalloc("flistout", nexport * sizeof(double));

  for(size_t i = 0; i < nimport; i++)
    {
      flistin[i] = 0;

      int slab_x         = partin[i].IntPos[0] / INTCELL;
      MyIntPosType rmd_x = partin[i].IntPos[0] % INTCELL;

      int slab_y         = partin[i].IntPos[1] / INTCELL;
      MyIntPosType rmd_y = partin[i].IntPos[1] % INTCELL;

      int slab_z         = partin[i].IntPos[2] / INTCELL;
      MyIntPosType rmd_z = partin[i].IntPos[2] % INTCELL;

      double dx = rmd_x * (1.0 / INTCELL);
      double dy = rmd_y * (1.0 / INTCELL);
      double dz = rmd_z * (1.0 / INTCELL);

      int slab_xx = slab_x + 1;
      int slab_yy = slab_y + 1;
      int slab_zz = slab_z + 1;

      if(slab_xx >= GRIDX)
        slab_xx = 0;
      if(slab_yy >= GRIDY)
        slab_yy = 0;
      if(slab_zz >= GRIDZ)
        slab_zz = 0;

#ifndef FFT_COLUMN_BASED
      if(myplan.slab_to_task[slab_x] == ThisTask)
        {
          slab_x -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += grid[FI(slab_x, slab_y, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FI(slab_x, slab_y, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz) +
                        grid[FI(slab_x, slab_yy, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) +
                        grid[FI(slab_x, slab_yy, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(myplan.slab_to_task[slab_xx] == ThisTask)
        {
          slab_xx -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += grid[FI(slab_xx, slab_y, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FI(slab_xx, slab_y, slab_zz)] * (dx) * (1.0 - dy) * (dz) +
                        grid[FI(slab_xx, slab_yy, slab_z)] * (dx) * (dy) * (1.0 - dz) +
                        grid[FI(slab_xx, slab_yy, slab_zz)] * (dx) * (dy) * (dz);
        }
#else
      int column0 = slab_x * GRIDY + slab_y;
      int column1 = slab_x * GRIDY + slab_yy;
      int column2 = slab_xx * GRIDY + slab_y;
      int column3 = slab_xx * GRIDY + slab_yy;

      if(column0 >= myplan.firstcol_XY && column0 <= myplan.lastcol_XY)
        {
          flistin[i] += grid[FC(column0, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FC(column0, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz);
        }
      if(column1 >= myplan.firstcol_XY && column1 <= myplan.lastcol_XY)
        {
          flistin[i] +=
              grid[FC(column1, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) + grid[FC(column1, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(column2 >= myplan.firstcol_XY && column2 <= myplan.lastcol_XY)
        {
          flistin[i] +=
              grid[FC(column2, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) + grid[FC(column2, slab_zz)] * (dx) * (1.0 - dy) * (dz);
        }

      if(column3 >= myplan.firstcol_XY && column3 <= myplan.lastcol_XY)
        {
          flistin[i] += grid[FC(column3, slab_z)] * (dx) * (dy) * (1.0 - dz) + grid[FC(column3, slab_zz)] * (dx) * (dy) * (dz);
        }
#endif
    }

  /* exchange the potential component data */
  int flag_big = 0, flag_big_all;
  for(int i = 0; i < NTask; i++)
    if(Sndpm_count[i] * sizeof(double) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, Communicator);

  /* exchange  data */
  myMPI_Alltoallv(flistin, Rcvpm_count, Rcvpm_offset, flistout, Sndpm_count, Sndpm_offset, sizeof(double), flag_big_all, Communicator);

  {
    size_t *send_count  = Sndpm_count;
    size_t *send_offset = Sndpm_offset;

    for(int j = 0; j < NTask; j++)
      send_count[j] = 0;

    int i_start = 0;

#ifdef ADDITIONAL_GRID
    i_start = Sp->NumICPart;
    mpi_printf("skipped %d IC particles\n", i_start);
#endif

    for(int i = i_start; i < Sp->NumPart; i++)
      {
        int slab_x  = Sp->P[i].IntPos[0] / INTCELL;
        int slab_xx = slab_x + 1;

        if(slab_xx >= GRIDX)
          slab_xx = 0;

#ifndef FFT_COLUMN_BASED
        int task0 = myplan.slab_to_task[slab_x];
        int task1 = myplan.slab_to_task[slab_xx];

        double value = flistout[send_offset[task0] + send_count[task0]++];

        if(task0 != task1)
          value += flistout[send_offset[task1] + send_count[task1]++];

#else
        int slab_y  = Sp->P[i].IntPos[1] / INTCELL;
        int slab_yy = slab_y + 1;

        if(slab_yy >= GRIDY)
          slab_yy = 0;

        int column0 = slab_x * GRIDY + slab_y;
        int column1 = slab_x * GRIDY + slab_yy;
        int column2 = slab_xx * GRIDY + slab_y;
        int column3 = slab_xx * GRIDY + slab_yy;

        int task0, task1, task2, task3;

        if(column0 < pivotcol)
          task0 = column0 / avg;
        else
          task0 = (column0 - pivotcol) / (avg - 1) + tasklastsection;

        if(column1 < pivotcol)
          task1 = column1 / avg;
        else
          task1 = (column1 - pivotcol) / (avg - 1) + tasklastsection;

        if(column2 < pivotcol)
          task2 = column2 / avg;
        else
          task2 = (column2 - pivotcol) / (avg - 1) + tasklastsection;

        if(column3 < pivotcol)
          task3 = column3 / avg;
        else
          task3 = (column3 - pivotcol) / (avg - 1) + tasklastsection;

        double value = flistout[send_offset[task0] + send_count[task0]++];

        if(task1 != task0)
          value += flistout[send_offset[task1] + send_count[task1]++];

        if(task2 != task1 && task2 != task0)
          value += flistout[send_offset[task2] + send_count[task2]++];

        if(task3 != task0 && task3 != task1 && task3 != task2)
          value += flistout[send_offset[task3] + send_count[task3]++];
#endif

        if(vel_option < 0)
          {
            Pdisp[i].deltapos[axis] += pfac * value;
            Sp->P[i].Vel[axis] += vfac * value;
          }

        if(vel_option == 0)
          {
#ifndef ADDITIONAL_GRID
            Pdisp[i].deltapos[axis] += pfac * value;
#else
            Pdisp[i].deltapos[axis] += 0.0;
#endif
          }

        if(vel_option == 1)
          {
#if defined (ADDITIONAL_GRID) && defined (THERMAL_VEL_IC)
            Sp->P[i].Vel[axis] += thermal_vel_vec[3*(i-Sp->NumICPart)+axis];
#endif
            Sp->P[i].Vel[axis] += vfac * value;
          }
      }
  }

  Mem.myfree(flistout);
  Mem.myfree(flistin);

}

void ngenic::ngenic_initialize_ffts(void)
{
#ifdef LONG_X
  if(LONG_X != (int)(LONG_X))
    Terminate("LONG_X must be an integer if used with PMGRID");
#endif

#ifdef LONG_Y
  if(LONG_Y != (int)(LONG_Y))
    Terminate("LONG_Y must be an integer if used with PMGRID");
#endif

#ifdef LONG_Z
  if(LONG_Z != (int)(LONG_Z))
    Terminate("LONG_Z must be an integer if used with PMGRID");
#endif

  /* Set up the FFTW-3 plan files. */
  int ndimx[1] = {GRIDX}; /* dimension of the 1D transforms */
  int ndimy[1] = {GRIDY}; /* dimension of the 1D transforms */
  int ndimz[1] = {GRIDZ}; /* dimension of the 1D transforms */

  int max_GRID2 = 2 * (std::max<int>(std::max<int>(GRIDX, GRIDY), GRIDZ) / 2 + 1);

  /* temporarily allocate some arrays to make sure that out-of-place plans are created */
  fft_real *DispGrid           = (fft_real *)Mem.mymalloc("DispGrid", max_GRID2 * sizeof(fft_real));
  fft_complex *fft_of_DispGrid = (fft_complex *)Mem.mymalloc("DispGrid", max_GRID2 * sizeof(fft_real));

#ifdef DOUBLEPRECISION_FFTW
  int alignflag = 0;
#else
  /* for single precision, the start of our FFT columns is presently only guaranteed to be 8-byte aligned */
  int alignflag = FFTW_UNALIGNED;
#endif

#ifndef FFT_COLUMN_BASED
  int stride = GRIDz;
#else
  int stride    = 1;
#endif

  myplan.backward_plan_xdir =
      FFTW(plan_many_dft)(1, ndimx, 1, (fft_complex *)DispGrid, 0, stride, GRIDz * GRIDX, fft_of_DispGrid, 0, stride, GRIDz * GRIDX,
                          FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_ydir =
      FFTW(plan_many_dft)(1, ndimy, 1, (fft_complex *)DispGrid, 0, stride, GRIDz * GRIDY, fft_of_DispGrid, 0, stride, GRIDz * GRIDY,
                          FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_zdir = FFTW(plan_many_dft_c2r)(1, ndimz, 1, (fft_complex *)DispGrid, 0, 1, GRIDz, (fft_real *)fft_of_DispGrid,
                                                      0, 1, GRID2, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.forward_plan_xdir = FFTW(plan_many_dft)(1, ndimx, 1, (fft_complex *)DispGrid, 0, stride, GRIDz * GRIDX, fft_of_DispGrid, 0,
                                                 stride, GRIDz * GRIDX, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.forward_plan_ydir = FFTW(plan_many_dft)(1, ndimy, 1, (fft_complex *)DispGrid, 0, stride, GRIDz * GRIDY, fft_of_DispGrid, 0,
                                                 stride, GRIDz * GRIDY, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.forward_plan_zdir = FFTW(plan_many_dft_r2c)(1, ndimz, 1, DispGrid, 0, 1, GRID2, (fft_complex *)fft_of_DispGrid, 0, 1, GRIDz,
                                                     FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  Mem.myfree(fft_of_DispGrid);
  Mem.myfree(DispGrid);

#ifndef FFT_COLUMN_BASED

  my_slab_based_fft_init(&myplan, GRIDX, GRIDY, GRIDZ);

  maxfftsize = std::max<int>(myplan.largest_x_slab * GRIDY, myplan.largest_y_slab * GRIDX) * ((size_t)GRID2);

#else

  my_column_based_fft_init(&myplan, GRIDX, GRIDY, GRIDZ);

  maxfftsize = myplan.max_datasize;

#endif
}

void ngenic::print_spec(void)
{
  if(ThisTask == 0)
    {
      char buf[3 * MAXLEN_PATH];
      sprintf(buf, "%s/inputspec_%s.txt", All.OutputDir, All.SnapshotFileBase);

      FILE *fd = fopen(buf, "w");

      double gf = ngenic_growth_factor(0.001, 1.0) / (1.0 / 0.001);

      double DDD = ngenic_growth_factor(All.cf_atime, 1.0);

      fprintf(fd, "%12g %12g\n", All.cf_redshift, DDD); /* print actual starting redshift and
                                                           linear growth factor for this cosmology */

      double kstart = 2 * M_PI / (1000.0 * (1e6 * PARSEC / All.UnitLength_in_cm));  /* 1000 Mpc/h */
      double kend   = 2 * M_PI / (0.001 * (3.1e6 * PARSEC / All.UnitLength_in_cm)); /* 0.001 Mpc/h */

      for(double k = kstart; k < kend; k *= 1.025)
        {
          double po = ngenic_power_spec(k);
          double dl = 4.0 * M_PI * k * k * k * po;

          double kf = 0.5;

          double po2 = ngenic_power_spec(1.001 * k * kf);
          double po1 = ngenic_power_spec(k * kf);
          double dnl = 0, knl = 0;

          if(po != 0 && po1 != 0 && po2 != 0)
            {
              double neff = (log(po2) - log(po1)) / (log(1.001 * k * kf) - log(k * kf));

              if(1 + neff / 3 > 0)
                {
                  double A     = 0.482 * pow(1 + neff / 3, -0.947);
                  double B     = 0.226 * pow(1 + neff / 3, -1.778);
                  double alpha = 3.310 * pow(1 + neff / 3, -0.244);
                  double beta  = 0.862 * pow(1 + neff / 3, -0.287);
                  double V     = 11.55 * pow(1 + neff / 3, -0.423) * 1.2;

                  dnl = fnl(dl, A, B, alpha, beta, V, gf);
                  knl = k * pow(1 + dnl, 1.0 / 3);
                }
            }

          fprintf(fd, "%12g %12g    %12g %12g\n", k, dl, knl, dnl);
        }
      fclose(fd);
    }
}

#endif
