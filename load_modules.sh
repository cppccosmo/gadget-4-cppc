#!/bin/bash

# The Intel-MPI option seems to be the best at this 
# point,since the --enable mpi of hdf5 and fftw3 have 
# been built with that one. OpenMPI might have
# some library issues, although this has to be checked.
# Alternative: module load openmpi/4.1.4

case $1 in 
    gadi)
        module purge
        module load gcc gsl hdf5 fftw3 intel-mpi intel-compiler python3
    ;;
    katana)
        module purge 
        module load hdf5/1.12.2
        module load fftw/3.3.10
        module load gcc
        module load gsl
        module load intel-mpi/2021.7.1
        module load intel-compilers/2022.2.1 
        module load python/3.10.8
    ;;
esac

echo ''
echo "Using following modules: "
echo $(module list | sed 's/Currently Loaded Modulefiles://g' | sed 's/Key://g' | sed 's/(symbolic-version)//g' | awk '{for(i=1; i      <=NF; i++) if ($i !~ /^[0-9]+\)/) print $i}')
echo ''
