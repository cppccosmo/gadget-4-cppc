

GADGET-4-CPPC
========

![](documentation/cover.png)

GADGET-4-CPPC is an N-body code to study neutrinos and any hot dark matter 
with linear response, multifluid and hybrid methods. It is based on the code 
GADGET-4, that is a massively parallel code for N-body/hydrodynamical
cosmological simulations and mainly written by Volker Springel.

The neutrino hybrid method implementations have been performed by Joe Chen, and a previous version of the 
code (for neutrinos only) can be found in [gadget4-hybrid_public](https://github.com/joechenUNSW/gadget4-hybrid_public).

The hot dark matter multifluid part integrated here has been developed by Amol Upadhye. The original repository can be found in [MuFLR-HDM](https://github.com/upadhye/MuFLR-HDM).



Documentation
=============

A documentation on the GADGET-4-CPPC additions and methods is under development.

For reference on the physics implemented, please refer to the main papers:

- *Joe Zhiyu Chen, Markus R. Mosbech, Amol Upadhye, Yvonne Y.Y. Wong*, Hybrid multi-fluid-particle simulations of the cosmic neutrino background,     [2210.16012](https://arxiv.org/abs/2210.16012) 
- *Joe Zhiyu Chen, Amol Upadhye, Yvonne Y.Y. Wong*, The cosmic neutrino background as a collection of fluids in large-scale structure simulations,     [2011.12503](https://arxiv.org/abs/2011.12503) 
- *Joe Zhiyu Chen, Amol Upadhye, Yvonne Y.Y. Wong*, One line to run them all: SuperEasy massive neutrino linear response in N-body simulations, [    2011.12504](https://arxiv.org/abs/2011.12504) 



Installation and Dependencies 
=============

To compile the code and run a typical simulation you will need the following libraries:

- Working C/C++ compilers (for GNU version 4.x or later)
- GNU Scientific Library ([GSL](https://www.gnu.org/software/gsl/))
- Fastest Fourier Transform in the West ([FFTW3](http://www.fftw.org/)) 
- Hierarchical Data Format ([HDF5](https://www.hdfgroup.org/solutions/hdf5/))
- Message Passing Interface (MPI), version 3.0 or higher. Examples are [Intel-MPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html), [OpenMPI](https://www.open-mpi.org/), [MPICH](https://www.mpich.org/) 

To install the code clone the repository 
```
git clone https://github.com/cppccosmo/gadget-4-cppc.git; cd gadget-4-cppc
```
Then to start a simulation folder choose between Linear Response (LR), Multifluid (MF) or Hybrid (HYB) and run  
```
bash hsetup.sh /path/to/sim/folder <system>
```
An example is
```
bash hsetup.sh /tests/hybrid_sim gadi
```


