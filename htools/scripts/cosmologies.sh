#!/bin/bash 

case $1 in
    
    nu05)

        COSMO_MHDM=0.155                 # HDM mass
        COSMO_DEG=3                      # HDM (neutrino) degeneracy for class
        COSMO_NUR=0.0                    # Ultra-relativistic dofs
        COSMO_THDMK=1.94546              # HDM component temperature today in K
        COSMO_FHDM=psd/f_fd0.dat         # HDM distribution function to pass for GLR
        COSMO_NS=0.9665                  # Primordial spectral index n_s
        COSMO_SIG8=0.711                 # Sigma8
        COSMO_H0=0.6766                  # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.3111                  # Omega_m for all matter
        COSMO_OB=0.04897                 # Omega_b for baryons
        COSMO_OH=0.0109                  # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    
    ;;

    ax01)
        
	COSMO_MHDM=0.57                  # HDM mass
        COSMO_DEG=1                      # HDM degeneracy for class
        COSMO_NUR=3.044                  # Ultra-relativistic dofs
        COSMO_THDMK=1.86                 # HDM component temperature today in K
        COSMO_FHDM=psd/f_ax01.dat        # HDM distribution function to pass for GLR
        COSMO_NS=0.9665                  # Primordial spectral index n_s
        COSMO_SIG8=0.752                 # Sigma8
        COSMO_H0=0.6766                  # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.3111                  # Omega_m for all matter
        COSMO_OB=0.04897                 # Omega_b for baryons
        COSMO_OH=0.00524                 # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    
    ;;

    ax02)
        
	COSMO_MHDM=0.29                  # HDM mass
        COSMO_DEG=1                      # HDM degeneracy for class
        COSMO_NUR=3.044                  # Ultra-relativistic dofs
        COSMO_THDMK=1.86                 # HDM component temperature today in K
        COSMO_FHDM=psd/f_ax02.dat        # HDM distribution function to pass for GLR
        COSMO_NS=0.9665                  # Primordial spectral index n_s
        COSMO_SIG8=0.791                 # Sigma8
        COSMO_H0=0.6766                  # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.3111                  # Omega_m for all matter
        COSMO_OB=0.04897                 # Omega_b for baryons
        COSMO_OH=0.00178                 # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    
    ;;

    ax03)
        
	COSMO_MHDM=0.19                  # HDM mass
        COSMO_DEG=1                      # HDM degeneracy for class
        COSMO_NUR=3.044                  # Ultra-relativistic dofs
        COSMO_THDMK=1.86                 # HDM component temperature today in K
        COSMO_FHDM=psd/f_ax03.dat        # HDM distribution function to pass for GLR
        COSMO_NS=0.9665                  # Primordial spectral index n_s
        COSMO_SIG8=0.8                   # Sigma8
        COSMO_H0=0.6766                  # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.3111                  # Omega_m for all matter
        COSMO_OB=0.04897                 # Omega_b for baryons
        COSMO_OH=0.000699                # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    ;;

esac


