#!/bin/bash 

case $1 in
    
    nu01)

        COSMO_MHDM=0.305519              # Nu/HDM mass
        COSMO_THDMK=1.95280830923360     # Nu/HDM component temperature today in K
        COSMO_FHDM=psd/f_fd0.dat         # Nu/HDM distribution function to pass for GLR
        COSMO_NS=0.9665                  # Primordial spectral index n_s
        COSMO_SIG8=0.7139                # Sigma8
        COSMO_H0=0.6766                  # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.2648284070620909      # Omega_m for all matter
        COSMO_OB=0.044792699861138666    # Omega_b for baryons
        COSMO_OH=0.019837333862328905    # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state 
        COSMO_WA=0.0                     
    ;;

    nu02)
        
    ;;
    
    nu03)
        
    ;;
    
    nu04)
        
    ;;
    
    nu05)
        
    ;;
esac


