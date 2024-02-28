#!/bin/bash 

case $1 in
    
    nu01)

        COSMO_MHDM=0.305519              # Nu/HDM mass
        COSMO_THDMK=1.95280830923360     # Nu/HDM component temperature today in K
        COSMO_FHDM=psd/f_fd0.dat         # Nu/HDM distribution function to pass for GLR
        COSMO_NS=0.963                   # Primordial spectral index n_s
        COSMO_SIG8=0.8                   # Sigma8
        COSMO_H0=0.71                    # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.2648284070620909      # Omega_m for all matter
        COSMO_OB=0.044792699861138666    # Omega_b for baryons
        COSMO_OH=0.019837333862328905    # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state 
        COSMO_WA=0.0                     
    
    ;;
    
    nu05)

        COSMO_MHDM=0.052135              # HDM mass
        COSMO_THDMK=1.95280830923360     # HDM component temperature today in K
        COSMO_FHDM=psd/f_fd0.dat         # HDM distribution function to pass for GLR
        COSMO_NS=0.96                    # Primordial spectral index n_s
        COSMO_SIG8=0.79                  # Sigma8
        COSMO_H0=0.67                    # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.3190020049            # Omega_m for all matter
        COSMO_OB=0.04900868              # Omega_b for baryons
        COSMO_OH=0.0038093               # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    
    ;;

    be01)
        
	COSMO_MHDM=0.93                  # HDM mass
        COSMO_THDMK=1.95280830923360     # HDM component temperature today in K
        COSMO_FHDM=psd/f_be0.dat         # HDM distribution function to pass for GLR
        COSMO_NS=0.963                   # Primordial spectral index n_s
        COSMO_SIG8=0.8                   # Sigma8
        COSMO_H0=0.71                    # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.2648284070620909      # Omega_m for all matter
        COSMO_OB=0.044792699861138666    # Omega_b for baryons
        COSMO_OH=0.0131995               # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    
    ;; 
    
    be05)
        
	COSMO_MHDM=0.15                  # HDM mass
        COSMO_THDMK=1.95280830923360     # HDM component temperature today in K
        COSMO_FHDM=psd/f_be0.dat         # HDM distribution function to pass for GLR
        COSMO_NS=0.963                   # Primordial spectral index n_s
        COSMO_SIG8=0.8                   # Sigma8
        COSMO_H0=0.71                    # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.2648284070620909      # Omega_m for all matter
        COSMO_OB=0.044792699861138666    # Omega_b for baryons
        COSMO_OH=0.00212897              # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    
    ;;
    
    ax01)
        
	COSMO_MHDM=0.57                  # HDM mass
        COSMO_THDMK=1.7034375            # HDM component temperature today in K
        COSMO_FHDM=psd/f_ax01.dat        # HDM distribution function to pass for GLR
        COSMO_NS=0.963                   # Primordial spectral index n_s
        COSMO_SIG8=0.8                   # Sigma8
        COSMO_H0=0.71                    # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.2648284070620909      # Omega_m for all matter
        COSMO_OB=0.044792699861138666    # Omega_b for baryons
        COSMO_OH=0.0037                  # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    
    ;;

    ax02)
        
	COSMO_MHDM=0.29                  # HDM mass
        COSMO_THDMK=1.8724185            # HDM component temperature today in K
        COSMO_FHDM=psd/f_ax02.dat        # HDM distribution function to pass for GLR
        COSMO_NS=0.963                   # Primordial spectral index n_s
        COSMO_SIG8=0.8                   # Sigma8
        COSMO_H0=0.71                    # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.2648284070620909      # Omega_m for all matter
        COSMO_OB=0.044792699861138666    # Omega_b for baryons
        COSMO_OH=0.0016                  # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    
    ;;

    ax03)
        
	COSMO_MHDM=0.19                  # HDM mass
        COSMO_THDMK=1.8724185            # HDM component temperature today in K
        COSMO_FHDM=psd/f_ax03.dat        # HDM distribution function to pass for GLR
        COSMO_NS=0.963                   # Primordial spectral index n_s
        COSMO_SIG8=0.8                   # Sigma8
        COSMO_H0=0.71                    # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.2648284070620909      # Omega_m for all matter
        COSMO_OB=0.044792699861138666    # Omega_b for baryons
        COSMO_OH=0.0006                  # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    


esac


