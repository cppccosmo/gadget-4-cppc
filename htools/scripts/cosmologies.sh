#!/bin/bash 

case $1 in
    
    nu05)

        COSMO_MHDM=0.155                 # HDM mass
        COSMO_DEG=3                      # HDM (neutrino) degeneracy for class
        COSMO_NUR=0.0                    # Ultra-relativistic dofs
        COSMO_THDMK=1.94546              # HDM component temperature today in K
        COSMO_FHDM=psd/f_fd0.dat         # HDM distribution function to pass for GLR
        COSMO_AS=2.100549e-09            # Primordial amplitude at pivot scale
        COSMO_NS=0.9660499               # Primordial spectral index n_s
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
        COSMO_FHDM=psd/f_ax01s.dat        # HDM distribution function to pass for GLR
        COSMO_AS=2.100549e-09            # Primordial amplitude at pivot scale
        COSMO_NS=0.9660499               # Primordial spectral index n_s
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
        COSMO_FHDM=psd/f_ax02s.dat        # HDM distribution function to pass for GLR
        COSMO_AS=2.100549e-09            # Primordial amplitude at pivot scale
        COSMO_NS=0.9660499               # Primordial spectral index n_s
        COSMO_SIG8=0.791                 # Sigma8
        COSMO_H0=0.6766                  # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.3111                  # Omega_m for all matter
        COSMO_OB=0.04897                 # Omega_b for baryons
        COSMO_OH=0.00178                 # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    
    ;;

    e1)
        
	COSMO_MHDM=0.04760               # HDM mass
        COSMO_DEG=1                      # HDM degeneracy for class
        COSMO_NUR=1.0146666              # Ultra-relativistic dofs
        COSMO_THDMK=1.87893              # HDM component temperature today in K
        COSMO_FHDM=psd/f_e1.dat          # HDM distribution function to pass for GLR
        COSMO_AS=2.100549e-09            # Primordial amplitude at pivot scale
        COSMO_NS=0.9660499               # Primordial spectral index n_s
        COSMO_SIG8=0.7876                # Sigma8
        COSMO_H0=0.6766                  # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.3111                  # Omega_m for all matter
        COSMO_OB=0.04897                 # Omega_b for baryons
        COSMO_OH=0.00246021              # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    ;;

    e2)
        
	COSMO_MHDM=0.053579              # HDM mass
        COSMO_DEG=1                      # HDM degeneracy for class
        COSMO_NUR=0.0                    # Ultra-relativistic dofs
        COSMO_THDMK=1.952463             # HDM component temperature today in K
        COSMO_FHDM=psd/f_e2.dat          # HDM distribution function to pass for GLR
        COSMO_AS=2.100549e-09            # Primordial amplitude at pivot scale
        COSMO_NS=0.9660499               # Primordial spectral index n_s
        COSMO_SIG8=0.786004              # Sigma8
        COSMO_H0=0.6766                  # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.3111                  # Omega_m for all matter
        COSMO_OB=0.04897                 # Omega_b for baryons
        COSMO_OH=0.00377395              # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    ;;

    e3)
        
	COSMO_MHDM=0.07382               # HDM mass
        COSMO_DEG=1                      # HDM degeneracy for class
        COSMO_NUR=1.0146666              # Ultra-relativistic dofs
        COSMO_THDMK=1.87893              # HDM component temperature today in K
        COSMO_FHDM=psd/f_e3.dat          # HDM distribution function to pass for GLR
        COSMO_AS=2.100549e-09            # Primordial amplitude at pivot scale
        COSMO_NS=0.9660499               # Primordial spectral index n_s
        COSMO_SIG8=0.751028              # Sigma8
        COSMO_H0=0.6766                  # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.3111                  # Omega_m for all matter
        COSMO_OB=0.04897                 # Omega_b for baryons
        COSMO_OH=0.00446593              # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    ;;

    e4)
        
	COSMO_MHDM=0.105075              # HDM mass
        COSMO_DEG=1                      # HDM degeneracy for class
        COSMO_NUR=0.0                    # Ultra-relativistic dofs
        COSMO_THDMK=1.952463             # HDM component temperature today in K
        COSMO_FHDM=psd/f_e4.dat          # HDM distribution function to pass for GLR
        COSMO_AS=2.100549e-09            # Primordial amplitude at pivot scale
        COSMO_NS=0.9660499               # Primordial spectral index n_s
        COSMO_SIG8=0.747084              # Sigma8
        COSMO_H0=0.6766                  # Hubble paramter today in units of 100 km/s/Mpc
        COSMO_OM=0.3111                  # Omega_m for all matter
        COSMO_OB=0.04897                 # Omega_b for baryons
        COSMO_OH=0.00740086              # Omega_h for the HDM component
        COSMO_TCMBK=2.7255               # CMB temperature today in K
        COSMO_W0=-1.0                    # Dark energy equation of state
        COSMO_WA=0.0                            
    ;;


esac


