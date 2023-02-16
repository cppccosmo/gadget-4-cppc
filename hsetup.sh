#!/bin/bash

# Modify your path here 
ClassPath='/home/z5278074/class_public'
MuFLRPath='/home/z5278074/MuFLR'

# Load your modules 

# Katana

module load hdf5/1.12.2
module load fftw/3.3.10
module load gcc
module load gsl
module load openmpi/4.1.4
module load intel-compilers/2022.2.1 
module load python/3.10.8

#############################################################################

cdir=$(pwd)

if [[ ! -d "$ClassPath" ]]; then
    echo 'Cannot find CLASS in the path you specified!'
    exit 1
fi

if [[ ! -d "$MuFLRPath" ]]; then
    echo 'Cannot find MuFLR in the path you specified!'
    exit 1
fi

if [[ -d $1 ]];
    then echo 'Directory already exists!'
    exit 1
else
    mkdir $1
fi

echo 'Checking the N_tau and N_mu parameters ...'
sleep 1

Ntau1=$(awk -F " |;" '{if ($3 == "N_tau") print $5}' $MuFLRPath/AU_cosmofunc.h)
Nmu1=$(awk -F " |;" '{if ($3 == "N_mu") print $5}'   $MuFLRPath/AU_cosmofunc.h)
Nk=$(sed -n '24p' $MuFLRPath/AU_fftgrid.h | awk '{print $3}' | sed 's/[^[:digit:]]/ /g; s/  */ /; s/^  *//')
Neff1a=$(sed -n '24p' $MuFLRPath/AU_cosmoparam.h | awk '{print $3}' | sed 's/[^[:digit:]]/ /g; s/  */ /; s/^  *//' | awk '{print $1}')
Neff1b=$(sed -n '24p' $MuFLRPath/AU_cosmoparam.h | awk '{print $3}' | sed 's/[^[:digit:]]/ /g; s/  */ /; s/^  *//' | awk '{print $2}')
Neffm1a=$(sed -n '25p' ../MuFLR/AU_cosmoparam.h | awk '{print $3}' | sed 's/[^[:digit:]]/ /g; s/  */ /; s/^  *//' | awk '{print $1}')
Neffm1b=$(sed -n '25p' ../MuFLR/AU_cosmoparam.h | awk '{print $3}' | sed 's/[^[:digit:]]/ /g; s/  */ /; s/^  *//' | awk '{print $2}')

Neff1=$Neff1a.$Neff1b
Neffm1=$Neffm1a.$Neffm1b


# This will be changed when N_mu and N_tau are config options
Ntau2=$(sed -n '13p' ./src/neutrinos/neutrinomflr.h | awk -F " |;" '{print $10}')
Nmu2=$(sed -n '14p' ./src/neutrinos/neutrinomflr.h | awk -F " |;" '{print $10}')
Neff2=$(sed -n '11p' ./src/neutrinos/neutrinomflr.h | awk -F " |;" '{print $9}')
Neffm2=$(sed -n '12p' ./src/neutrinos/neutrinomflr.h | awk -F " |;" '{print $9}')


if [[ $Ntau1 != $Ntau2 ]];
then echo 'MuFLR has N_tau =' $Ntau1 ', while gadget-4 has N_tau =' $Ntau2
 echo 'Set the parameters equally and recompile!' 
 exit 1
fi

if [[ $Nmu1 != $Nmu2 ]];
then echo 'MuFLR has N_mu =' $Nmu1 ', while gadget-4 has N_mu =' $Nmu2
 echo 'Set the parameters equally and recompile!'
 exit 1 
fi

if [[ $Neff1 != $Neff2 ]];
then echo 'MuFLR has N_nu_eff =' $Neff1 ', while gadget-4 has N_nu_eff =' $Neff2
 echo 'Set the parameters equally and recompile!' 
 exit 1
fi

if [[ $Neffm1 != $Neffm2 ]];
then echo 'MuFLR has N_nu_massive =' $Neffm1 ', while gadget-4 has N_nu_massive =' $Neffm2
 echo 'Set the parameters equally and recompile!' 
 exit 1
fi

echo 'Preparing simulation folder in' $1
echo 'You are compiling with N_tau =' $Ntau1 ', N_mu =' $Nmu1 'and N_k =' $Nk 
read -p 'Do you wish to continue? [y/n] ' ANS

if [[ $ANS != "y" ]];
    then exit 1
fi

echo 'Copying htools in' $1 '...'
cp $ClassPath/class $1
mkdir $1/scripts
mkdir $1/start_files
cp htools/htools $1; chmod u+x $1/htools
cp htools/nu.ini $1
cp htools/param.txt $1
cp htools/param_restart.txt $1
cp htools/scripts/* $1/scripts/
cp htools/job $1

echo 'Compiling MuFLR ...'

cd $MuFLRPath
make
rsync -avh MuFLR $1
cd $cdir

Ntypes=$(awk -F"\t" '/NTYPES/ {print $0}' htools/Config.sh | awk -F "=" '{print $2}')
PMGrid=$(awk -F"\t" '/PMGRID/ {print $0}' htools/Config.sh | awk -F "=" '{print $2}')
Ngenic=$(awk -F"\t" '/NGENIC/ {print $0}' htools/Config.sh | awk -F "=" '{print $2}')

NtypesM=$(awk -F"\t" '/NTYPES/ {print $0}' htools/Config_mf_restart.sh | awk -F "=" '{print $2}')
PMGridM=$(awk -F"\t" '/PMGRID/ {print $0}' htools/Config_mf_restart.sh | awk -F "=" '{print $2}')

NtypesR=$(awk -F"\t" '/NTYPES/ {print $0}' htools/Config_hybrid_restart.sh | awk -F "=" '{print $2}')
PMGridR=$(awk -F"\t" '/PMGRID/ {print $0}' htools/Config_hybrid_restart.sh | awk -F "=" '{print $2}')
NgenicR=$(awk -F"\t" '/NGENIC/ {print $0}' htools/Config_hybrid_restart.sh | awk -F "=" '{print $2}')

if [[ $Ntypes != $NtypesR ]];
    then echo 'NTYPES has to match in Config.sh, Config_mf_restart.sh and Config_hybrid_restart.sh '
    sleep 1 
    exit 1
fi

if [[ $Ntypes != $NtypesM ]];
    then echo 'NTYPES has to match in Config.sh, Config_mf_restart.sh and Config_hybrid_restart.sh '
    sleep 1 
    exit 1
fi

if [[ $NtypesM != $NtypesR ]];
    then echo 'NTYPES has to match in Config.sh, Config_mf_restart.sh and Config_hybrid_restart.sh '
    sleep 1 
    exit 1
fi 

if [[ $PMGrid != $PMGridR ]];
    then echo 'PMGRID has to match in Config.sh, Config_mf_restart and Config_hybrid_restart.sh '
    sleep 1 
    exit 1
fi

if [[ $PMGrid != $PMGridM ]];
    then echo 'PMGRID has to match in Config.sh, Config_mf_restart and Config_hybrid_restart.sh '
    sleep 1 
    exit 1
fi

if [[ $PMGridM != $PMGridR ]];
    then echo 'PMGRID has to match in Config.sh, Config_mf_restart and Config_hybrid_restart.sh '
    sleep 1 
    exit 1
fi

echo 'Compiling Hybrid for start (Config.sh) with'
echo ''
echo 'NTYPES =' $Ntypes  
echo 'PMGRID =' $PMGrid 
echo 'NGENIC =' $Ngenic
echo ''

sleep 1

make clean
make
rsync -avh Gadget4 $1

echo 'Compiling Hybrid for mf restart (Config_mf_restart.sh) with'
echo ''
echo 'NTYPES =' $NtypesM 
echo 'PMGRID =' $PMGridM
echo ''

sleep 1 

make clean
make CONFIG=htools/Config_mf_restart.sh EXEC=$1/Gadget4-mf-restart


echo 'Compiling Hybrid for hybrid restart (Config_hybrid_restart.sh) with'
echo ''
echo 'NTYPES =' $NtypesR  
echo 'PMGRID =' $PMGridR 
echo 'NGENIC =' $NgenicR
echo ''

sleep 1

make clean
make DIR=$1

echo 'Writing the log.txt ...'
sleep 1

echo 'Compilation options'  > $1/log.txt
echo ''                    >> $1/log.txt
echo 'N_tau =' $Ntau1      >> $1/log.txt
echo 'N_mu =' $Nmu1        >> $1/log.txt
echo 'N_k =' $Nk           >> $1/log.txt
echo ''                    >> $1/log.txt
echo 'Config.sh had'       >> $1/log.txt
echo ''                    >> $1/log.txt

echo $(sed '/^#/d' htools/Config.sh)                >> $1/log.txt 
echo ''                                             >> $1/log.txt
echo 'Config_mf_restart.sh had'                     >> $1/log.txt
echo ''                                             >> $1/log.txt
echo $(sed '/^#/d' htools/Config_mf_restart.sh)     >> $1/log.txt 
echo ''                                             >> $1/log.txt
echo 'Config_hybrid_restart.sh had'                 >> $1/log.txt
echo ''                                             >> $1/log.txt
echo $(sed '/^#/d' htools/Config_hybrid_restart.sh) >> $1/log.txt 

echo 'Ready to start in folder ' $1
sleep 1



