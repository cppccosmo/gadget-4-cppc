#!/bin/bash

set -e 

cdir=$(pwd)

NTAU1=$(grep "const int N_tau = " MF/AU_cosmofunc.h | sed -e s/"const int N_tau = "//g -e s/";.*$"//g)
NMU1=$(grep "const int N_mu = " MF/AU_cosmofunc.h | sed -e s/"const int N_mu = "//g -e s/";.*$"//g)
NTAU2=$(grep "N_TAU" Config_mf_restart.sh | awk -F'=' '{print $2}')
NMU2=$(grep "N_MU" Config_mf_restart.sh | awk -F'=' '{print $2}')

if [[ $NTAU1 != $NTAU2 ]];
then echo 'MuFLR has N_tau =' $NTAU1 ', while gadget-4 has N_tau =' $NTAU2
 echo 'Set the parameters equally and recompile!'
 exit 1
fi

if [[ $NMU1 != $NMU2 ]];
then echo 'MuFLR has N_mu =' $NMU1 ', while gadget-4 has N_mu =' $NMU2
 echo 'Set the parameters equally and recompile!'
 exit 1
fi


echo 'Copying htools in' $1 
cp $ClassPath/class $1

mkdir $1/scripts
mkdir $1/info
mkdir $1/class_data
mkdir $1/mf_data
mkdir $1/psd

cp htools/htools $1; chmod u+x $1/htools
cp htools/scripts/*.py $1/scripts/
cp htools/scripts/cosmofiles.sh $1/scripts/
cp htools/scripts/cosmologies.sh $1/scripts/
cp load_modules.sh $1/scripts
cp htools/gadget-tools $1
cp htools/gadget-run $1
cp htools/param_template.txt $1
cp htools/class_template_LR.ini $1 
cp htools/class_template_HR.ini $1 
cp MF/distribution_functions/* $1/psd

echo 'Compiling MuFLR ...'

cd $MuFLRPath
make > $1/info/compilation.txt
cp MuFLR $1; cd $cdir
echo 'Done!'
echo ''

Ntypes=$(awk -F"\t" '/NTYPES/ {print $0}' Config.sh | awk -F "=" '{print $2}')
PMGrid=$(awk -F"\t" '/PMGRID/ {print $0}' Config.sh | awk -F "=" '{print $2}')
Ngenic=$(awk -F"\t" '/NGENIC/ {print $0}' Config.sh | awk -F "=" '{print $2}')

echo 'Using Config.sh with'
echo ''
echo 'NTYPES =' $Ntypes  
echo 'PMGRID =' $PMGrid 
echo 'NGENIC =' $Ngenic
echo 'N_TAU  =' $NTAU1
echo 'N_MU   =' $NMU1
echo ''

sleep 1

echo 'Compiling Gadget4 ...'
make CONFIG=Config.sh EXEC=$1/Gadget4 BUILD_DIR=$1/build >> $1/info/compilation.txt 2> $1/info/issues.txt 
echo 'Done!'

sleep 1

echo 'Compiling Gadget4 for MultiFluid restart...'
make CONFIG=Config_mf_restart.sh EXEC=$1/Gadget4_mf_restart BUILD_DIR=$1/build >> $1/info/compilation.txt 2> $1/info/issues.txt 
echo 'Done!'

echo 'Writing the log.txt ...'

sleep 1

NK=$(grep "#define NK " MF/AU_fftgrid.h | awk '{print $3}' | sed 's/(//' | sed 's/)//')

echo 'Compilation options'     > $1/info/log.txt
echo ''                       >> $1/info/log.txt
echo 'NK    = ' $NK           >> $1/info/log.txt
echo 'N_TAU = ' $NTAU2        >> $1/info/log.txt
echo 'N_MU  = ' $NMU2         >> $1/info/log.txt
echo ''                       >> $1/info/log.txt
echo 'Config.sh had'          >> $1/info/log.txt
echo ''                       >> $1/info/log.txt
echo $(sed '/^#/d' Config.sh) >> $1/info/log.txt 
echo ''                       >> $1/info/log.txt

echo 'Ready to start in folder ' $1
sleep 1

