#!/bin/bash

set -e 

cdir=$(pwd)

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
cp htools/templates/param_template.txt $1
cp htools/templates/class_template.ini $1 
cp MultiFluid/distribution_functions/* $1/psd

case $2 in
    katana)
	cp htools/jobs/katana-run $1
	cp htools/jobs/katana-tools $1
    ;;
    gadi)
	cp htools/jobs/gadi-run $1
	cp htools/jobs/gadi-tools $1
    ;;
esac

echo 'Compiling MuFLR ...'

cd $MuFLRPath
make > $1/info/compilation.txt
cp MultiFluid $1; cd $cdir
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
echo ''

sleep 1

echo 'Building Gadget4 ...'
make -j2 CONFIG=Config.sh EXEC=$1/Gadget4 BUILD_DIR=$1/build >> $1/info/compilation.txt 2> $1/info/issues.txt

sleep 1

NK=$(grep "#define NK " MultiFluid/utils/fftgrid.h | awk '{print $3}' | sed 's/(//' | sed 's/)//')
N_TAU=$(grep "const int N_tau = " MultiFluid/utils/cosmofunc.h | sed -e s/"const int N_tau = "//g -e s/";.*$"//g)
N_MU=$(grep "const int N_mu = " MultiFluid/utils/cosmofunc.h | sed -e s/"const int N_mu = "//g -e s/";.*$"//g)
echo 'Compilation options'     > $1/info/log.txt
echo ''                       >> $1/info/log.txt
echo 'NK    = ' $NK           >> $1/info/log.txt
echo 'N_TAU = ' $N_TAU        >> $1/info/log.txt
echo 'N_MU  = ' $N_MU         >> $1/info/log.txt
echo ''                       >> $1/info/log.txt
echo 'Config.sh had'          >> $1/info/log.txt
echo ''                       >> $1/info/log.txt
echo $(sed '/^#/d' Config.sh) >> $1/info/log.txt 
echo ''                       >> $1/info/log.txt

