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

# Move scripts on the simulation folder
cp htools/htools $1; chmod u+x $1/htools
cp htools/scripts/*.py $1/scripts/
cp htools/scripts/cosmofiles.sh $1/scripts/
cp htools/scripts/katana_modules.sh $1/scripts/
cp htools/gadget-pre $1
cp htools/gadget-run $1
cp htools/param.txt $1
cp htools/class_template.ini $1 
cp htools/class_template_LR.ini $1 
cp htools/class_template_HR.ini $1 
cp htools/pk_ref.pre $1 
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
echo ''

sleep 1

echo 'Compiling Gadget4 ...'
make CONFIG=Config.sh EXEC=$1/Gadget4 BUILD_DIR=$1/build >> $1/info/compilation.txt 2> $1/info/issues.txt 
echo 'Done!'

sleep 1

echo 'Writing the log.txt ...'

sleep 1

NK=$(grep "#define NK " MF/AU_fftgrid.h | awk '{print $3}' | sed 's/(//' | sed 's/)//')
echo 'Compilation options'     > $1/info/log.txt
echo ''                       >> $1/info/log.txt
echo 'NK = ' $NK              >> $1/info/log.txt
echo 'Config.sh had'          >> $1/info/log.txt
echo ''                       >> $1/info/log.txt
echo $(sed '/^#/d' Config.sh) >> $1/info/log.txt 
echo ''                       >> $1/info/log.txt

echo 'Ready to start in folder ' $1
sleep 1

