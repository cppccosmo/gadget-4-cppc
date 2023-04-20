#!/bin/bash

set -e 
cdir=$(pwd)
export MuFLRPath=$cdir/MF 

########################################################

# Modify your path here 
export ClassPath=~/class_public

# Load modules
source load_modules.sh $2

########################################################

if [[ ! -d "$ClassPath" ]]; then
    echo 'Cannot find CLASS in the path you specified!'
    exit 1
fi

if [[ -d $1 ]];
    then echo 'Directory' $1 'already exists!'
    exit 1
else
    mkdir $1
fi

while true; do
  read -p 'Do you want to restart a MultiFluid simulation? [y/n] ' RES_MF
  case $RES_MF in
    [yY]) break;;
    [nN]) break;;
    *) echo "Please enter only 'y' or 'n'.";;
  esac
done

while true; do
  read -p 'Do you want to restart a Hybrid simulation? [y/n] ' RES_HYB
  case $RES_HYB in
    [yY]) break;;
    [nN]) break;;
    *) echo "Please enter only 'y' or 'n'.";;
  esac
done

case $RES_HYB in
    y)
        case $RES_MF in 
            y)
                echo 'Preparing LR, MF and Hybrid'
                source htools/scripts/build_full.sh # Hybrid restart and MF restart
            ;;
            n) 
                echo 'Preparing LR and Hybrid'
                source htools/scripts/build_HYB.sh # Hybrid restart only
            ;;
        esac
    ;;
    n)
        case $RES_MF in 
            y)
                echo 'Preparing LR and MF'
                source htools/scripts/build_MF.sh # MF restart only
            ;;
            n) 
                echo 'Preparing LR'
                source htools/scripts/build_LR.sh # LR only
            ;;
        esac
    ;;
esac


