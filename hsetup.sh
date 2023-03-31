#!/bin/bash

set -e 
cdir=$(pwd)
export MuFLRPath=$cdir/MF 

########################################################

# Modify your path here 
export ClassPath=/home/z5278074/class_public

# Load modules
source load_modules.sh 

########################################################

if [[ ! -d "$ClassPath" ]]; then
    echo 'Cannot find CLASS in the path you specified!'
    exit 1
fi

if [[ -d $2 ]];
    then echo 'Directory' $2 'already exists!'
    exit 1
else
    mkdir $2
fi

case $1 in 
    LR)
        source htools/scripts/build_LR.sh $2
    ;;
    MF)
        source htools/scripts/build_MF.sh $2
    ;;
    GLR)
        source htools/scripts/build_GLR.sh $2
    ;;
    HYB)
        source htools/scripts/build_HYB.sh $2
    ;;
esac

