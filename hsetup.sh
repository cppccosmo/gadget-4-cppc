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

NTAU1=$(grep "const int N_tau = " MF/AU_cosmofunc.h | sed -e s/"const int N_tau = "//g -e s/";.*$"//g)
NMU1=$(grep "const int N_mu = " MF/AU_cosmofunc.h | sed -e s/"const int N_mu = "//g -e s/";.*$"//g)
NTAU2=$(grep "N_TAU" Config.sh | awk -F'=' '{print $2}')
NMU2=$(grep "N_MU" Config.sh | awk -F'=' '{print $2}')

if [[ $NTAU1 != $NTAU2 ]];
then echo 'MuFLR has N_tau =' $NTAU1 ', while gadget-4 has N_tau =' $NTAU2
 echo 'Set the parameters equally and recompile! (Either in MF/AU_cosmofunc.h or Config.sh)'
 rm -r $1
 exit 1
fi

if [[ $NMU1 != $NMU2 ]];
then echo 'MuFLR has N_mu =' $NMU1 ', while gadget-4 has N_mu =' $NMU2
 echo 'Set the parameters equally and recompile! (Either in MF/AU_cosmofunc.h or Config.sh)'
 rm -r $1
 exit 1
fi

while true; do
  read -p 'Do you want to restart a Hybrid simulation? [y/n] ' RES_HYB
  case $RES_HYB in
    [yY]) break;;
    [nN]) break;;
    *) echo "Please enter only 'y' or 'n'.";;
  esac
done

if [[ $RES_HYB == "n" ]];
then 
    while true; do
      read -p 'Do you want to restart a MultiFluid simulation? [y/n] ' RES_MF
      case $RES_MF in
        [yY]) break;;
        [nN]) break;;
        *) echo "Please enter only 'y' or 'n'.";;
      esac
    done
fi 

case $RES_HYB in
    y)
        source htools/scripts/build_LR.sh $1 $2 
        source htools/scripts/build_MF.sh $1
        source htools/scripts/build_HY.sh $1
        echo ''
        echo 'Ready to start in folder ' $1
        sleep 1
    ;;
    n)
        case $RES_MF in 
            y)
                echo 'Preparing LR and MF'
                source htools/scripts/build_LR.sh $1 $2 
                source htools/scripts/build_MF.sh $1
                echo ''
                echo 'Ready to start in folder ' $1
            ;;
            n) 
                echo 'Preparing LR'
                source htools/scripts/build_LR.sh $1 $2 
                echo ''
                echo 'Ready to start in folder ' $1
            ;;
        esac
    ;;
esac


