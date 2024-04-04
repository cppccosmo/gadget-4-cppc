#!/bin/bash

set -e 

NTAU1=$(grep "N_TAU" Config.sh | awk -F'=' '{print $2}')
NMU1=$(grep "N_MU" Config.sh | awk -F'=' '{print $2}')
NTAU2=$(grep "N_TAU" Config_hybrid_restart.sh | awk -F'=' '{print $2}')
NMU2=$(grep "N_MU" Config_hybrid_restart.sh | awk -F'=' '{print $2}')

if [[ $NTAU1 != $NTAU2 ]]; then 
    echo ''
    echo 'ERROR: Config.sh has N_tau =' $NTAU1 ', while Config_hybrid_restart has N_tau =' $NTAU2
    echo 'Set the parameters equally and recompile!'
    rm -r $1
    exit 1
fi

if [[ $NMU1 != $NMU2 ]]; then
    echo ''
    echo 'ERROR: Config.sh has N_tau =' $NMU1 ', while Config_hybrid_restart has N_tau =' $NMU2
    echo 'Set the parameters equally and recompile!'
    rm -r $1
    exit 1
fi

cp htools/param_restart_template.txt $1

echo 'Compiling Gadget4 hybrid ...'
make -j CONFIG=Config_hybrid_restart.sh EXEC=$1/Gadget4-hybrid-restart BUILD_DIR=$1/build_hy >> $1/info/compilation.txt 2> $1/info/issues.txt 

sleep 1

