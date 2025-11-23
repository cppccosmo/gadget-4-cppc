#!/bin/bash

set -e 

NTAU1=$(grep "N_TAU" Config.sh | awk -F'=' '{print $2}')
NMU1=$(grep "N_MU" Config.sh | awk -F'=' '{print $2}')
NTAU2=$(grep "N_TAU" Config_mf_restart.sh | awk -F'=' '{print $2}')
NMU2=$(grep "N_MU" Config_mf_restart.sh | awk -F'=' '{print $2}')

if [[ $NTAU1 != $NTAU2 ]]; then 
    echo ''
    echo 'ERROR: Config.sh has N_tau =' $NTAU1 ', while Config_mf_restart has N_tau =' $NTAU2
    echo 'Set the parameters equally and recompile!'
    rm -r $1
    exit 1
fi

if [[ $NMU1 != $NMU2 ]]; then
    echo ''
    echo 'ERROR: Config.sh has N_tau =' $NMU1 ', while Config_mf_restart has N_tau =' $NMU2
    echo 'Set the parameters equally and recompile!'
    rm -r $1
    exit 1
fi

cp htools/templates/param_mf_template.txt $1

echo 'Building Gadget4-mf-restart ...'
make -j2 CONFIG=Config_mf_restart.sh EXEC=$1/Gadget4-mf-restart BUILD_DIR=$1/build_mf >> $1/info/compilation_mf.txt 2> $1/info/issues_mf.txt

sleep 1

