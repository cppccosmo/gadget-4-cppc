#!/bin/bash

set -e

DERIV_COSMO_OC=$(echo "${COSMO_OM} - ${COSMO_OB} - ${COSMO_OH}" | bc -l)
DERIV_COSMO_THDM_TCMB=$(echo "${COSMO_THDMK} / ${COSMO_TCMBK}" | bc -l)
DERIV_NZ=$(echo ${OUTPUT_Z} | wc -w)

case $1 in 
    0)
    
    CLASS_TEMPLATE=class_template_LR.ini
    echo 'Selecting the low-res class template'
    ;;
    1)

    CLASS_TEMPLATE=class_template_HR.ini
    echo 'Selecting the high-res class template'
    ;;
esac



echo 'Adjusting class_template.ini with the cosmology selected'
sed -e s/TEMPLATE_CLASS_H0/${COSMO_H0}/g \
	-e s/TEMPLATE_CLASS_TAU/${COSMO_TAU}/g \
	-e s/TEMPLATE_CLASS_TCMBK/${COSMO_TCMBK}/g \
	-e s/TEMPLATE_CLASS_OB/${COSMO_OB}/g \
	-e s/TEMPLATE_CLASS_OC/${DERIV_COSMO_OC}/g \
	-e s+TEMPLATE_CLASS_F_DISTRIB+"${COSMO_FHDM}"+g \
	-e s/TEMPLATE_CLASS_OH/${COSMO_OH}/g \
	-e s/TEMPLATE_CLASS_THDM_TCMB/${DERIV_COSMO_THDM_TCMB}/g \
	-e s/TEMPLATE_CLASS_DEG_HDM/${COSMO_GHDM}/g \
	-e s/TEMPLATE_CLASS_W0/${COSMO_W0}/g \
	-e s/TEMPLATE_CLASS_WA/${COSMO_WA}/g \
	-e s/TEMPLATE_CLASS_SIG8/${COSMO_SIG8}/g \
	-e s/TEMPLATE_CLASS_NS/${COSMO_NS}/g \
	${CLASS_TEMPLATE} > input_class_${COSMO_MODEL}.ini
echo 'Created file' class_${COSMO_MODEL}.ini
echo ''
sleep 1

echo 'Creating the input file for multifluid with the cosmology selected'
echo ${COSMO_NS} > params_MuFLR-HDM.dat
echo ${COSMO_SIG8} >> params_MuFLR-HDM.dat
echo ${COSMO_H0} >> params_MuFLR-HDM.dat
echo ${COSMO_OM} >> params_MuFLR-HDM.dat
echo ${COSMO_OB} >> params_MuFLR-HDM.dat
echo ${COSMO_OH} >> params_MuFLR-HDM.dat
echo ${COSMO_TCMBK} >> params_MuFLR-HDM.dat
echo ${COSMO_W0} >> params_MuFLR-HDM.dat
echo ${COSMO_WA} >> params_MuFLR-HDM.dat
echo ${OUTPUT_NL} >> params_MuFLR-HDM.dat
echo "0" >> params_MuFLR-HDM.dat
echo ${OUTPUT_PRINT} >> params_MuFLR-HDM.dat
echo "0" >> params_MuFLR-HDM.dat
echo "200" >> params_MuFLR-HDM.dat
echo ${DERIV_NZ} >> params_MuFLR-HDM.dat
echo "${OUTPUT_Z}" >> params_MuFLR-HDM.dat
echo "${CLASS_ROOT}/out_tk.dat" >> params_MuFLR-HDM.dat
echo "1" >> params_MuFLR-HDM.dat
echo ${COSMO_THDMK} >> params_MuFLR-HDM.dat
echo ${COSMO_GHDM} >> params_MuFLR-HDM.dat
echo ${COSMO_FHDM} >> params_MuFLR-HDM.dat

#N_TAU=$(grep "const int N_tau = " AU_cosmofunc.h | sed -e s/"const int N_tau = "//g -e s/";.*$"//g)
#N_MU=$(grep "const int N_mu = " AU_cosmofunc.h | sed -e s/"const int N_mu = "//g -e s/";.*$"//g)

