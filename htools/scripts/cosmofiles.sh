#!/bin/bash

set -e

DERIV_COSMO_OC=$(echo "${COSMO_OM} - ${COSMO_OB} - ${COSMO_OH}" | bc -l)
DERIV_COSMO_THDM_TCMB=$(echo "${COSMO_THDMK} / ${COSMO_TCMBK}" | bc -l)
DERIV_NZ=$(echo ${OUTPUT_Z} | wc -w)

CLASS_TEMPLATE=class_template.ini

echo 'Adjusting class_template.ini with the cosmology selected'
sed -e s/TEMPLATE_CLASS_H0/${COSMO_H0}/g \
	-e s/TEMPLATE_CLASS_TCMBK/${COSMO_TCMBK}/g \
	-e s/TEMPLATE_CLASS_OB/${COSMO_OB}/g \
	-e s/TEMPLATE_CLASS_OC/${DERIV_COSMO_OC}/g \
	-e s+TEMPLATE_CLASS_F_DISTRIB+"${COSMO_FHDM}"+g \
	-e s/TEMPLATE_CLASS_OH/${COSMO_OH}/g \
	-e s/TEMPLATE_CLASS_THDM_TCMB/${DERIV_COSMO_THDM_TCMB}/g \
	-e s/TEMPLATE_CLASS_M_HDM/${COSMO_MHDM}/g \
	-e s/TEMPLATE_CLASS_NUR/${COSMO_NUR}/g \
	-e s/TEMPLATE_CLASS_DEGH/${COSMO_DEG}/g \
	-e s/TEMPLATE_CLASS_W0/${COSMO_W0}/g \
	-e s/TEMPLATE_CLASS_WA/${COSMO_WA}/g \
	-e s/TEMPLATE_CLASS_NS/${COSMO_NS}/g \
	-e s/TEMPLATE_CLASS_AS/${COSMO_AS}/g \
	${CLASS_TEMPLATE} > input_class_${COSMO_MODEL}.ini
echo 'Created file' class_${COSMO_MODEL}.ini
echo ''
sleep 1

MF_FILE=params_MF.dat

echo 'Creating the input file for multifluid with the cosmology selected'
echo ${COSMO_NS}                 > $MF_FILE 
echo ${COSMO_SIG8}              >> $MF_FILE
echo ${COSMO_H0}                >> $MF_FILE 
echo ${COSMO_OM}                >> $MF_FILE 
echo ${COSMO_OB}                >> $MF_FILE 
echo ${COSMO_OH}                >> $MF_FILE 
echo ${COSMO_TCMBK}             >> $MF_FILE
echo ${COSMO_W0}                >> $MF_FILE
echo ${COSMO_WA}                >> $MF_FILE
echo ${OUTPUT_NL}               >> $MF_FILE
echo "0"                        >> $MF_FILE
echo ${OUTPUT_PRINT}            >> $MF_FILE 
echo "0"                        >> $MF_FILE 
echo "200"                      >> $MF_FILE 
echo ${DERIV_NZ}                >> $MF_FILE 
echo "${OUTPUT_Z}"              >> $MF_FILE 
echo "${CLASS_ROOT}/out_tk.dat" >> $MF_FILE 
echo "1"                        >> $MF_FILE 
echo ${COSMO_THDMK}             >> $MF_FILE 
echo ${COSMO_MHDM}              >> $MF_FILE 
echo ${COSMO_FHDM}              >> $MF_FILE 


