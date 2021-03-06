#!/bin/bash

##########
#name:          04b_select_model_split.sh
#description:   calls plotting of model selection criteria split
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          August 12, 2019
##########


source ../set_variables_hg19.sh


SIGNAL_GENERATOR_DIR_PROM=$OPEN_REGIONS_DIR_SPLIT/Signal_Generator_promoters
SIGNAL_GENERATOR_DIR_ENH=$OPEN_REGIONS_DIR_SPLIT/Signal_Generator_enhancers
TIMELESS_DIR_PROM=$OPEN_REGIONS_DIR_SPLIT/timeless_promoters
TIMELESS_DIR_ENH=$OPEN_REGIONS_DIR_SPLIT/timeless_enhancers


cd $TIMELESS_DIR_PROM

cat 2_30/AIC_*.mat |sort -k1,1n | cut -f2 > AIC.txt
cat 2_30/BIC_*.mat |sort -k1,1n | cut -f2 > BIC.txt
cat 2_30/likelihood_*.mat |sort -k1,1n | cut -f2 > likelihood.txt
cat 2_30/numPar_*.mat |sort -k1,1n | cut -f2 > numPar.txt

Rscript $SCRIPT_DIR/open_regions_fullset/plot_AIC_BIC.r AIC.txt BIC.txt \
    $TIMELESS_DIR_PROM likelihood.txt numPar.txt

exit
cd $TIMELESS_DIR_ENH

cat 2_30/AIC_*.mat |sort -k1,1n | cut -f2 > AIC.txt
cat 2_30/BIC_*.mat |sort -k1,1n | cut -f2 > BIC.txt
cat 2_30/likelihood_*.mat |sort -k1,1n | cut -f2 > likelihood.txt
cat 2_30/numPar_*.mat |sort -k1,1n | cut -f2 > numPar.txt

Rscript $SCRIPT_DIR/open_regions_fullset/plot_AIC_BIC.r AIC.txt BIC.txt \
    $TIMELESS_DIR_ENH likelihood.txt numPar.txt


exit

