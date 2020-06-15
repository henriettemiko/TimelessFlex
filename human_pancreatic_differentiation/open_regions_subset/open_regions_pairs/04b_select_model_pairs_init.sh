#!/bin/bash

##########
#name:          04b_select_model_pairs_init.sh
#description:   calls plotting of model selection criteria for pairs init
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 17, 2019
##########


source ../../set_variables_hg19.sh

SIGNAL_GENERATOR_DIR_INIT_PROM_ENH=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_init_prom-enh

SIGNAL_GENERATOR_DIR_INIT_PROM_PROM=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_init_prom-prom

SIGNAL_GENERATOR_DIR_INIT_ENH_ENH=$OPEN_REGIONS_DIR_PAIRS/\
Signal_Generator_init_enh-enh


TIMELESS_DIR_INIT_PROM_ENH=$OPEN_REGIONS_DIR_PAIRS/timeless_init_prom-enh
TIMELESS_DIR_INIT_PROM_PROM=$OPEN_REGIONS_DIR_PAIRS/timeless_init_prom-prom
TIMELESS_DIR_INIT_ENH_ENH=$OPEN_REGIONS_DIR_PAIRS/timeless_init_enh-enh


cd $TIMELESS_DIR_INIT_PROM_ENH

cat 2_30/AIC_*.mat |sort -k1,1n | cut -f2 > AIC.txt
cat 2_30/BIC_*.mat |sort -k1,1n | cut -f2 > BIC.txt
cat 2_30/likelihood_*.mat |sort -k1,1n | cut -f2 > likelihood.txt
cat 2_30/numPar_*.mat |sort -k1,1n | cut -f2 > numPar.txt

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/plot_AIC_BIC.r \
    AIC.txt BIC.txt $TIMELESS_DIR_INIT_PROM_ENH likelihood.txt numPar.txt


cd $TIMELESS_DIR_INIT_PROM_PROM

cat 2_30/AIC_*.mat |sort -k1,1n | cut -f2 > AIC.txt
cat 2_30/BIC_*.mat |sort -k1,1n | cut -f2 > BIC.txt
cat 2_30/likelihood_*.mat |sort -k1,1n | cut -f2 > likelihood.txt
cat 2_30/numPar_*.mat |sort -k1,1n | cut -f2 > numPar.txt

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/plot_AIC_BIC.r \
    AIC.txt BIC.txt $TIMELESS_DIR_INIT_PROM_PROM likelihood.txt numPar.txt


cd $TIMELESS_DIR_INIT_ENH_ENH

cat 2_30/AIC_*.mat |sort -k1,1n | cut -f2 > AIC.txt
cat 2_30/BIC_*.mat |sort -k1,1n | cut -f2 > BIC.txt
cat 2_30/likelihood_*.mat |sort -k1,1n | cut -f2 > likelihood.txt
cat 2_30/numPar_*.mat |sort -k1,1n | cut -f2 > numPar.txt

Rscript $SCRIPT_DIR/open_regions_subset/open_regions_pairs/plot_AIC_BIC.r \
    AIC.txt BIC.txt $TIMELESS_DIR_INIT_ENH_ENH likelihood.txt numPar.txt


exit

