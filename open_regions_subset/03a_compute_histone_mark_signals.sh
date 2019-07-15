#!/bin/bash

##########
#name:          03a_compute_histone_mark_signals.sh
#description:   compute histone mark signals at feature regions subset
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 12, 2019
##########


source ../set_variables_hg19.sh


FEATURE_REGIONS_DIR=$OUTPUT_DIR/open_regions_subset/get_feature_regions
PROM_REGIONS=$FEATURE_REGIONS_DIR/all_ext_500_final_promoters.bed
ENH_REGIONS=$FEATURE_REGIONS_DIR/all_ext_500_final_enhancers.bed

echo "prom regions files: $PROM_REGIONS"
echo `wc -l $PROM_REGIONS`
echo "enh regions files: $ENH_REGIONS"
echo `wc -l $ENH_REGIONS`

#call Signal Generator for each mark and for each time point
#for promoter and enhancer feature regions

for f in $OUTPUT_DIR/ChIP/H3K*/D*
do
    echo $f

    if [[ ! $f =~ "D7" ]]
        then

        MARK=$(echo $f | rev | cut -d "/" -f 2 | rev)
        TIME=$(echo $f | rev | cut -d "/" -f 1 | rev)
        echo $MARK
        echo $TIME
    
#        SIGNAL_GENERATOR_DIR_PROM=$OPEN_REGIONS_DIR_SUB/\
#Signal_Generator_promoters
#        mkdir -p $SIGNAL_GENERATOR_DIR_PROM
#
#        qsub -l os=centos7 -V -j y \
#            -o "$SIGNAL_GENERATOR_DIR_PROM/$MARK/$TIME/\
#call_Signal_Generator_prom.txt" -cwd -pe smp 1 \
#            -l mem_free=50G,h_vmem=50G \
#            $SCRIPT_DIR/open_regions_subset/call_Signal_Generator.sh \
#            $MARK $TIME $SIGNAL_GENERATOR_DIR_PROM $OUTPUT_DIR $PROM_REGIONS \
#            $CHR_SIZES
#

        SIGNAL_GENERATOR_DIR_ENH=$OPEN_REGIONS_DIR_SUB/\
Signal_Generator_enhancers
        mkdir -p $SIGNAL_GENERATOR_DIR_ENH

        qsub -l os=centos7 -V -j y \
            -o "$SIGNAL_GENERATOR_DIR_ENH/$MARK/$TIME/\
call_Signal_Generator_enh.txt" -cwd -pe smp 1 \
            -l mem_free=50G,h_vmem=50G \
            $SCRIPT_DIR/open_regions_subset/call_Signal_Generator.sh \
            $MARK $TIME $SIGNAL_GENERATOR_DIR_ENH $OUTPUT_DIR $ENH_REGIONS \
            $CHR_SIZES


    fi
done


exit

