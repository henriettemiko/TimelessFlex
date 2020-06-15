#!/bin/bash

##########
#name:          03a_compute_histone_mark_signals_split.sh
#description:   compute histone mark signals at feature regions split
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          August 9, 2019
##########


source ../set_variables_hg19.sh


FEATURE_REGIONS_DIR=$OUTPUT_DIR/open_regions_split/get_feature_regions
PROM_REGIONS=$FEATURE_REGIONS_DIR/all_ext_500_final_promoters.bed
ENH_REGIONS=$FEATURE_REGIONS_DIR/all_ext_500_final_enhancers.bed

echo "prom regions files: $PROM_REGIONS"
echo `wc -l $PROM_REGIONS`
echo "enh regions files: $ENH_REGIONS"
echo `wc -l $ENH_REGIONS`

#call Signal Generator for each mark and for each time point
#for promoter and enhancer feature regions

#for f in $OUTPUT_DIR/ChIP/H3K*/D*
#do
#    echo $f
#
#    MARK=$(echo $f | rev | cut -d "/" -f 2 | rev)
#    TIME=$(echo $f | rev | cut -d "/" -f 1 | rev)
#    echo $MARK
#    echo $TIME
#
#    SIGNAL_GENERATOR_DIR_PROM=$OPEN_REGIONS_DIR_SPLIT/Signal_Generator_promoters
#    mkdir -p $SIGNAL_GENERATOR_DIR_PROM
#
#    qsub -l os=centos7 -V -j y \
#        -o "$SIGNAL_GENERATOR_DIR_PROM/$MARK/$TIME/\
#call_Signal_Generator_prom.txt" -cwd -pe smp 1 \
#        -l mem_free=50G,h_vmem=50G \
#        $SCRIPT_DIR/open_regions_fullset/call_Signal_Generator.sh \
#        $MARK $TIME $SIGNAL_GENERATOR_DIR_PROM $OUTPUT_DIR $PROM_REGIONS \
#        $CHR_SIZES
#
#
#    SIGNAL_GENERATOR_DIR_ENH=$OPEN_REGIONS_DIR_SPLIT/Signal_Generator_enhancers
#    mkdir -p $SIGNAL_GENERATOR_DIR_ENH
#
#    qsub -l os=centos7 -V -j y \
#        -o "$SIGNAL_GENERATOR_DIR_ENH/$MARK/$TIME/\
#call_Signal_Generator_enh.txt" -cwd -pe smp 1 \
#        -l mem_free=50G,h_vmem=50G \
#        $SCRIPT_DIR/open_regions_fullset/call_Signal_Generator.sh \
#        $MARK $TIME $SIGNAL_GENERATOR_DIR_ENH $OUTPUT_DIR $ENH_REGIONS \
#        $CHR_SIZES
#
#done


#mv $OPEN_REGIONS_DIR_SPLIT/ChIP/H3K27ac/D8 $OUTPUT_DIR/ChIP/H3K27ac
#mv $OPEN_REGIONS_DIR_SPLIT/ChIP/H3K27me3/D8 $OUTPUT_DIR/ChIP/H3K27me3
#mv $OPEN_REGIONS_DIR_SPLIT/ChIP/H3K4me1/D8 $OUTPUT_DIR/ChIP/H3K4me1
#mv $OPEN_REGIONS_DIR_SPLIT/ChIP/H3K4me3/D8 $OUTPUT_DIR/ChIP/H3K4me3


#D8
for f in $OUTPUT_DIR/ChIP/H3K*/D8
do
    echo $f

    MARK=$(echo $f | rev | cut -d "/" -f 2 | rev)
    TIME=$(echo $f | rev | cut -d "/" -f 1 | rev)
    echo $MARK
    echo $TIME

    SIGNAL_GENERATOR_DIR_PROM=$OPEN_REGIONS_DIR_SPLIT/Signal_Generator_promoters
    mkdir -p $SIGNAL_GENERATOR_DIR_PROM

    qsub -l os=centos7 -V -j y \
        -o "$SIGNAL_GENERATOR_DIR_PROM/$MARK/$TIME/\
call_Signal_Generator_prom.txt" -cwd -pe smp 1 \
        -l mem_free=50G,h_vmem=50G \
        $SCRIPT_DIR/open_regions_fullset/call_Signal_Generator.sh \
        $MARK $TIME $SIGNAL_GENERATOR_DIR_PROM $OUTPUT_DIR $PROM_REGIONS \
        $CHR_SIZES


    SIGNAL_GENERATOR_DIR_ENH=$OPEN_REGIONS_DIR_SPLIT/Signal_Generator_enhancers
    mkdir -p $SIGNAL_GENERATOR_DIR_ENH

    qsub -l os=centos7 -V -j y \
        -o "$SIGNAL_GENERATOR_DIR_ENH/$MARK/$TIME/\
call_Signal_Generator_enh.txt" -cwd -pe smp 1 \
        -l mem_free=50G,h_vmem=50G \
        $SCRIPT_DIR/open_regions_fullset/call_Signal_Generator.sh \
        $MARK $TIME $SIGNAL_GENERATOR_DIR_ENH $OUTPUT_DIR $ENH_REGIONS \
        $CHR_SIZES

done

#mv $OUTPUT_DIR/ChIP/H3K27ac/D8 $OPEN_REGIONS_DIR_SPLIT/ChIP/H3K27ac
#mv $OUTPUT_DIR/ChIP/H3K27me3/D8 $OPEN_REGIONS_DIR_SPLIT/ChIP/H3K27me3
#mv $OUTPUT_DIR/ChIP/H3K4me1/D8 $OPEN_REGIONS_DIR_SPLIT/ChIP/H3K4me1
#mv $OUTPUT_DIR/ChIP/H3K4me3/D8 $OPEN_REGIONS_DIR_SPLIT/ChIP/H3K4me3


exit

