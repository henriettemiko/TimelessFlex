#!/bin/bash

##########
#name:          02b_plot_feature_regions_overlap.sh
#description:   plot overlap of feature regions with histone marks
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 10, 2019
##########


source ../set_variables_hg19.sh

FEATURE_REGIONS_DIR=$OPEN_REGIONS_DIR_FULL/get_feature_regions

cd $FEATURE_REGIONS_DIR

#compute the overlaps

hist=("H3K27ac" "H3K27me3" "H3K4me1" "H3K4me3")
time=("D0" "D2" "D5" "D7" "D10")

for h in "${hist[@]}"
do 

    for t in "${time[@]}"
    do

        sort -k1,1 -k2,2n \
            $OUTPUT_DIR/ChIP/$h/$t/JAMM_peaks/peaks/filtered.peaks.narrowPeak \
            > ${h}_${t}_sorted.bed

    done
done


len=(500 1000)

for l in "${len[@]}"
do

    sort -k1,1 -k2,2n all_ext_${l}_final.bed > all_ext_${l}_final_sorted.bed
    sort -k1,1 -k2,2n all_ext_${l}_final_promoters.bed > \
        all_ext_${l}_final_promoters_sorted.bed
    sort -k1,1 -k2,2n all_ext_${l}_final_enhancers.bed > \
        all_ext_${l}_final_enhancers_sorted.bed

    for h in "${hist[@]}"
    do

        #all regions with at least 1bp overlap with any hist peak
        paste all_ext_${l}_final_sorted.bed \
            <( bedops-typical --everything ${h}_*_sorted.bed | \
            bedmap-typical --indicator --echo-overlap-size --delim '\t' \
            all_ext_${l}_final_sorted.bed - ) > \
            dist_${l}_all_indicator_${h}_1bp.bed
        
        paste all_ext_${l}_final_promoters_sorted.bed \
            <( bedops-typical --everything ${h}_*_sorted.bed | \
            bedmap-typical --indicator --echo-overlap-size --delim '\t'\
            all_ext_${l}_final_promoters_sorted.bed - ) > \
            dist_${l}_all_promoters_indicator_${h}_1bp.bed

        paste all_ext_${l}_final_enhancers_sorted.bed \
            <( bedops-typical --everything ${h}_*_sorted.bed | \
            bedmap-typical --indicator --echo-overlap-size --delim '\t' \
            all_ext_${l}_final_enhancers_sorted.bed - ) > \
            dist_${l}_all_enhancers_indicator_${h}_1bp.bed

        awk 'BEGIN{OFS="\t"} {if($12=="1"){print $13}}' \
            dist_${l}_all_indicator_${h}_1bp.bed | sed 's/;/\t/g' > \
            feature_regions_${l}_${h}_overlaps.bed
        awk 'BEGIN{OFS="\t"} {if($12=="1"){print $13}}' \
            dist_${l}_all_promoters_indicator_${h}_1bp.bed | sed 's/;/\t/g' > \
            feature_regions_${l}_promoters_${h}_overlaps.bed
        awk 'BEGIN{OFS="\t"} {if($12=="1"){print $13}}' \
            dist_${l}_all_enhancers_indicator_${h}_1bp.bed | sed 's/;/\t/g' > \
            feature_regions_${l}_enhancers_${h}_overlaps.bed

    done

    Rscript "$SCRIPT_DIR/open_regions_fullset/\
plot_feature_region_widths_overlaps.r" $FEATURE_REGIONS_DIR ${l} "" "all"
    Rscript "$SCRIPT_DIR/open_regions_fullset/\
plot_feature_region_widths_overlaps.r" $FEATURE_REGIONS_DIR ${l} \
        "promoters_" "promoter"
    Rscript "$SCRIPT_DIR/open_regions_fullset/\
plot_feature_region_widths_overlaps.r" $FEATURE_REGIONS_DIR ${l} \
        "enhancers_" "enhancer"

done


exit

