#!/bin/bash

##########
#name:          02a_compute_feature_regions.sh
#description:   compute feature region windows around open regions subset
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 12, 2019
##########


source ../set_variables_hg19.sh

CLASSIFY_OPEN_REGIONS_DIR=$OPEN_REGIONS_DIR_SUB/classify_open_regions

FEATURE_REGIONS_DIR=$OPEN_REGIONS_DIR_SUB/get_feature_regions
mkdir -p $FEATURE_REGIONS_DIR


#feature regions should be windows 
#for all D0, D2, D5, D10 IDR peaks: extend  window 500 bp from downstream and
#upstream edge

cd $FEATURE_REGIONS_DIR

len=(500 1000)

for l in "${len[@]}"
do

    awk -v l=${l} 'BEGIN{OFS="\t"} {print $1,($2-l),($3+l),$4".ext",
    $5,$6,$7,$8,$9,$10,$11}' \
        $CLASSIFY_OPEN_REGIONS_DIR/peaks_enhancers_intergenic.bed > \
        peaks_enhancers_intergenic_ext_${l}.bed

    awk -v l=${l} 'BEGIN{OFS="\t"} {print $1,($2-l),($3+l),$4".ext",
    $5,$6,$7,$8,$9,$10,$11}' \
        $CLASSIFY_OPEN_REGIONS_DIR/peaks_enhancers_intragenic.bed > \
        peaks_enhancers_intragenic_ext_${l}.bed

    awk -v l=${l} 'BEGIN{OFS="\t"} {print $1,($2-l),($3+l),$4".ext",
    $5,$6,$7,$8,$9,$10,$11}' \
        $CLASSIFY_OPEN_REGIONS_DIR/peaks_unidirectional.bed > \
        peaks_unidirectional_ext_${l}.bed

    awk -v l=${l} 'BEGIN{OFS="\t"} {print $1,($2-l),($3+l),$4".ext",
    $5,$6,$7,$8,$9,$10,$11}' \
        $CLASSIFY_OPEN_REGIONS_DIR/peaks_bidirectional.bed > \
        peaks_bidirectional_ext_${l}.bed

    cat peaks_unidirectional_ext_${l}.bed \
        peaks_bidirectional_ext_${l}.bed \
        peaks_enhancers_intergenic_ext_${l}.bed \
        peaks_enhancers_intragenic_ext_${l}.bed | \
        sort -k1,1 -k2,2n > all_ext_${l}.bed

    #regions are book-ended, they should not be merged
    #only regions with at least one overlap should be merged
    bedtools merge -d -1 -c 4,5,6,7,8,9,10,11 -o distinct \
        -i all_ext_${l}.bed > all_ext_${l}_merged.bed

    awk 'OFS="\t" {if ($4!~",") print $0}' all_ext_${l}_merged.bed > \
        all_ext_${l}_final.bed

    #split in uni,bidir,intra,inter

    awk 'OFS="\t" {if ($8=="." && $9=="." && $10=="." && $11==".") print $0}' \
        all_ext_${l}_final.bed > all_ext_${l}_final_enhancers_intergenic.bed

    awk 'OFS="\t" {if ($8=="." && $9=="." && $10!="." && $11!=".") print $0}' \
        all_ext_${l}_final.bed > all_ext_${l}_final_enhancers_intragenic.bed

    awk 'OFS="\t" {if ($8!="." && $9!="." && $10!="." && $11!=".") print $0}' \
        all_ext_${l}_final.bed > all_ext_${l}_final_bidirectional.bed

    awk 'OFS="\t" {if ($8!="." && $9!="." && $10=="." && $11==".") print $0}' \
        all_ext_${l}_final.bed > all_ext_${l}_final_unidirectional.bed

    awk 'OFS="\t" {if ($8!="." && $9!=".") print $0}' all_ext_${l}_final.bed \
        > all_ext_${l}_final_promoters.bed

    awk 'OFS="\t" {if ($8=="." && $9==".") print $0}' all_ext_${l}_final.bed \
        > all_ext_${l}_final_enhancers.bed

    Rscript $SCRIPT_DIR/open_regions_subset/plot_feature_region_widths.r \
        all_ext_${l}_final.bed all_ext_${l}_final_promoters.bed \
        all_ext_${l}_final_enhancers.bed ${l} 

done


exit

