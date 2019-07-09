#!/bin/bash

##########
#name:          compute_fragment_lengths.sh
#description:   compute fragment lengths of ATAC data
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 9, 2019
##########

TIME=$1
QUALITY_DIR=$2
SCRIPT_DIR=$3


echo "#####################"
echo $0 started on `hostname` at `date` with parameters $*
echo "Started processing:"
echo "#####################"


f1=$OUTPUT_DIR/ATAC/$TIME/mapping/*_R1_*filtered_final_1bp.bedpe
f2=$OUTPUT_DIR/ATAC/$TIME/mapping/*_R2_*filtered_final_1bp.bedpe

echo $f1
echo $f2

#plot done on bedpe file because we need paired end information 
#about fragment size
#count lines in bedpe  files for later normalization

NUM_LINES1=$(wc -l $f1 | cut -d ' ' -f 1)
echo $NUM_LINES1
NUM_LINES2=$(wc -l $f2 | cut -d ' ' -f 1)
echo $NUM_LINES2


#SHORT_NAME=${NAME:17}
#SHORTER_NAME=${SHORT_NAME::${#SHORT_NAME}-26}
#echo $SHORT_NAME
#echo $SHORTER_NAME

#get fragment lengths
#count how often they occur

awk 'BEGIN {OFS="\t"} {print ($6-$2)}' $f1 > \
    $QUALITY_DIR/ATAC_${TIME}_${NAME}_fragment_lengths_rep1.txt

sort $QUALITY_DIR/ATAC_${TIME}_${NAME}_fragment_lengths_rep1.txt | uniq -c | \
    sort -k2n,2n > \
    $QUALITY_DIR/ATAC_${TIME}_${NAME}_fragment_lengths_counts_rep1.txt


awk 'BEGIN {OFS="\t"} {print ($6-$2)}' $f2 > \
    $QUALITY_DIR/ATAC_${TIME}_${NAME}_fragment_lengths_rep2.txt

sort $QUALITY_DIR/ATAC_${TIME}_${NAME}_fragment_lengths_rep2.txt | uniq -c | \
    sort -k2n,2n > \
    $QUALITY_DIR/ATAC_${TIME}_${NAME}_fragment_lengths_counts_rep2.txt

Rscript ${SCRIPT_DIR}/ATAC/plot_fragment_lengths.r $QUALITY_DIR/ATAC_${TIME}_${NAME}_fragment_lengths_counts_rep1.txt $QUALITY_DIR/ATAC_${TIME}_${NAME}_fragment_lengths_counts_rep2.txt $NAME ${TIME} $NUM_LINES1 $NUM_LINES2 $QUALITY_DIR


exit

