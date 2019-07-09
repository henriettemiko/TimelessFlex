#!/bin/bash

##########
#name:          compute_fragment_lengths.sh
#description:   compute fragment lengths of ATAC data
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 9, 2019
##########

f=$1
NAME=$2
TIME=$3
QUALITY_DIR=$4
SCRIPT_DIR=$5


echo "#####################"
echo $0 started on `hostname` at `date` with parameters $*
echo "Started processing:"
echo "current file name: $NAME"
echo "#####################"



#plot done on bedpe file because we need paired end information 
#about fragment size
#count lines in bedpe  files for later normalization

NUM_LINES=$(wc -l $f | cut -d ' ' -f 1)
echo $NUM_LINES


#SHORT_NAME=${NAME:17}
#SHORTER_NAME=${SHORT_NAME::${#SHORT_NAME}-26}
#echo $SHORT_NAME
#echo $SHORTER_NAME

#get fragment lengths
#count how often they occur

awk 'BEGIN {OFS="\t"} {print ($6-$2)}' $f > \
    $QUALITY_DIR/ATAC_${TIME}_${NAME}_fragment_lengths.txt

sort $QUALITY_DIR/ATAC_${TIME}_${NAME}_fragment_lengths.txt | uniq -c | \
    sort -k2n,2n > \
    $QUALITY_DIR/ATAC_${TIME}_${NAME}_fragment_lengths_counts.txt

Rscript ${SCRIPT_DIR}/ATAC/plot_fragment_lengths.r $NAME ${TIME} $NUM_LINES \
    $QUALITY_DIR


exit

