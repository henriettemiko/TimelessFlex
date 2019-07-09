#!/bin/bash

##########
#name:          01c_compute_fragment_lengths_ATAC.sh
#description:   call computing of fragment lengths
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 9, 2019
##########

source ../set_variables_hg19.sh


QUALITY_DIR=$OUTPUT_DIR/quality_ATAC
mkdir -p $QUALITY_DIR

for f in $OUTPUT_DIR/ATAC/D*
do

    echo $f

    TIME=$(echo $f | rev | cut -d "/" -f 1 | rev )
    echo $TIME

     qsub -hold_jid "process_ATAC_PE_fastq_*" -N compute_fragment_lengths_${TIME} \
	     -V -j y -o $QUALITY_DIR/compute_fragment_lengths_ATAC_${TIME}.txt \
	     -cwd -pe smp 1 \
	     -l mem_free=1G,h_vmem=1G $SCRIPT_DIR/ATAC/compute_fragment_lengths.sh \
	     $TIME $QUALITY_DIR $SCRIPT_DIR

done


exit

