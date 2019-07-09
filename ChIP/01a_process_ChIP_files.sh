#!/bin/bash

##########
#name:          01a_process_ChIP_files.sh
#description:   call processing for each ChIP fastq file
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


source ../set_variables_hg19.sh


for f in $INPUT_DIR/ChIP/*/D*/*
do
    echo $f

    NAME=$(basename $f .fastq.gz)
    echo $NAME

    MARK=$(echo $f | rev | cut -d "/" -f 3 | rev)
    TIME=$(echo $f | rev| cut -d "/" -f 2 | rev)

    echo $MARK
    echo $TIME

    CHIP_DIR=$OUTPUT_DIR/ChIP/$MARK/$TIME/
    mkdir -p $CHIP_DIR

    qsub -N process_ChIP_fastq_${NAME} -V -j y \
        -o $CHIP_DIR/process_ChIP_fastq_${NAME}.txt -cwd -pe smp 1 \
        -l mem_free=50G,h_vmem=50G $SCRIPT_DIR/ChIP/process_ChIP_fastq_SE.sh \
        $f $NAME $CHIP_DIR $BOWTIE2_INDEX

done


exit

