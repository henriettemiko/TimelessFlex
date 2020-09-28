#!/bin/bash

##########
#name:          01a_process_RNA.sh
#description:   call processing for each RNA fastq file D8
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          August 8, 2019
##########


source ../set_variables_hg19.sh


for f in $INPUT_DIR/RNA/D8/*

do
    echo $f

    NAME=$(basename $f .fastq.gz)
    echo $NAME

    TIME=$(echo $f | rev| cut -d "/" -f 2 | rev)

    echo $TIME

    RNA_DIR=$OUTPUT_DIR/RNA/$TIME/
    mkdir -p $RNA_DIR

    qsub -N process_RNA_fastq_${NAME} \
        -V -j y -o $RNA_DIR/process_RNA_fastq_${NAME}.txt -cwd -l os=centos7 \
        -l mem_free=100G,h_vmem=100G $SCRIPT_DIR/RNA/process_RNA_fastq.sh $f \
        $NAME $RNA_DIR $RSEM_REF

done


exit

