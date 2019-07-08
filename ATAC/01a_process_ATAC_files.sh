#!/bin/bash

##########
#name:          01a_process_ATAC_files.sh
#description:   call processing of ATAC fastq files
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


source ../set_variables_hg19.sh


for f in $INPUT_DIR/ATAC/D*/*
do

    if [[ "$f" =~ "_R1_" ]]
    then

        #make sure paired end read files are named the same
        #here first read named R1, second read R2     

        TIME=$(echo $f | rev | cut -d "/" -f 2 | rev)
        echo $TIME
        echo $f

        f1=$f
        f2=${f1/_R1_/_R2_}
        echo $f2

        NAME1=$(basename $f1 .fastq.gz)
        NAME2=$(basename $f2 .fastq.gz)
        echo $NAME1
        echo $NAME2

        ATAC_DIR=$OUTPUT_DIR/ATAC/$TIME/
        mkdir -p $ATAC_DIR

        qsub -V -j y -o $ATAC_DIR/process_ATAC_PE_fastq_${NAME1}.txt -cwd \
            -pe smp 1 -l h_rt=24:00:00 -l mem_free=30G,h_vmem=30G \
            $SCRIPT_DIR/ATAC/process_ATAC_PE_fastq.sh $f1 $f2 $NAME1 $NAME2 \
            $TIME $OUTPUT_DIR $BOWTIE2_INDEX $CHR_SIZES $SCRIPT_DIR $ATAC_DIR

    fi

done


exit

