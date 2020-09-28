#!/bin/bash

##########
#name:          process_RNA_fastq.sh
#description:   trim Illumina adapters and call RSEM
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


echo $0 started on `hostname` at `date` with parameters $*

f=$1
NAME=$2
RNA_DIR=$3
RSEM_REF=$4

echo "#####################"
echo "Started processing:"
echo "current file name: $NAME"
echo "#####################"

FASTQC_DIR=$RNA_DIR/fastqc
TRIMMING_DIR=$RNA_DIR/trimming
RSEM_DIR=$RNA_DIR/RSEM


mkdir -p $FASTQC_DIR
mkdir -p $TRIMMING_DIR
mkdir -p $RSEM_DIR


##run fastqc to get information about encoding
fastqc -o $FASTQC_DIR $f
unzip -q $FASTQC_DIR/${NAME}_fastqc.zip -d $FASTQC_DIR

##find out encoding of fastq file from fastqc output
ENCODING=$(grep Encoding $FASTQC_DIR/${NAME}_fastqc/fastqc_data.txt)
echo $ENCODING
if [[ $ENCODING =~ "Illumina 1.5" ]]
then
    QUALITY_SCORE=--phred64
    echo $QUALITY_SCORE
elif [[ $ENCODING =~ "Illumina 1.9" ]]
then
    QUALITY_SCORE=--phred33
    echo $QUALITY_SCORE
else
    echo "unkown encoding: $ENCODING"
    exit
fi

##trim Illumina adapters
trim_galore --illumina --stringency 3 $QUALITY_SCORE --fastqc \
    -o $TRIMMING_DIR --gzip $f


mkdir -p $RSEM_DIR/${NAME}_tmp

#forward prob needs to be determined by RSeQc
rsem-calculate-expression -p 8 $QUALITY_SCORE-quals --forward-prob=0 \
    --temporary-folder $RSEM_DIR/${NAME}_tmp --calc-ci \
    <( zcat ${TRIMMING_DIR}/${NAME}_trimmed.fq.gz ) $RSEM_REF \
    $RSEM_DIR/${NAME}_trimmed


exit

