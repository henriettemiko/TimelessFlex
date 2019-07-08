#!/bin/bash

##########
#name:          01a_process_ChIP_files.sh
#description:   trim illumina adapters, map with bowtie2, filter, 
#               convert to bed
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


f=$1 #filepath
NAME=$2
CHIP_DIR=$3
BOWTIE2_INDEX=$4

echo "#####################"
echo "Started processing:"
echo $0 started on `hostname` at `date` with parameters $*
echo $f
echo $NAME
echo $CHIP_DIR
echo $BOWTIE2_INDEX
echo "#####################"


##create directory structure

FASTQC_DIR=$CHIP_DIR/fastqc
TRIMMING_DIR=$CHIP_DIR/trimming
MAPPING_DIR=$CHIP_DIR/mapping
BEDGRAPHS_DIR=$CHIP_DIR/bedgraphs

mkdir -p $FASTQC_DIR
mkdir -p $TRIMMING_DIR
mkdir -p $MAPPING_DIR
mkdir -p $BEDGRAPHS_DIR


##run fastqc to get information about encoding
fastqc -o $FASTQC_DIR $f
unzip -q $FASTQC_DIR/${NAME}_fastqc.zip -d $FASTQC_DIR

#find out encoding of fastq file from fastqc output
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

##map with bowtie2
bowtie2 $QUALITY_SCORE -p 8 -x $BOWTIE2_INDEX \
    -U ${TRIMMING_DIR}/${NAME}_trimmed.fq.gz \
    -S ${MAPPING_DIR}/${NAME}_trimmed.sam

##filter for uniquely mapped with at most 2 mismatches
samtools view -@ 8 -Sh ${MAPPING_DIR}/${NAME}_trimmed.sam | \
    grep -e "^@" -e "XM:i:[012][^0-9]" | \
    grep -v "XS:i:" > ${MAPPING_DIR}/${NAME}_trimmed_filtered.sam

##convert to bam
samtools view -@ 8 -S -b ${MAPPING_DIR}/${NAME}_trimmed_filtered.sam \
    > ${MAPPING_DIR}/${NAME}_trimmed_filtered.bam

##sort
samtools sort -@ 8 -o ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted.bam \
    ${MAPPING_DIR}/${NAME}_trimmed_filtered.bam

##remove duplicates
samtools rmdup -s ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted.bam \
    ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup.bam 

##index
samtools index ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup.bam

##convert to bed
bedtools bamtobed -i ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup.bam \
    > ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup.bed

##make bedgraph and bigwig
bamCoverage -b ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup.bam \
    -o ${BEDGRAPHS_DIR}/${NAME}.bedGraph --outFileFormat bedgraph \
    --normalizeUsing RPKM

bamCoverage -b ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup.bam \
    -o ${BEDGRAPHS_DIR}/${NAME}.bigWig --outFileFormat bigwig \
    --normalizeUsing RPKM


exit

