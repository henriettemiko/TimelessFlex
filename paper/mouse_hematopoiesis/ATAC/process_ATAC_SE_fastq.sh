#!/bin/bash

##########
#name:          process_ATAC_SE_fastq.sh
#description:   process SE ATAC fastq files
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          March 27, 2020
##########


f=$1
NAME=$2
TIME=$3
OUTPUT_DIR=$4
BOWTIE2_INDEX=$5
CHR_SIZES=$6
SCRIPT_DIR=$7
ATAC_DIR=$8

echo "#####################"
echo $0 started on `hostname` at `date` with parameters $*
echo "Started processing:"
echo "current time point: $TIME"
echo "current file name: $NAME"
echo "#####################"


##create directory structure

FASTQC_DIR=$ATAC_DIR/fastqc
TRIMMING_DIR=$ATAC_DIR/trimming
MAPPING_DIR=$ATAC_DIR/mapping
BEDGRAPHS_DIR=$ATAC_DIR/bedgraphs


mkdir -p $FASTQC_DIR
mkdir -p $TRIMMING_DIR
mkdir -p $MAPPING_DIR
mkdir -p $BEDGRAPHS_DIR


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


##trim nextera adapters
trim_galore --nextera --stringency 3 $QUALITY_SCORE --fastqc \
    -o $TRIMMING_DIR --gzip $f

##map with bowtie2
bowtie2 $QUALITY_SCORE -p 8 \
    -x $BOWTIE2_INDEX \
    -U ${TRIMMING_DIR}/${NAME}_trimmed.fq.gz \
    -S ${MAPPING_DIR}/${NAME}_trimmed.sam

##filter for uniquely mapped with at most 2 mismatches
samtools view -@ 8 -Sh ${MAPPING_DIR}/${NAME}_trimmed.sam | \
    grep -e "^@" -e "XM:i:[012][^0-9]" | \
    grep -v "XS:i:" > ${MAPPING_DIR}/${NAME}_trimmed_filtered.sam

##convert to bam
samtools view -@ 8 -S -b ${MAPPING_DIR}/${NAME}_trimmed_filtered.sam > \
    ${MAPPING_DIR}/${NAME}_trimmed_filtered.bam

##sort
samtools sort -@ 8 -o ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted.bam \
    ${MAPPING_DIR}/${NAME}_trimmed_filtered.bam

##remove duplicates
samtools rmdup -s ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted.bam \
    ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup.bam

##index
samtools index ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup.bam

##make bedgraph and bigwig
bamCoverage -b ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup.bam \
    -o ${BEDGRAPHS_DIR}/${NAME}.bedGraph --outFileFormat bedgraph \
    --normalizeUsing RPKM

bamCoverage -b ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup.bam \
    -o ${BEDGRAPHS_DIR}/${NAME}.bigWig --outFileFormat bigwig \
    --normalizeUsing RPKM

##convert to bed
bedtools bamtobed \
    -i ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup.bam > \
    ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup.bed

##filtering
grep -v Un_gl ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup.bed | \
    grep -v random | grep -v chrM | grep -v hap > \
    ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup_final.bed

##cutting to 5' ends
awk 'OFS="\t" {if ($6=="+"){print $1,$2,$2+1,$4,$5,$6;} else print $1,$3-1,$3,$4,$5,$6;}'  \
    ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup_final.bed > \
    ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup_final_1bp.bed

##convert to bam
bedtools bedtobam \
    -i ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup_final_1bp.bed \
    -g $CHR_SIZES > \
    ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup_final_1bp_frombed.bam

##sorting
samtools sort \
    -o ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup_final_1bp_frombed_sorted.bam \
    ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup_final_1bp_frombed.bam

##indexing
samtools index \
    ${MAPPING_DIR}/${NAME}_trimmed_filtered_sorted_nodup_final_1bp_frombed_sorted.bam


exit

