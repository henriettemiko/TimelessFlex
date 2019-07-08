#!/bin/bash

##########
#name:          process_ATAC_PE_fastq.sh
#description:   process ATAC fastq files
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


f1=$1
f2=$2
NAME1=$3
NAME2=$4
TIME=$5
OUTPUT_DIR=$6
BOWTIE2_INDEX=$7
CHR_SIZES=$8
SCRIPT_DIR=$9
ATAC_DIR=${10}

echo "#####################"
echo $0 started on `hostname` at `date` with parameters $*
echo "Started processing:"
echo "current time point: $TIME"
echo "current file name: $NAME1"
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
##only on one replicate
fastqc -o $FASTQC_DIR $f1
unzip -q $FASTQC_DIR/${NAME1}_fastqc.zip -d $FASTQC_DIR

##find out encoding of fastq file from fastqc output
ENCODING=$(grep Encoding $FASTQC_DIR/${NAME1}_fastqc/fastqc_data.txt)
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
trim_galore --paired --nextera --stringency 3 $QUALITY_SCORE --fastqc \
    -o $TRIMMING_DIR --gzip $f1 $f2

##map with bowtie2
bowtie2 $QUALITY_SCORE -p 8 --no-discordant --no-mixed -X 2000 \
    -x $BOWTIE2_INDEX -1 ${TRIMMING_DIR}/${NAME1}_val_1.fq.gz \
    -2 ${TRIMMING_DIR}/${NAME2}_val_2.fq.gz \
    -S ${MAPPING_DIR}/${NAME1}_trimmed.sam

##filtering by mapping quality and converting to bam
samtools view -S -b -q 20 ${MAPPING_DIR}/${NAME1}_trimmed.sam > \
    ${MAPPING_DIR}/${NAME1}_trimmed_q20.bam

samtools flagstat ${MAPPING_DIR}/${NAME1}_trimmed_q20.bam

##remove duplicates with Picard
java -Xmx2G -jar ${GUIX_PROFILE}/share/java/picard.jar SortSam \
    INPUT=${MAPPING_DIR}/${NAME1}_trimmed_q20.bam \
    OUTPUT=${MAPPING_DIR}/${NAME1}_trimmed_q20_picardsorted.bam \
    SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT > \
    ${MAPPING_DIR}/${NAME1}_trimmed_q20_picardsort.txt 2>&1

java -Xmx8G -jar ${GUIX_PROFILE}/share/java/picard.jar MarkDuplicates \
    INPUT=${MAPPING_DIR}/${NAME1}_trimmed_q20_picardsorted.bam \
    OUTPUT=${MAPPING_DIR}/${NAME1}_trimmed_q20_picardsorted_nodup.bam \
    METRICS_FILE=${MAPPING_DIR}/${NAME1}_trimmed_q20_picardsorted_nodup.txt \
    REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT > \
    ${MAPPING_DIR}/${NAME1}_trimmed_q20_picardmark.txt 2>&1

##indexing
samtools index ${MAPPING_DIR}/${NAME1}_trimmed_q20_picardsorted_nodup.bam

##namesorting
samtools sort -n \
    -o ${MAPPING_DIR}/${NAME1}_trimmed_q20_picardsorted_nodup_namesorted.bam \
    ${MAPPING_DIR}/${NAME1}_trimmed_q20_picardsorted_nodup.bam

##converting to bedpe
bedtools bamtobed -bedpe \
    -i ${MAPPING_DIR}/${NAME1}_trimmed_q20_picardsorted_nodup_namesorted.bam \
    > ${MAPPING_DIR}/${NAME1}.bedpe

#filtering
#chrUn for fly chromosomes
grep -v Un_gl ${MAPPING_DIR}/${NAME1}.bedpe | grep -v random | grep -v chrM | \
    grep -v hap | grep -v chrUn > ${MAPPING_DIR}/${NAME1}_filtered.bedpe

##at least 38 bp
awk 'OFS="\t" {if(($6-$2)>=38) print $0,$6-$2}' \
    ${MAPPING_DIR}/${NAME1}_filtered.bedpe > \
    ${MAPPING_DIR}/${NAME1}_filtered_final.bedpe

##convert to 1 bp
awk 'OFS="\t" {print $1,$2,$2+1,$4,$6-1,$6,$7,$8,$9,$10,$11}' \
    ${MAPPING_DIR}/${NAME1}_filtered_final.bedpe > \
    ${MAPPING_DIR}/${NAME1}_filtered_final_1bp.bedpe

##convert to bed
awk 'OFS="\t" {print $1,$2,$3,$7,$8,$9"\n"$4,$5,$6,$7,$8,$10}' \
    ${MAPPING_DIR}/${NAME1}_filtered_final_1bp.bedpe | sort -k1,1 -k2,2n > \
    ${MAPPING_DIR}/${NAME1}_filtered_final_1bp.bed

##converting to bam
bedtools bedtobam -i ${MAPPING_DIR}/${NAME1}_filtered_final_1bp.bed \
    -g $CHR_SIZES > ${MAPPING_DIR}/${NAME1}_filtered_final_1bp_frombed.bam

##sorting
samtools sort \
    -o ${MAPPING_DIR}/${NAME1}_filtered_final_1bp_frombed_sorted.bam \
    ${MAPPING_DIR}/${NAME1}_filtered_final_1bp_frombed.bam

##indexing
samtools index ${MAPPING_DIR}/${NAME1}_filtered_final_1bp_frombed_sorted.bam

##make bedgraph
bamCoverage -b ${MAPPING_DIR}/${NAME1}_trimmed_q20_picardsorted_nodup.bam \
    -o ${BEDGRAPHS_DIR}/${NAME1}.bedGraph --outFileFormat bedgraph \
    --normalizeUsing RPKM

##make bigwig
bamCoverage -b ${MAPPING_DIR}/${NAME1}_trimmed_q20_picardsorted_nodup.bam \
    -o ${BEDGRAPHS_DIR}/${NAME1}.bigWig --outFileFormat bigwig \
    --normalizeUsing RPKM


exit

