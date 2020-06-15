#!/bin/bash

##########
#name:          set_variables_mm10.sh
#description:   set variables for mm10 and download missing files
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          April 4, 2020
##########


echo $0 started on `hostname` at `date` with parameters $*


#set variables for directories

INPUT_DIR=/fast/AG_Ohler/henriette/hematopoiesis/input
START_DIR=/fast/AG_Ohler/henriette/hematopoiesis/

SCRIPT_DIR=$START_DIR/scripts

source $SCRIPT_DIR/set_guix_profile.sh


OUTPUT_DIR=$START_DIR/output_mm10
mkdir -p $OUTPUT_DIR

OPEN_REGIONS_DIR=$OUTPUT_DIR/open_regions
mkdir -p $OPEN_REGIONS_DIR


MACS2_ORGANISM=mm


##########
#ANNOTATION
ANNOTATION_DIR=$INPUT_DIR/annotation
mkdir -p $ANNOTATION_DIR

if [ ! -f $ANNOTATION_DIR/gencode.vM24.annotation.gtf ]
then
    #Gencode annotation https://www.gencodegenes.org
    wget -P $ANNOTATION_DIR \
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/"\
"gencode.vM24.annotation.gtf.gz"
    gunzip $ANNOTATION_DIR/gencode.vM24.annotation.gtf.gz
#else
    #echo "annotation is already there"
fi
ANNOTATION=$ANNOTATION_DIR/gencode.vM24.annotation.gtf 


#BLACKLIST
BLACKLIST_DIR=$INPUT_DIR/blacklist
mkdir -p $BLACKLIST_DIR
if [ ! -f $BLACKLIST_DIR/mm10-blacklist.v2.bed ]
then
    #blacklist v2 from 
    #https://sites.google.com/site/anshulkundaje/projects/blacklists
    wget -P $BLACKLIST_DIR \
"https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz"
    gunzip $BLACKLIST_DIR/mm10-blacklist.v2.bed.gz
#else
    #echo "blacklist is already there"
fi
BLACKLIST=$BLACKLIST_DIR/mm10-blacklist.v2.bed
##########


##########
#CHROMOSOME SIZES
CHR_SIZES_DIR=$INPUT_DIR/chr_sizes
mkdir -p $CHR_SIZES_DIR
if [ ! -f $CHR_SIZES_DIR/mm10_chr_sizes.txt ]
then
    wget -O $CHR_SIZES_DIR/mm10_chromInfo.txt.gz \
        "http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/"\
"chromInfo.txt.gz"
    gunzip $CHR_SIZES_DIR/mm10_chromInfo.txt.gz
    grep -v hap $CHR_SIZES_DIR/mm10_chromInfo.txt | grep -v Un | \
        grep -v random | cut -f1-2 > $CHR_SIZES_DIR/mm10_chr_sizes.txt
#else
    #echo "chr sizes file is already there"
fi
CHR_SIZES=$CHR_SIZES_DIR/mm10_chr_sizes.txt
##########


##########
#BOWTIE2 INDEX
BOWTIE2_INDEX_DIR=$INPUT_DIR/Bowtie2_index
mkdir -p $BOWTIE2_INDEX_DIR

if [ ! -f $BOWTIE2_INDEX_DIR/mm10.1.bt2 ]
then

    wget -P $BOWTIE2_INDEX_DIR \
        ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip
    unzip -d $BOWTIE2_INDEX_DIR $BOWTIE2_INDEX_DIR/mm10.zip

#else
    #echo "Bowtie2 index is already there"
fi

BOWTIE2_INDEX=$BOWTIE2_INDEX_DIR/mm10
##########


##########
##PRIMARY ASSEMBLY AND RSEM REFERENCE
#PRIMARY_ASSEMBLY_DIR=$INPUT_DIR/primary_assembly
#mkdir -p $PRIMARY_ASSEMBLY_DIR
#
#if [ ! -f $PRIMARY_ASSEMBLY_DIR/GRCm38.primary_assembly.genome.fa ]
#then
#    wget -P $PRIMARY_ASSEMBLY_DIR \
#        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/"\
#"GRCm38.primary_assembly.genome.fa.gz"
#    gunzip $PRIMARY_ASSEMBLY_DIR/GRCm38.primary_assembly.genome.fa.gz
##else
#    #echo "primary assembly is already there"
#fi
#
#RSEM_REF_DIR=$INPUT_DIR/RSEM_reference
#mkdir -p $RSEM_REF_DIR
#
#if [ ! -f $RSEM_REF_DIR/RSEM_ref.1.ebwt ]
#then
#    rsem-prepare-reference --gtf $ANNOTATION \
#        --bowtie $PRIMARY_ASSEMBLY_DIR/GRCm38.primary_assembly.genome.fa \
#        $RSEM_REF_DIR/RSEM_ref
##else
#    #echo "RSEM reference is already there"
#fi
#
#RSEM_REF=$RSEM_REF_DIR/RSEM_ref
###########

