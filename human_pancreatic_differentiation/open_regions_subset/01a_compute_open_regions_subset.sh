#!/bin/bash

##########
#name:          01a_compute_open_regions_subset.sh
#description:   compute open regions for subset (D0, D2, D5, D10)
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 12, 2019
##########


source ../set_variables_hg19.sh


#combine ATAC peaks to open regions
#loop through time points

COMBINE_PEAKS_DIR=$OPEN_REGIONS_DIR_SUB/combine_peaks
mkdir -p $COMBINE_PEAKS_DIR

ANNOTATION_DIR=$OPEN_REGIONS_DIR_SUB/annotation
mkdir -p $ANNOTATION_DIR

CLASSIFY_OPEN_REGIONS_DIR=$OPEN_REGIONS_DIR_SUB/classify_open_regions
mkdir -p $CLASSIFY_OPEN_REGIONS_DIR

#HERE ONLY TAKE D0, D2, D5 and D10
cp $OUTPUT_DIR/ATAC/D0/MACS2_peaks/\
D0_pooled_p005_peaks_IDR_summitpos_final.narrowPeak $COMBINE_PEAKS_DIR

cp $OUTPUT_DIR/ATAC/D2/MACS2_peaks/\
D2_pooled_p005_peaks_IDR_summitpos_final.narrowPeak $COMBINE_PEAKS_DIR

cp $OUTPUT_DIR/ATAC/D5/MACS2_peaks/\
D5_pooled_p005_peaks_IDR_summitpos_final.narrowPeak $COMBINE_PEAKS_DIR

cp $OUTPUT_DIR/ATAC/D10/MACS2_peaks/\
D10_pooled_p005_peaks_IDR_summitpos_final.narrowPeak $COMBINE_PEAKS_DIR


#for each time point take peak calls from MACS2 and merge over time points
#keep summit positions and compute median of summit positions (not offsets!)
#overlap of 101 so we don't have transitive overlaps over time points
#(minimum peak length is 200)
#if peaks overlap more than half of minimal size they get merged, 
#other overlaps stay
#col11 is summit position in peak files
#after merge column 5 is median summit position

cd $COMBINE_PEAKS_DIR

cat *_pooled_p005_peaks_IDR_summitpos_final.narrowPeak | sort -k1,1 -k2,2n > \
    all.narrowPeak

bedtools merge -c 4,11 -o distinct,median -i all.narrowPeak -d -101 > \
    all_merged.bed

#median could not be integer, round with awk int
awk 'OFS="\t" {print $1,$2,$3,$4,int($5)}' all_merged.bed > \
    peaks_beforeblacklist.bed

#check if they overlap with blacklist regions
bedtools intersect -v -wa -a peaks_beforeblacklist.bed -b $BLACKLIST > \
    peaks.bed


#############

#now Hi-C comes into play
#compute fully overlaps of Hi-C bins and open chromatin regions/peaks
#full peak must be overlapped by bin (if not then difficult to decide which 
#Hi-C information is used for this peak
bedtools intersect -wo -F 1 -a $HIC_DIR/all_bins_unique.bed -b peaks.bed > \
all_bins_unique_overlapping_peaks.bed
#col1-6 Hi-C bin, col 7- peak information, last col overlap
#all overlaps of unique bins with peaks
#for each overlap with a peak the bin occurs once in file

#number of occurences of bin is number of overlaps
#count occurences and plot
cut -f 1-3 all_bins_unique_overlapping_peaks.bed | sort | uniq -c > \
all_bins_unique_overlapping_peaks_counts.txt
#NUMBER OF UNIQUE BINS THAT OVERLAP AT LEAST ONE PEAK FULLY (OPEN CHROMATIN REGION)

awk '{print $1}' all_bins_unique_overlapping_peaks_counts.txt | sort | uniq -c > overlap_counts.txt

#plot how many peaks Hi-C bins are overlapping fully
Rscript $SCRIPT_DIR/open_regions_subset/plot_bins_overlaps.r \
all_bins_unique_overlapping_peaks_counts.txt


#get unique peaks
cut -f 7-11 all_bins_unique_overlapping_peaks.bed | sort | uniq | \
sort -k1,1 -k2,2n > peaks_overlapped.bed

#get unique bins
cut -f 1-6 all_bins_unique_overlapping_peaks.bed | sort | uniq | \
sort -k1,1 -k2,2n > bins_overlapping.bed
#peaks_overlapped.bed is then divided into promoters and enhancers


##########

#for later: notunique bins
bedtools intersect -wo -F 1 -a $HIC_DIR/all_bins_notunique.bed -b peaks.bed > all_bins_notunique_overlapping_peaks.bed

#out of curiosity: how many bins are not overlapping
bedtools intersect -v -a $HIC_DIR/all_bins_unique.bed -b peaks.bed > \
all_bins_unique_nonoverlapping_peaks.bed

#########


exit


cd $ANNOTATION_DIR

#get annotation files
#get ALL transcripts from gencode
awk '{if($3=="transcript"){print $0}}' $ANNOTATION > all_transcripts.gtf

#convert to bed
awk 'OFS="\t" {print $1,$4-1,$5,$10,$14,$7}' all_transcripts.gtf | \
    tr -d '";' > all_transcripts.bed

#get TSSs
awk 'BEGIN{OFS="\t"} {if($6=="+"){print $1,$2,($2+1),$4,$5,$6} 
else {print $1,($3-1),$3,$4,$5,$6}}' all_transcripts.bed > \
    all_transcript_TSSs.bed

grep ^chr all_transcript_TSSs.bed | sort -k1,1 -k2,2n | uniq > \
    all_transcript_TSSs_chr_uniq.bed

grep ^chr all_transcripts.bed | sort -k1,1 -k2,2n > all_transcripts_chr.bed


cd $CLASSIFY_OPEN_REGIONS_DIR

#find closest TSS to each peak on plus and minus strand
#classify peaks based on overlapping TSSs

#look at peaks once as if on plus strand and once as if on minus strand
awk 'OFS="\t" {print $1,$2,$3,$1"."$2,"0","+";}' \
    $COMBINE_PEAKS_DIR/peaks_overlapped.bed > peaks_plus.bed
awk 'OFS="\t" {print $1,$2,$3,$1"."$2,"0","-";}' \
    $COMBINE_PEAKS_DIR/peaks_overlapped.bed > peaks_minus.bed

#find closest TSS for each peak on plus strand and on minus strand
#keep all ties
bedtools closest -s -D "a" -t "all" -a peaks_plus.bed \
    -b $ANNOTATION_DIR/all_transcript_TSSs_chr_uniq.bed > \
    peaks_plus_closestTSS_all.bed

bedtools closest -s -D "a" -t "all" -a peaks_minus.bed \
    -b $ANNOTATION_DIR/all_transcript_TSSs_chr_uniq.bed > \
    peaks_minus_closestTSS_all.bed

#take minimum of TSSs col8 and maximum of TES col9
bedtools merge -c 4,5,6,7,8,9,10,11,12,13 \
    -o "distinct,distinct,distinct,distinct,min,max,distinct,distinct,\
distinct,absmin" -i peaks_plus_closestTSS_all.bed > \
    peaks_plus_closestTSS_merged.bed

bedtools merge -c 4,5,6,7,8,9,10,11,12,13 \
    -o "distinct,distinct,distinct,distinct,min,max,distinct,distinct,\
distinct,absmin" -i peaks_minus_closestTSS_all.bed > \
    peaks_minus_closestTSS_merged.bed

#comma in col4: multiple peaks
#comma in col10: multiple genes
#only want to keep peaks with no comma in col4 and no comma in col10
#no comma in col10 because we only want to have one gene (multiple TSSs from 
#one gene are ok)
#what about comma in col4 but no comma in col10? then peak overlaps other peak 
#and both overlap a TSS from same gene
#what about comma in col 4 and no comma in col10: overlapping peaks that 
#overlap a TSS from same gene: discard (for now)

#take only peaks that have closest TSSs from same gene and that don't 
#intersect with other peaks that have closest TSSs from same gene
awk 'OFS="\t" {if($10!~"," && $4!~","){print $0}}' \
    peaks_plus_closestTSS_merged.bed > peaks_plus_closestTSS_filtered.bed

awk 'OFS="\t" {if($10!~"," && $4!~","){print $0}}' \
    peaks_minus_closestTSS_merged.bed > peaks_minus_closestTSS_filtered.bed

#no overlapping peaks anymore and no peaks that overlap TSSs from different 
#genes
#multiple TSSs from same gene is okay, smallest distance is kept


#for each peak: write peak, then closest TSS from plus strand, 
#then closest TSS from minus strand
#now for all peaks: write closest TSS from plus strand to left and closest TSS 
#on minus strand to right part of file
#here one TSS from plus strand and one from minus strand
#for each peak in PLUS combine with MINUS
#all peaks have plus strand here but strand of closest TSSs is correct, always 
#one plus and one minus, first one plus, second one minus
join -t $'\t' -j 4 \
    -o "1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.1,2.2,2.3,\
2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13" \
    <( sort -k4,4 peaks_plus_closestTSS_filtered.bed) \
    <(sort -k4,4 peaks_minus_closestTSS_filtered.bed) | \
    sort -k1,1 -k2,2n > peaks_closestTSS.bed

#peaks_plus_* files: col1-4 correct, col5 is 0, col6 is +
#peaks_minus_* files: col1-4 correct, col5 is 0, col6 is -
#this means col5 and 6 in these files don't mean anything
#1.5,1.6,2.5,2.6 should be ignored: col5,6,18,19

#needed for join
awk 'OFS="\t" {print $1"."$2,$0}' $COMBINE_PEAKS_DIR/peaks_overlapped.bed > \
    peaks_identifier.bed


# get unidirectional peaks
#one TSS on one strand only

#direction of TSS is fw
#distance 0
awk 'OFS="\t" {if($13==0 && $26!=0){print $0}}' peaks_closestTSS.bed > \
    peaks_closestTSS_unidirectional_plus.bed
awk 'OFS="\t" {if($13!=0 && $26==0){print $0}}' peaks_closestTSS.bed > \
    peaks_closestTSS_unidirectional_minus.bed

#keep peak name 1.5 and idenfier 2.4 here
#name 1.5 kept in last column, then known from which cell types peak came
#2.10 gene name, 2.11 gene type, 2.4 identifier,2.6 is strand that had closest
#TSS
#1.1 peak identifier, lose here, 1.2 chr, 1.3 start, 1.4 end, 1.5 name, 
#1.6 summit
join -t $'\t' -1 1 -2 4 -o 1.2,1.3,1.4,2.4,1.6,2.6,1.5,2.10,2.11 \
    <( sort -k1,1 peaks_identifier.bed)  \
    <( sort -k4,4 peaks_closestTSS_unidirectional_plus.bed) | \
    sort -k1,1 -k2,2n > peaks_unidirectional_plus.bed

#2.23 gene name, 2.24 gene type, 2.4 identifier,2.19 is strand that had 
#closest TSS
#1.1 peak identifier, lose here, 1.2 chr, 1.3 start, 1.4 end, 1.5 name, 
#1.6 summit
join -t $'\t' -1 1 -2 4 -o 1.2,1.3,1.4,2.4,1.6,2.19,1.5,2.23,2.24 \
    <( sort -k1,1 peaks_identifier.bed)  \
    <( sort -k4,4 peaks_closestTSS_unidirectional_minus.bed) | \
    sort -k1,1 -k2,2n > peaks_unidirectional_minus.bed

cat peaks_unidirectional_plus.bed peaks_unidirectional_minus.bed | \
    sort -k1,1 -k2,2n | awk 'OFS="\t" {print $0,".",".";}' > \
    peaks_unidirectional.bed


# enhancer candidates, no overlapping TSS on any strand!!!
#here no TSS pcg or lnc within 0 bp, that means overlapping
#strand information does not matter

##distance 0
awk 'OFS="\t" {if($13!=0 && $26!=0){print $0}}' peaks_closestTSS.bed > \
    peaks_closestTSS_nooverlappingTSS.bed

# 2.4 identifier,2.6 is strand that had closest TSS, has no meaning here
#1.1 peak identifier, lose here, 1.2 chr, 1.3 start, 1.4 end, 
#1.5 name, lose here, 1.6 summit
#remove strand information
join -t $'\t' -1 1 -2 4 -o 1.2,1.3,1.4,2.4,1.6,2.6,1.5 \
    <( sort -k1,1 peaks_identifier.bed)  \
    <( sort -k4,4 peaks_closestTSS_nooverlappingTSS.bed) | \
    sort -k1,1 -k2,2n | \
    awk 'OFS="\t" {print $1,$2,$3,$4,$5,".",$7,".",".",".",".";}' > \
    peaks_enhancers.bed

#divide enhancer candidates into intergenic and intragenic
bedtools intersect -v -wa -a peaks_enhancers.bed \
    -b $ANNOTATION_DIR/all_transcripts_chr.bed > peaks_enhancers_intergenic.bed


#strand, gene, type is from overlapping transcript
#strand is col17
bedtools intersect -wa -wb -a peaks_enhancers.bed \
    -b $ANNOTATION_DIR/all_transcripts_chr.bed | cut -f1-5,7,15-17 | \
    bedtools merge -c 4,5,6,7,8,9 -o distinct -i stdin | \
    awk 'OFS="\t" {print $1,$2,$3,$4,$5,$9,$6,".",".",$7,$8;}' > \
    peaks_enhancers_intragenic.bed


# bidirectional

#two overlapping TSSs, one frome each strand
#distance 0
awk 'OFS="\t" {if($13==0 && $26==0){print $0}}' peaks_closestTSS.bed > \
    peaks_closestTSS_bidirectional.bed

#filter out "bidirectional peaks" where TSS from plus strand has smaller 
#coordinate than TSS from minus strand and the genes are overlapping kind of
#$6 is always plus
#TSS from plus strand is $8, TSS from minus strand is $21
awk 'OFS="\t" {if($8>=$21){print $0}}' peaks_closestTSS_bidirectional.bed > \
    peaks_closestTSS_bidirectional_filtered.bed

#one TSS from plus strand and one TSS from minus strand
#2.10 gene name, 2.11 gene type, 2.4 identifier,2.6 is strand that had 
#closest TSS
#2.23 gene name, 2.24 gene type,2.19 strand
#1.1 peak identifier, lose here, 1.2 chr, 1.3 start, 1.4 end, 1.5 name, 
#1.6 summit
#TSSs from both strands, write "."
join -t $'\t' -1 1 -2 4 \
    -o 1.2,1.3,1.4,2.4,1.6,1.5,2.10,2.11,2.23,2.24 \
    <( sort -k1,1 peaks_identifier.bed)  \
    <( sort -k4,4 peaks_closestTSS_bidirectional_filtered.bed) | \
    sort -k1,1 -k2,2n | \
    awk 'OFS="\t" {print $1,$2,$3,$4,$5,".",$6,$7,$8,$9,$10;}'> \
    peaks_bidirectional.bed


cat peaks_unidirectional.bed peaks_bidirectional.bed \
    peaks_enhancers_intragenic.bed peaks_enhancers_intergenic.bed | \
    sort -k1,1 -k2,2n > peaks_allclassified.bed

cat peaks_unidirectional.bed peaks_bidirectional.bed | \
    sort -k1,1 -k2,2n > peaks_allclassified_promoters.bed

cat peaks_enhancers_intragenic.bed peaks_enhancers_intergenic.bed | \
    sort -k1,1 -k2,2n > peaks_allclassified_enhancers.bed
#FINAL SET OF OPEN CHROMATIN REGIONS

Rscript $SCRIPT_DIR/open_regions_subset/plot_classified_open_regions_widths.r 


#for decision about distance from open region and closest TSS
#I chose distance 0 (overlapping TSS)
#write distances from peaks to TSSs to file and plot
cut -f13 peaks_closestTSS.bed > peaks_closestTSS_plus_dist.txt
cut -f26 peaks_closestTSS.bed > peaks_closestTSS_minus_dist.txt

Rscript $SCRIPT_DIR/open_regions_subset/plot_closestTSS_dist.r 


#check distances to neighbor regions
#to decide how long feature regions are
#I remove overlapping feature regions later
paste <(head -n -1 peaks_allclassified.bed) \
    <(tail -n +2 peaks_allclassified.bed) > peaks_allclassified_neighbors.bed

Rscript $SCRIPT_DIR/open_regions_subset/plot_neighbor_dist.r 


exit


##########

#all peak files have:
#chr start end identifier summit strand name gene1 type1 gene2 type2

#UNIDIRECTIONAL: one overlapping TSS (any TSS) on plus or minus strand
#name gene type . .
#subtype overlapping pc TSS (peaks_unidirectional_pcg.bed)

#ENHANCER CANDIDATES: no overlapping TSS (any TSS)
#intragenic: name . . gene type
#intergenic: name . . . .
#subtype intragenic overlapping only one transcript 
#(peaks_enhancers_intragenic_onlyonetranscript.bed)

#BIDIRECTIONAL: two overlapping TSSs (any TSS)
#name gene1 type1 gene2 type2
#subtype bidirectional and at least one is pc TSS (peaks_bidirectional_pcg.bed)

##########

