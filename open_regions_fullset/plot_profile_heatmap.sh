#!/bin/bash

##########
#name:          plot_profile_heatmap.sh
#description:   call computeMatrix and plotHeatmap from deeptools
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 10, 2019
##########


echo $0 started on `hostname` at `date` with parameters $*

PEAKS_UNI=$1
PEAKS_BI=$2
PEAKS_INTER=$3
PEAKS_INTRA=$4
OUTPUT_DIR=$5
HEATMAPS_DIR=$6
MARK=$7

cd $HEATMAPS_DIR

#extend the regions around summit so center of new regions is summit
awk 'BEGIN{OFS="\t"} {print $1,($5-2000),($5+2001),$4".ext",
$5,$6,$7,$8,$9,$10,$11}' $PEAKS_UNI > peaks_uni_2000.bed

awk 'BEGIN{OFS="\t"} {print $1,($5-2000),($5+2001),$4".ext",
$5,$6,$7,$8,$9,$10,$11}' $PEAKS_BI > peaks_bi_2000.bed

awk 'BEGIN{OFS="\t"} {print $1,($5-2000),($5+2001),$4".ext",
$5,$6,$7,$8,$9,$10,$11}' $PEAKS_INTER > peaks_inter_2000.bed

awk 'BEGIN{OFS="\t"} {print $1,($5-2000),($5+2001),$4".ext",
$5,$6,$7,$8,$9,$10,$11}' $PEAKS_INTRA > peaks_intra_2000.bed

if [ $MARK == "ATAC" ]
then
    FILES=($OUTPUT_DIR/ATAC/D*/bedgraphs/*.bigWig)
    echo "${FILES[@]}"
else 
    FILES=($OUTPUT_DIR/ChIP/$MARK/D*/bedgraphs/*.bigWig)
    echo "${FILES[@]}"
fi

computeMatrix reference-point --referencePoint center -a 2000 -b 2000 \
    -R "peaks_uni_2000.bed peaks_bi_2000.bed peaks_inter_2000.bed 
peaks_intra_2000.bed" -S "${FILES[@]}" -o matrix_${MARK}.gz -p 16 > \
    out.txt 2>&1

plotHeatmap -m matrix_${MARK}.gz -out heatmap_${MARK}.pdf --perGroup \
    --plotTitle "${MARK}" --startLabel "-2000" --endLabel "+2000" \
    --refPointLabel "summit" \
    --samplesLabel "D0 rep1" "D0 rep2" "D10 rep1" "D10 rep2" "D2 rep1" 
"D2 rep2" "D5 rep1" "D5 rep2" "D7 rep1" "D7 rep2"  \
    --regionsLabel "unidirectional" "bidirectional" "intergenic" "intragenic" \
    --legendLocation center-left

computeMatrix reference-point --referencePoint center -a 2000 -b 2000 \
    -R peaks_uni_2000.bed peaks_bi_2000.bed -S "${FILES[@]}" \
    -o matrix_${MARK}_prom.gz -p 16 > out_prom.txt 2>&1

computeMatrix reference-point --referencePoint center -a 2000 -b 2000 \
    -R peaks_inter_2000.bed peaks_intra_2000.bed -S "${FILES[@]}" \
    -o matrix_${MARK}_enh.gz -p 16 > out_enh.txt 2>&1

plotHeatmap -m matrix_${MARK}_prom.gz -out heatmap_${MARK}_prom.pdf \
    --perGroup --plotTitle "${MARK}" --startLabel "-2000" --endLabel "+2000" \
    --refPointLabel "summit" \
    --samplesLabel "D0 rep1" "D0 rep2" "D10 rep1" "D10 rep2" "D2 rep1" 
"D2 rep2" "D5 rep1" "D5 rep2" "D7 rep1" "D7 rep2"  \
    --regionsLabel "unidirectional" "bidirectional" \
    --legendLocation center-left

plotHeatmap -m matrix_${MARK}_enh.gz -out heatmap_${MARK}_enh.pdf --perGroup \
    --plotTitle "${MARK}" --startLabel "-2000" --endLabel "+2000" \
    --refPointLabel "summit" \
    --samplesLabel "D0 rep1" "D0 rep2" "D10 rep1" "D10 rep2" "D2 rep1" 
"D2 rep2" "D5 rep1" "D5 rep2" "D7 rep1" "D7 rep2"  \
    --regionsLabel "intergenic" "intragenic" --legendLocation center-left


exit

