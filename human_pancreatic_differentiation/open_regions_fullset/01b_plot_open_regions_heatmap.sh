#!/bin/bash

##########
#name:          01b_plot_open_regions_heatmap.sh
#description:   plot heatmaps for open regions
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 10, 2019
##########


source ../set_variables_hg19.sh

CLASSIFY_OPEN_REGIONS_DIR=$OPEN_REGIONS_DIR_FULL/classify_open_regions

HEATMAPS_DIR=$OPEN_REGIONS_DIR_FULL/heatmaps
mkdir -p $HEATMAPS_DIR

qsub -l os=centos7 -V -cwd -j y \
    -o $HEATMAPS_DIR/plot_profile_heatmap_H3K4me1.txt \
    -l mem_free=10G,h_vmem=10G \
    $SCRIPT_DIR/open_regions_fullset/plot_profile_heatmap.sh \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_unidirectional.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_bidirectional.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_enhancers_intergenic.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_enhancers_intragenic.bed \
    $OUTPUT_DIR $HEATMAPS_DIR H3K4me1

qsub -l os=centos7 -V -cwd -j y \
    -o $HEATMAPS_DIR/plot_profile_heatmap_H3K27ac.txt \
    -l mem_free=10G,h_vmem=10G \
    $SCRIPT_DIR/open_regions_fullset/plot_profile_heatmap.sh \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_unidirectional.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_bidirectional.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_enhancers_intergenic.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_enhancers_intragenic.bed \
    $OUTPUT_DIR $HEATMAPS_DIR H3K27ac

qsub -l os=centos7 -V -cwd -j y \
    -o $HEATMAPS_DIR/plot_profile_heatmap_H3K4me3.txt \
    -l mem_free=10G,h_vmem=10G \
    $SCRIPT_DIR/open_regions_fullset/plot_profile_heatmap.sh \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_unidirectional.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_bidirectional.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_enhancers_intergenic.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_enhancers_intragenic.bed \
    $OUTPUT_DIR $HEATMAPS_DIR H3K4me3

qsub -l os=centos7 -V -cwd -j y \
    -o $HEATMAPS_DIR/plot_profile_heatmap_H3K27me3.txt \
    -l mem_free=10G,h_vmem=10G \
    $SCRIPT_DIR/open_regions_fullset/plot_profile_heatmap.sh \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_unidirectional.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_bidirectional.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_enhancers_intergenic.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_enhancers_intragenic.bed \
    $OUTPUT_DIR $HEATMAPS_DIR H3K27me3

qsub -l os=centos7 -V -cwd -j y \
    -o $HEATMAPS_DIR/plot_profile_heatmap_ATAC.txt \
    -l mem_free=10G,h_vmem=10G \
    $SCRIPT_DIR/open_regions_fullset/plot_profile_heatmap.sh \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_unidirectional.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_bidirectional.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_enhancers_intergenic.bed \
    $CLASSIFY_OPEN_REGIONS_DIR/peaks_enhancers_intragenic.bed \
    $OUTPUT_DIR $HEATMAPS_DIR ATAC


exit

