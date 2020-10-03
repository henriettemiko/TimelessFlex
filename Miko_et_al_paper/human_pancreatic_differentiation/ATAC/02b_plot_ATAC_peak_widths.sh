#!/bin/bash

##########
#name:          02b_plot_ATAC_peak_widths.sh
#description:   plot widths of peaks called with MACS2
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


source ../set_variables_hg19.sh

QUALITY_DIR=$OUTPUT_DIR/quality_ATAC

mkdir -p $QUALITY_DIR
cd $QUALITY_DIR


qsub -hold_jid "call_MACS2_ATAC_*" -N plot_peaks \
    -V -j y \
    -o $QUALITY_DIR/call_plot_peaks.txt -cwd -pe smp 1 \
    -l mem_free=10G,h_vmem=10G -l h_rt=24:00:00 \
    $SCRIPT_DIR/ATAC/call_plot_ATAC_MACS2_peak_widths.sh \
    $OUTPUT_DIR $SCRIPT_DIR $QUALITY_DIR


exit

