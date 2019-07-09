#!/bin/bash

##########
#name:          run_all.sh
#description:   runs all scripts
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 9, 2019
##########


source set_variables_hg19.sh


#ChIP-seq data
source ChIP/01a_process_ChIP_files.sh
source ChIP/01b_call_plotFingerprint_ChIP.sh
source ChIP/02a_compute_ChIP_JAMM_peaks.sh
source ChIP/02b_plot_JAMM_peak_widths.sh


#ATAC-seq data
source ATAC/01a_process_ATAC_files.sh
source ATAC/01b_call_plotFingerprint_ATAC.sh
source ATAC/02a_compute_ATAC_MACS2_peaks.sh
source ATAC/02b_plot_ATAC_peak_widths.sh
source ATAC/01c_compute_fragment_lengths_ATAC.sh


#RNA-seq data
source RNA/01a_process_RNA.sh



exit

