#!/bin/bash

##########
#name:          call_MACS2_ATAC.sh
#description:   call peaks with MACS2 from ATAC-seq data
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 8, 2019
##########


MACS2_DIR=$1
TIME=$2
OUTPUT_DIR=$3
MACS2_ORGANISM=$4
CHR_SIZES=$5


echo "#####################"
echo $0 started on `hostname` at `date` with parameters $*
echo "current time point: $TIME"
echo "#####################"


#create list of bed files as input for macs2
BED_FILES=($OUTPUT_DIR/ATAC/$TIME/mapping/*_filtered_final_1bp.bed)
NUMBER_BED_FILES="${#BED_FILES[@]}"

if [[ $NUMBER_BED_FILES == 1 ]]
then

    NAME=$(basename $BED_FILES _filtered_final_1bp.bed)

    #if there is no replicate call MACS2 with strict value q=0.01
    macs2 callpeak -t "${BED_FILES[@]}" -f BED -g $MACS2_ORGANISM \
        --outdir ${MACS2_DIR} -n ${TIME}_${NAME}_q001 --nomodel --shift -100 \
        --extsize 200 --keep-dup all -q 0.01


    #store summit position (not offset!) in last column
    awk 'BEGIN{OFS="\t"} {print $0,($2+$10)}' \
        ${MACS2_DIR}/${TIME}_${NAME}_q001_peaks.narrowPeak > \
        ${MACS2_DIR}/${TIME}_${NAME}_q001_peaks_summitpos_final.narrowPeak


elif [[ $NUMBER_BED_FILES == 2 ]]
then

    #call peaks with MACS2 for pooled

    macs2 callpeak -t "${BED_FILES[@]}" -f BED -g $MACS2_ORGANISM \
        --outdir ${MACS2_DIR} -n ${TIME}_pooled_p005 --nomodel --shift -100 \
        --extsize 200 --keep-dup all -p 0.05

    #call MACS2 with lose value p=0.05 for each replicate
    #and then IDR (deprecated version)

    for i in "${BED_FILES[@]}"
    do
        echo $i

        NAME=$(basename $i _filtered_final_1bp.bed)
        macs2 callpeak -t $i -f BED -g $MACS2_ORGANISM --outdir ${MACS2_DIR} \
            -n ${TIME}_${NAME}_p005 --nomodel --shift -100 --extsize 200 \
            --keep-dup all -p 0.05

        #for IDR take top 100000 peaks, sort according to p-value
        sort -k 8nr,8nr ${MACS2_DIR}/${TIME}_${NAME}_p005_peaks.narrowPeak | \
            head -n 100000 > \
            ${MACS2_DIR}/${TIME}_${NAME}_p005_top100000.narrowPeak

    done

    #idrCode/batch-consistency-analysis.r
    #https://sites.google.com/site/anshulkundaje/projects/idr/deprecated
    #take top 100.000 peaks from rep1 and rep2, sort according to p.value, 
    #run IDR
    IDR_FILES=(${MACS2_DIR}/*top100000.narrowPeak)
    cp $CHR_SIZES genome_table.txt


    Rscript ${GUIX_PROFILE}/share/idr-legacy/batch-consistency-analysis.r \
        "${IDR_FILES[@]}" -1 ${MACS2_DIR}/idr 0 F p.value

    THRESHOLD=$(awk '$11 <= 0.01 {print $0}' \
        ${MACS2_DIR}/idr-overlapped-peaks.txt | wc -l)
    echo $THRESHOLD

    #select top peaks from pooled file
    sort -k 8nr,8nr ${MACS2_DIR}/${TIME}_pooled_p005_peaks.narrowPeak | \
        head -n $THRESHOLD > \
        ${MACS2_DIR}/${TIME}_pooled_p005_peaks_IDR.narrowPeak

    #store summit position (not offset!) in last column
    awk 'BEGIN{OFS="\t"} {print $0,($2+$10)}' \
        ${MACS2_DIR}/${TIME}_pooled_p005_peaks_IDR.narrowPeak > \
        ${MACS2_DIR}/${TIME}_pooled_p005_peaks_IDR_summitpos_final.narrowPeak

else

    echo "number of bedfiles ís $NUM_BED_FILES but must be 1 or 2"
fi


exit

