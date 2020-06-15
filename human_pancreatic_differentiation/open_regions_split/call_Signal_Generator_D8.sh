#!/bin/bash

##########
#name:          call_Signal_Generator.sh
#description:   call JAMM Signal Generator
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 10, 2019
##########


echo $0 started on `hostname` at `date` with parameters $*

MARK=$1
TIME=$2
SIGNAL_GENERATOR_DIR=$3
OUTPUT_DIR=$4
REGIONS=$5
CHR_SIZES=$6

echo "#####################"
echo "Started processing:"
echo "open regions: $REGIONS"
echo "current histone marks: $MARK"
echo "current time point: $TIME"
echo "#####################"


JAMM_OUT_FILE=$OUTPUT_DIR/ChIP/$MARK/$TIME/JAMM_peaks/call_JAMM.txt

#read in fragment lengths from JAMM_OUT_FILE
frags=$(grep "Fragment Length:" "${JAMM_OUT_FILE}" | cut -d":" -f2);
echo $frags
fragments=""; while read -r line; do
fragments="$fragments,$line"; done <<< "$frags"; 
fragments=${fragments#","};
echo $fragments;

#for D8 Input fragment length estimation failed, use 100

bash SignalGenerator.sh -s $OUTPUT_DIR/ChIP/$MARK/$TIME/mapping \
    -c $OUTPUT_DIR/ChIP/Input/$TIME/mapping -g $CHR_SIZES -r $REGIONS \
    -f "$fragments,100" -b 1 -n depth -o $SIGNAL_GENERATOR_DIR/$MARK/$TIME;


CUR_DIR=$PWD
cd $CUR_DIR

for f in $SIGNAL_GENERATOR_DIR/$MARK/$TIME/signal/
do 
    echo $f
    cd $f
    rm chrM.bedGraph
    cat *.bedGraph | awk '{print $1"\t"$2"\t"$3"\t"".""\t"$4}' | \
        sort-bed-typical - > ../all.bg
    cd $CUR_DIR
done

#get maximum counts over regions
for f in $SIGNAL_GENERATOR_DIR/$MARK/$TIME/
do
    cd $f
    echo $f
    sort-bed-typical $REGIONS | bedmap-typical --delim "\t" --max - all.bg | \
        sed 's/NAN/0/' > max_counts.txt
    cd $CUR_DIR
done

#get mean counts over regions
for f in $SIGNAL_GENERATOR_DIR/$MARK/$TIME/
do
    cd $f
    echo $f
    sort-bed-typical $REGIONS | bedmap-typical --delim "\t" --mean - all.bg | \
        sed 's/NAN/0/' > mean_counts.txt
    cd $CUR_DIR

done


exit

