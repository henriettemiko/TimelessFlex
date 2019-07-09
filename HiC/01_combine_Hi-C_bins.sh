#!/bin/bash 

##########
#name:          01_combine_Hi-C_bins.sh
#description:   combine Hi-C bins from multiple time points
#author:        Henriette Miko (henriette.miko@mdc-berlin.de)
#date:          July 9, 2019
##########


source ../set_variables_hg19.sh

cd $HIC_DIR
mkdir -p D0 D2 D5 D10


#Hi-C files before collapsing bins
#store timepoint information in files
#compute 10kb bin bedfiles for visualization in the browser

for file in $START_DIR/input/HiC/D0/D0.chr*; do

    t=$(basename $file ".BP.10000.txt.removed" | cut -d. -f1)
    c=$(basename $file ".BP.10000.txt.removed" | cut -d. -f2)
    echo $t
    echo $c

    awk -v c="$c" -v t="$t" 'OFS="\t" {
    print c,$1,($1+10000),"D0."c"."$1"."$2".left",".","D0\n"
    c,$2,($2+10000),"D0."c"."$1"."$2".right",".","D0"}' \
        <(tail -n +2 $file) | sort -k1,1 -k2,2n > D0/D0_$c.bed

done

cat D0/D0_chr*.bed | sort -k1,1 -k2,2n > D0/D0.bed


for file in $START_DIR/input/HiC/D2/D2.chr*; do

    t=$(basename $file ".BP.10000.txt.removed" | cut -d. -f1)
    c=$(basename $file ".BP.10000.txt.removed" | cut -d. -f2)
    echo $t
    echo $c

    awk -v c="$c" -v t="$t" 'OFS="\t" {
    print c,$1,($1+10000),"D2."c"."$1"."$2".left",".","D2\n"
    c,$2,($2+10000),"D2."c"."$1"."$2".right",".","D2"}' \
        <(tail -n +2 $file) | sort -k1,1 -k2,2n > D2/D2_$c.bed

done

cat D2/D2_chr*.bed | sort -k1,1 -k2,2n > D2/D2.bed


for file in $START_DIR/input/HiC/D5/D5.chr*; do

    t=$(basename $file ".BP.10000.txt.removed" | cut -d. -f1)
    c=$(basename $file ".BP.10000.txt.removed" | cut -d. -f2)
    echo $t
    echo $c
    awk -v c="$c" -v t="$t" 'OFS="\t" {
    print c,$1,($1+10000),"D5."c"."$1"."$2".left",".","D5\n"
    c,$2,($2+10000),"D5."c"."$1"."$2".right",".","D5"}' \
        <(tail -n +2 $file) | sort -k1,1 -k2,2n > D5/D5_$c.bed
done

cat D5/D5_chr*.bed | sort -k1,1 -k2,2n > D5/D5.bed


for file in $START_DIR/input/HiC/D10/D10.chr*; do

    t=$(basename $file ".BP.10000.txt.removed" | cut -d. -f1)
    c=$(basename $file ".BP.10000.txt.removed" | cut -d. -f2)
    echo $t
    echo $c

    awk -v c="$c" -v t="$t" 'OFS="\t" {
    print c,$1,($1+10000),"D10."c"."$1"."$2".left",".","D10\n"
    c,$2,($2+10000),"D10."c"."$1"."$2".right",".","D10"}' \
        <(tail -n +2 $file) | sort -k1,1 -k2,2n > D10/D10_$c.bed

done

cat D10/D10_chr*.bed | sort -k1,1 -k2,2n > D10/D10.bed


#combine Hi-C bins over time points
cat D0/D0.bed D2/D2.bed D5/D5.bed D10/D10.bed | sort -k1,1 -k2,2n > \
    all_bins_notunique.bed
#here a lot of bins can occur multiple times, multiple times in same timepoint and multiple times at different timepoints

#get unique bins
bedtools merge -c 4,5,6 -o distinct -d -10000 -i all_bins_notunique.bed > \
    all_bins_unique.bed
#6 cols: chr start end name(chr*.*.*.left) strand(".") timepoint 

echo "number of unique Hi-C bins:"
cut -f 1-3 all_bins_unique.bed | sort | uniq | wc -l


exit

