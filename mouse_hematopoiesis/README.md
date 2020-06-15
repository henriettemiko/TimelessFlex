
# TimelessFlex - A flexible framework for investigating chromatin state trajectories

by Henriette Miko (henriette.miko@mdc-berlin.de), [Ohler lab](
https://github.com/ohlerlab) at BIMSB/MDC, June 15, 2020

## Mouse hematopoiesis data

All code for processing and analyzing data from mouse hematopoiesis at 6 time points (CMP, MEP, EryA, GMP, Granu, Mono) can be found here.
The structure of this directory is as follows:

[ATAC](./ATAC): scripts for processing of ATAC-seq data and peak calling with MACS2 and IDR

[ChIP](./ChIP): scripts for processing of ChIP-seq data for histone modifications H3K27ac/me3 and H3K4me1/3 and peak calling with JAMM

[HiC](./HiC): script for combining of processed Hi-C data

[RNA](./RNA): scripts for processing of RNA-seq data with RSEM

[open_regions_fullset](./open_regions_fullset): scripts for a linear time-course of 5 time points (D0 -> D2 -> D5 -> D7 -> D10)

[open_regions_split](./open_regions_split): scripts for a branched set with 5 time points of linear differentiation and a split from second time point (D0 -> D2 -> D5 -> D7 -> D10 and D8 split from D5)

[open_regions_subset](./open_regions_subset): scripts for a linear time course with 4 t
