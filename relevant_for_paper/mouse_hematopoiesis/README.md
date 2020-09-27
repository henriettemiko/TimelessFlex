
# TimelessFlex - A flexible framework for investigating chromatin state trajectories

by Henriette Miko (henriette.miko@mdc-berlin.de), [Ohler lab](
https://github.com/ohlerlab) at BIMSB/MDC, June 15, 2020

## Mouse hematopoiesis data

All code for processing and analyzing data from mouse hematopoiesis at 6 time points (CMP, MEP, EryA, GMP, Granu, Mono) can be found here.
The structure of this directory is as follows:

[ATAC](./ATAC): scripts for processing of ATAC-seq data and peak calling with MACS2

[ChIP](./ChIP): scripts for processing of ChIP-seq data for histone modifications H3K27ac and H3K4me1/2/3 and peak calling with JAMM

[open_regions](./open_regions): scripts for a branched time-course of 6 time points (CMP -> MEP -> EryA, CMP -> GMP -> Granu/Mono)
