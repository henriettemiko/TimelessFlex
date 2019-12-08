# TimelessFlex - A flexible framework for investigating chromatin state trajectories

by Henriette Miko (henriette.miko@mdc-berlin.de) at [Ohler lab](
https://github.com/ohlerlab) at BIMSB/MDC, December 2019


[ATAC](./ATAC): processing of ATAC-seq data and peak calling with MACS2 and IDR

[ChIP](./ChIP): processing of ChIP-seq data for histone modifications H3K27ac/me3 and H3K4me1/3 and peak calling with JAMM

[RNA](./RNA): processing of RNA-seq data with RSEM

[HiC](./ChIP): combining of processed Hi-C data

[open_regions_fullset](./open_regions_fullset): scripts for a linear time-course of 5 time points (D0, D2, D5, D7, D10)

[open_regions_split](./open_regions_split): scripts for a branched set with 5 time points of linear differentiation and a split from second time point (D0, D2, D5, D7, D10, D8 (split from D5))

[open_regions_subset](./open_regions_subset): scripts for a linear time course with 4 time points, for which Hi-C data is available (D0, D2, D5, D10)




This framework adapts and extends [Timeless](https://github.com/mahmoudibrahim/timeless). 


The [Bayes Net Toolbox for Matlab](https://github.com/bayesnet/bnt) is required.


