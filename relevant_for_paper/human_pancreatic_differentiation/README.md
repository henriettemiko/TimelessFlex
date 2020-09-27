

# TimelessFlex - A flexible framework for investigating chromatin state trajectories

by Henriette Miko (henriette.miko@mdc-berlin.de), [Ohler lab](
https://github.com/ohlerlab) at BIMSB/MDC, 2020

## Human pancreatic differentiation data

All code for processing and analyzing data from human pancreatic differentiation at 6 time points D0, D2, D5, D7, D10 and D8 can be found here.
The structure of this directory is as follows:

[ATAC](./ATAC): scripts for processing of ATAC-seq data and peak calling with MACS2 and IDR

[ChIP](./ChIP): scripts for processing of ChIP-seq data for histone modifications H3K27ac/me3 and H3K4me1/3 and peak calling with JAMM

[HiC](./HiC): script for combining of processed Hi-C data

[RNA](./RNA): scripts for processing of RNA-seq data with RSEM

[open_regions_subset](./open_regions_subset): scripts for a linear time course with 4 time points, for which Hi-C data is available (D0 -> D2 -> D5 -> D10)

[open_regions_subset/open_regions_pairs](./open_regions_subset/open_regions_pairs): scripts for Hi-C interaction pairs (promoter-enhancer, promoter-promoter and enhancer-enhancer pairs)
