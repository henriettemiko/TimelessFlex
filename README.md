# TimelessFlex - A flexible framework for investigating chromatin state trajectories

by Henriette Miko (henriette.miko@mdc-berlin.de) at [Ohler lab](
https://github.com/ohlerlab) at BIMSB/MDC, December 8, 2019




### Framework description

TimelessFlex is a flexible framework for investigating chromatin state trajectories 
during linear and branched tissue differentiation. 
This framework adapts and extends [Timeless](https://github.com/mahmoudibrahim/timeless), a Bayesian network for co-clustering of multiple histone modifications at promoter and enhancer feature regions. 
TimelessFlex is flexible depending on the available genomic data. The basic required data is ChIP-seq for histone modifications with at least three time points and a set of regions of interest.


The framework consists of three parts:
1. Regions definition: 
   - combine ATAC-seq peaks over multiple time points into one set of open chromatin regions
   - define promoters and enhancers from open chromatin regions
   - if Hi-C data is available, assign promoters and enhancers to Hi-C interaction pairs
   - define feature regions around promoters and enhancer   
2. Clustering:
   - compute histone mark signals over feature regions
   - cluster histone mark signals with Bayesian network
3. Validation and interpretation of clusters
   - validate clusters with genomic data not used in clustering
   - functional analyses of clusters



### Implementation details

This framework is written in R and bash and it makes use of multiple publicly available bioinformatics software. The clustering is performed in MATLAB with the [Bayes Net Toolbox (BNT)](https://github.com/bayesnet/bnt).



### Repository structure

[ATAC](./ATAC): processing of ATAC-seq data and peak calling with MACS2 and IDR

[ChIP](./ChIP): processing of ChIP-seq data for histone modifications H3K27ac/me3 and H3K4me1/3 and peak calling with JAMM

[RNA](./RNA): processing of RNA-seq data with RSEM

[HiC](./ChIP): combining of processed Hi-C data

[open_regions_fullset](./open_regions_fullset): scripts for a linear time-course of 5 time points (D0 -> D2 -> D5 -> D7 -> D10)

[open_regions_split](./open_regions_split): scripts for a branched set with 5 time points of linear differentiation and a split from second time point (D0 -> D2 -> D5 -> D7 -> D10 and D8 split from D5)

[open_regions_subset](./open_regions_subset): scripts for a linear time course with 4 time points, for which Hi-C data is available (D0 -> D2 -> D5 -> D10)

[open_regions_subset/open_regions_pairs](./open_regions_subset/open_regions_pairs): scripts for Hi-C interaction pairs
