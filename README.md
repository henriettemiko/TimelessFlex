# TimelessFlex - A flexible framework for investigating chromatin state trajectories

by Henriette Miko (henriette.miko@mdc-berlin.de), [Ohler lab](
https://github.com/ohlerlab) at BIMSB/MDC, 2020


## Framework description

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
3. Validation and interpretation:
   - validate clusters with genomic data not used in clustering
   - functional analyses of clusters


## Implementation details

This framework is written in R and bash and it makes use of multiple publicly available bioinformatics software. The clustering is performed in MATLAB with the [Bayes Net Toolbox (BNT)](https://github.com/bayesnet/bnt).


## Repository structure

### Code for data from mouse hematopoiesis

All code for processing and analyzing data from mouse hematopoiesis at 6 time points (CMP, MEP, EryA, GMP, Granu, Mono) can be found in [mouse_hematopoiesis](./mouse_hematopoiesis).


### Code for data from human pancreatic differentiation

All code for processing and analyzing data from human pancreatic differentiation at 5 time points D0, D2, D5, D7 and D10 can be found in [human_pancreatic_differentiation](./human_pancreatic_differentiation).


### Relevant code for paper

All relevant code for Miko et al., "Inferring time-course chromatin states for promoter-enhancer pairs based on Hi-C data" can be found in [paper](./paper). This code is basically a subset of [mouse_hematopoiesis](./mouse_hematopoiesis) and [human_pancreatic_differentiation](./human_pancreatic_differentiation).
