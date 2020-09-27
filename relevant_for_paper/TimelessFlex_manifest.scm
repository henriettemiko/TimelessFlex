;; This is the custom BIMSB module containing definitions for non-free
;; software, such as the old bowtie.
(use-modules (bimsb packages bioinformatics-nonfree))

(define packages
  (list "glibc-locales" ; for locales
        "nss-certs"     ; for certificates

        "wget"
	"sra-tools"			;fastq-dump, prefetch
	"fastqc"
	"trim-galore"
	"cutadapt"
	"bowtie1"			; non-free!
	"bowtie@2"
	"macs"
	"samtools"
	"bedtools"
	"bedops"
	"java-picard"		  ; NOTE: there are no wrapper scripts
	"rsem"
	"deeptools"
	"idr"			;idrCode/batch-consistency-analysis.r)
	"python2"
        "perl"
	
	"bash"
	"grep"
	"gawk"
	"gzip"
	"coreutils"

	;; R with packages:
	"r-minimal"
	"r-limma"
	"r-mclust"
	"r-psych"
	"r-ggplot2"
	"r-reshape2"
	"r-factoextra"
	"r-nbclust"
	"r-signal"
	"r-pheatmap"
	"r-biomart"
	"r-gprofiler"
	"r-deseq2"
	"r-tximport"
	"r-gplots"
	"r-rgraphviz"
	"r-rcolorbrewer"
	"r-dplyr"
	"r-plyr"

	"rseqc"
	"htseq"

	"jamm"
	"kentutils"
	"kentutils-nonfree"
	"idr-legacy"))

(specifications->manifest packages)
