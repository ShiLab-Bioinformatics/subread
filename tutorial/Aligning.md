# Tutorial of using the Subread aligner for read mapping

Rsubread provides a hughly accurate and efficient read aligner called ``subread-align``. This tutorial demonstrates a workflow for mapping HTS (high-thoughtput sequencing) data using our read aligner. 

## Prerequisites

You need an R environment installed on your computer. It can be a Windows PC, a Linux server or a macOS computer, but it must be a 64-bit system. A computer built in or after 2012 is nearly certain a 64-bit system. Your computer also needs at least 100GB of disk space and 16GB of RAM. 

If you are using Windows or macOS, R studio is recommended for that it has everything ready for you installing our package. If you run Linux, you need to have a modern version of gcc installed with libz-devel.

## Installing Rsubread

After entrining the R environment, you can enter commands. The commands for installing the latest version of R package is:

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("Rsubread")
```
You will be asked a few questions about where to download the related files and (sometime) whether you want to upgrade related packages. After you select the best server for downloading and ugrading "all" related packages, it takes a while to get everything ready. No error should be reported in all the installation process but it is not uncommon to have some errors. You may report the error through [the supporting site of Bioconductor](https://support.bioconductor.org/).

After successful installation of Rsubread, you will be able to load the package into R:

```R
library(Rsubread)
sessionInfo()
```
You will be able to see Rsubread in the output from sessionInfo().

## Downloading the data
We provide a set of publically available data.


## Building an index
Read aligners need to pre-process the reference genome sequences to efficiently map reads to them. This step is usually called index building. A function is provided in Rsubread for building an index on the given reference genome. 
```R
buildindex("hg38-index", "hg38.fasta.gz", indexSplit=FALSE, gappedIndex=TRUE)
```

If your computer has 32GB or more memory, you can set ``gappedIndex=FALSE`` to build a full index. A full index uses around 3 times more memory than a gapped index, but can largely accelerate read mapping and deliver higher mapping quality. 

## Read mapping
The Rsubread package provides two functions for read mapping: ``align()`` and ``subjunc()``. The ``align()`` function can map both DNA-seq and RNA-seq data while the ``subjunc()`` function is dedicated for mapping RNA-seq reads, but with all exons detected in reads that contain exon-exon junctions.

Although it only maps at most one exon in each read, the ``align()`` function is usually good enough for gene Differential-Expression (DE) analyses. The ``subjunc()`` function is for specific analyses that involve calling of known/de novo exon-exon junctions.
```R
align.sum <- align(
  "hg38-index", "SEQC-A_R1.fastq.gz", readfile2="SEQC-A_R2.fastq.gz", output_file="SEQC-A-align.bam",
  nthreads=4, useAnnotation=TRUE, annot.inbuilt="hg38"
)

junc.sum <- subjunc(
  "hg38-index", "SEQC-A_R1.fastq.gz", readfile2="SEQC-A_R2.fastq.gz", output_file="SEQC-A-align.bam",
  nthreads=4, useAnnotation=TRUE, annot.inbuilt="hg38"
)

print(align.sum)
print(junc.sum)
```



## 
