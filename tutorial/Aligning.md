## Tutorial of using the Subread aligner for read mapping

Rsubread provides a hughly accurate and efficient read aligner called ``subread-align``. This tutorial demonstrates a workflow for mapping HTS (high-thoughtput sequencing) data using our read aligner. 

# Prerequisites

You need an R environment installed on your computer. It can be a Windows PC, a Linux server or a macOS computer, but it must be a 64-bit system. A computer built in or after 2012 is nearly certain a 64-bit system. Your computer also needs at least 100GB of disk space and 16GB of RAM. 

If you are using Windows or macOS, R studio is recommended for that it has everything ready for you installing our package. If you run Linux, you need to have a modern version of gcc installed with libz-devel.

# Installing Rsubread

After entrining the R environment, you can enter commands. The commands for installing the latest version of R package is:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rsubread")
```
You will be asked a few questions about where to download the related files and (sometime) if you want to upgrade related packages. After you select the best server for downloading and ugrading "all" related packages, it takes a while to get everything ready. No error should be reported in all the installation process but it is not uncommon to have some errors. You may report the error through [the supporting site of Bioconductor](https://support.bioconductor.org/).

After successful installation of Rsubread, you will be able to load the package into R:

```R
library(Rsubread)
sessionInfo()
```
You will be able to see Rsubread in the output from sessionInfo().
