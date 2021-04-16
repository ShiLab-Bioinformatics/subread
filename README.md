# Subread
The Subread software package is a tool kit for processing next-gen sequencing data. It includes Subread aligner, Subjunc exon-exon junction detector and featureCounts read summarization program.

## Installation
The latest releases can be downloaded from the [release page](https://github.com/ShiLab-Bioinformatics/subread/releases).

### Installation from a binary package
The easist way to installing Subread on Linux, Windows and macOS is to directly download the binary packages on our Release page. Simply decompress the package and the programs will be in the "/bin" directory.

### Installation in the R environment
We also provide an R version of our package, [Rsubread](http://bioconductor.org/packages/Rsubread), on Bioconductor. You can follow the instructions on the Rsubread webpage to install it in R.

### Installation from the source code
An experienced user may also try building the binary programs from source code. To this end, some programs and libraries are necessary.

1. A C language compiler. It can be gcc or clang or anything. Intel CC should work well but we have not tried it on our source code.
2. Libraries including zlib (for gzip), libpthread (for multi-threading) and libm (for math). We tried to reduce our dependency as much as possible so no fancy libraries are needed.
3. GNUmake and Shell.

If you use Windows, you may consider to install [Mingw-w64](http://mingw-w64.org/doku.php) for that it provides all the required programs and libraries in one place. You do not need Cygwin or the Linux subsystem in Windows for building Subread.

Compiling the source code is simple. 
```sh
$ cd src
$ make -f Makefile.Linux (for Linux)
$ make -f Makefile.MacOS (for macOS)
$ make -f Makefile.Windows (for Windows)
```
The executable programs will be moved to the "/bin" directory. Because we want to minimise the dependency to other packages, we do not use autoconf to generate the Makefiles. 

### Testing the installation
After installation, you may test the programs to see if it works.

The Subread package incorporates many small testcases that cover most of its functions. No matter if Subread is installed from the source code or a binary package, you can run "test_all.sh" in the "/test" directory. This assumes that you have a Shell program and a Python2 interpreter in PATH.

```sh
$ cd test
$ sh test_all.sh
```

## Usage
The usages of the programs in this package can be found in the users-guide in the "/doc" directory.

## Citation
We have published papers on our Subread/Subjunc read aligners and featureCounts read quantifiers.

1. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote, ***Y Liao, GK Smyth, W Shi***, Nucleic acids research, 2013 [PMID:23558742](https://pubmed.ncbi.nlm.nih.gov/23558742/)

2. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features, ***Y Liao, GK Smyth, W Shi***, Bioinformatics, 2014 [PMID:24227677](https://pubmed.ncbi.nlm.nih.gov/24227677/)

3. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads, ***Y Liao, GK Smyth, W Shi***,  Nucleic acids research, 2019 [PMID:30783653](https://pubmed.ncbi.nlm.nih.gov/30783653/)

## PhD projects
PhD projects are available for further development of the Subread package, including the development of new methods for analyzing single-cell sequencing data. For any inquiries, please contact [Prof Wei Shi](https://www.onjcri.org.au/about-us/wei-shi/).
