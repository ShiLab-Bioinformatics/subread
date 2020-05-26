# Subread
The Subread software package is a tool kit for processing next-gen sequencing data. It includes Subread aligner, Subjunc exon-exon junction detector and featureCounts read summarization program.

# Installation
## Installation from a binary package
The easist way to installing Subread on Linux, Windows and macOS is to directly download the binary packages on our Release page. Simply decompress the package and the programs will be in the "/bin" directory.

## Installation from the source code
An experienced user may also try building the binary programs from source code. To this end, some programs and libraries are necessary.

1. A C language compiler. It can be gcc or clang or anything. Intel CC should work well but we have not tried it on our source code.
2. Libraries including zlib (for gzip), libpthread (for multi-threading) and libm (for math). We tried to reduce our dependency as much as possible so no fancy libraries are needed.
3. GNUmake and Shell.

If you use Windows, you may consider to install [Mingw-w64](http://mingw-w64.org/doku.php) for that it provides all the required programs and libraries in one place.

Compiling the source code is simple. 
```sh
$ make -f Makefile.Linux (for Linux)
$ make -f Makefile.MacOS (for macOS)
$ make -f Makefile.Windows (for Windows)
```
The products will be moved to the "/bin" directory. Because our dependency is minimum, we do not use autoconf to generate the Makefiles. 

## Testing the installation
After installation, you may test the programs to see if it can work.

The Subread package incoperates many small test cases that cover most of its functions. No matter if Subread is installed from the source code or binary package, you can run "test_all.sh" in the "/test" directory. This assumes that you have a Shell program and a Python2 interpreter in PATH.

# Usage
The usages of the programs in this package can be found in the users-guide in the "/doc" directory.

# Citation
We have published papers on our Subread/Subjunc reada aligners and featureCounts read quantifiers:

1. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote ***Y Liao, GK Smyth, W Shi*** - Nucleic acids research, 2013

2. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features ***Y Liao, GK Smyth, W Shi*** - Bioinformatics, 2014
