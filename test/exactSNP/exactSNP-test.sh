# This is a minimum runnable test case.
# The options are:
#    -g : reference sequence.
#    -o : output VCF.
#    -i : input read alignment; SAM by default, but we use a BAM here.
#    -b : the input read alignment file is a BAM file.
# More options are available.
mkdir -p result
echo ../../bin/exactSNP -g ../chr901.fa -o result/test-out.VCF -i data/test-in.BAM -b
../../bin/exactSNP -g ../chr901.fa -o result/test-out.VCF -i data/test-in.BAM -b
