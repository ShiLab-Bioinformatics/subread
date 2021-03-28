mkdir -p result
cat <<EOF > /dev/null
 The options in the command lines below are:
  -a : annotation file (GTF by default)
  -F : annotation format (GTF by default)
  -A : chromosome alias file
  -i : input file for reads (SAM by default)
  -b : input file is in the BAM format
  -o : output file
  -p : paired-end assignment
  -f : feature-level (exon level) assignment
  -O : allowing a read to overlap with multi features 
  -S : resorting the input SAM or BAM file
  -T : threads for assignment
  -d and -D : minimum and maximum allowed template lengths
  -B : both ends must be mapped in a paired-end fragment
  -C : no chimeric fragments are allowed
  -M : multi-mapping reads reported by the aligner are allowed

 More options are available. Reference to the user guide for the full option list.

EOF


echo
echo "================================================================================"
printf " FeatureCounts Basic Test\n http://subread.sourceforge.net/\n"
echo "================================================================================"
echo



rm -f data/test-minimum.log
# ================================================================================
# The minimum runnable test
../../bin/featureCounts -p -a data/test-minimum.GTF -o result/test-minimum.FC data/test-minimum.sam 

echo "================================================================================"
echo "Basic Test finished."
echo "The results are in result/test-minimum.FC"
echo "================================================================================"
echo
echo
echo
