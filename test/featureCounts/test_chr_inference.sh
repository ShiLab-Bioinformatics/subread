SH_CMD=bash
mkdir -p result
echo
echo "================================================================================"
printf " FeatureCounts Chromosome Name Inference Tests\n http://subread.sourceforge.net/\n"
echo "================================================================================"
echo

# ================================================================================
# Testing incomplete chromosome names in the annotations and in the SAM file
$SH_CMD data/compare.sh data/test-chrname.sam data/test-chrname.ora data/test-chrname.SAF  "-F SAF -p " "automatic inference of chromosome names"

echo
