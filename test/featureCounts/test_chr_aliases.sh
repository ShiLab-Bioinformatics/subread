SH_CMD=bash
mkdir -p result
echo
echo "================================================================================"
printf " FeatureCounts Chromosome Name Aliases Tests\n http://subread.sourceforge.net/\n"
echo "================================================================================"
echo

# ================================================================================
# Testing alias file to convert chromosome names in the annotation file
# The alias file has each line defining an alias: Chro_Name_in_Annotation,Chro_Name_in_SAM
$SH_CMD data/compare.sh data/test-chralias.sam data/test-chralias.ora data/test-chralias.SAF  "-F SAF -p -A data/test-chralias.txt " "chromosome aliases"
echo
