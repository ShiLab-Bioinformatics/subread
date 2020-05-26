mkdir -p result
echo
echo "================================================================================"
printf " FeatureCounts Corner Case Tests\n http://subread.sourceforge.net/\n"
echo "================================================================================"
echo



SH_CMD=bash
$SH_CMD data/compare.sh data/corner-INDEL.sam data/corner-INDEL.ora data/test-minimum.GTF "-p" "indel reads"
$SH_CMD data/compare.sh data/corner-JUNC.sam data/corner-JUNC.ora data/test-minimum.GTF "-p" "junction reads"
$SH_CMD data/compare.sh data/corner-ONEEND.sam data/corner-ONEEND.ora data/test-minimum.GTF "-p" "paired-end reads (fragment counting)"
$SH_CMD data/compare.sh data/corner-ONEEND.sam data/corner-ONEEND-BOTH.ora data/test-minimum.GTF "-p -B " "paired-end reads (fragment counting, both ends mapped)"
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum-O.ora data/test-minimum.GTF "-p -O " "multi-overlapping reads"
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum-FL.ora data/test-minimum.GTF "-p -f " "feature-level summarization" FL
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum.ora data/test-minimum.GTF "-p " "gene-level summarization"
$SH_CMD data/compare.sh data/corner-NH.sam data/corner-NH.ora data/test-minimum.GTF "-p" "multi-mapping reads"
$SH_CMD data/compare.sh data/corner-NH.sam data/corner-NH-PM.ora data/test-minimum.GTF "-p --primary -M " "multi-mapping reads (primary only)"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-BothEnds.ora data/test-minimum.SAF "-p -F SAF -B " "both ends mapped"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-Chimeric.ora data/test-minimum.SAF "-p -F SAF -C " "disallowing chimeric fragments"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-MultiMapping.ora data/test-minimum.SAF "-p -F SAF -M " "Allowing multi-mapped reads"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-DoNotSort.ora data/test-minimum.SAF " -p -F SAF --donotsort " "not sorting input file"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-MinOverlap.ora data/test-minimum.SAF " --minOverlap 125 -p -F SAF " "minimum overlapping length"
$SH_CMD data/compare.sh data/test-fracOverlap.sam data/test-fracOverlap.ora data/corner-fractions.SAF " --fracOverlap 0.62 -O -p -F SAF " "minimum overlapping fraction"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-LargestOverlap.ora data/test-minimum.SAF "-p -F SAF --largestOverlap" "Largest Overlapping"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-PEdist.ora data/test-minimum.SAF " -p -F SAF -B -C -P -d 130 -D 770 " "paired-end distance"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-Read2Pos5.ora data/test-minimum.SAF " -p -F SAF --read2pos 5 " "Read to position (5' end)"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-Read2Pos3.ora data/test-minimum.SAF " -p -F SAF --read2pos 3 " "Read to position (3' end)"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-Extend3.ora data/test-minimum.SAF " -p -F SAF --readExtension3 1000 " "Read extension to the 3' end"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-Extend5.ora data/test-minimum.SAF " -p -F SAF --readExtension5 1000 " "Read extension to the 5' end"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-MaxOPs.ora data/test-minimum.SAF " -p -F SAF --maxMOp 2 " "Low maxOPs value"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-MinMAPQ.ora data/test-minimum.SAF " -p -F SAF -Q 58" "minimum mapping quality"
$SH_CMD data/compare.sh data/test-dup.sam data/corner-IgnoreDup.ora data/test-minimum.SAF "-p -F SAF --ignoreDup " "Ignoring duplicated reads"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-Fraction.ora data/test-minimum.SAF "-p -F SAF --fraction -M " "Fraction counting"
$SH_CMD data/compare.sh data/corner-fractions.sam data/corner-fractions.ora data/corner-fractions.SAF "  -O -M -F SAF --fraction  " "Advanced fractions"
$SH_CMD data/compare.sh data/test-junc.sam data/corner-Jcounts.ora data/test-minimum.SAF "-p -F SAF -J " "Junction counting" JC

if test -f /usr/local/work/work/liao/subread/chromosomes/all_34_alt.fa
then
	$SH_CMD data/compare.sh data/test-junc.sam data/corner-Jcounts-FA.ora data/test-minimum.SAF "-p -F SAF -J --genome /usr/local/work/work/liao/subread/chromosomes/all_34_alt.fa " "Junction counting (with genome) " JC
else
	echo "Skipping Junction counting (with genome)."
fi


# by default it is in GTF.
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum.ora data/test-minimum.GTF "-p " "GTF format annotations"
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum.ora data/test-minimum.SAF "-p -F SAF " "SAF format annotations"

# by default it is in SAM.
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum.ora data/test-minimum.GTF "-p " "SAM format input"
$SH_CMD data/compare.sh data/test-minimum.bam data/test-minimum.ora data/test-minimum.GTF "-p " "BAM format input" 

# by default it is non-strand specific.
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum.ora data/test-minimum.GTF "-p -s 0 " "unstranded read summarization"
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum-STR.ora data/test-minimum.GTF "-p -s 1 " "stranded read summarization"
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum-UNSTR.ora data/test-minimum.GTF "-p -s 2 " "reversely stranded read summarization"

# test 5' and 3' end extension
$SH_CMD data/compare.sh data/test-chrname.sam data/test-minimum-dup.ora data/test-minimum.GTF " -p --ignoreDup " "Ignoring duplicate fragments" 
$SH_CMD data/compare.sh data/corner-JUNC.sam data/corner-JUNC-ONLY.ora data/test-minimum.GTF " --splitOnly -O -f " "Junction reads only" FL
$SH_CMD data/compare.sh data/corner-JUNC.sam data/corner-EXON-ONLY.ora data/test-minimum.GTF " --nonSplitOnly -p " "Exonic reads only" 

echo
