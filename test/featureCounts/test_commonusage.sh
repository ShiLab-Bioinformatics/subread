echo
echo "================================================================================"
printf " FeatureCounts Common Scenario Tests\n http://subread.sourceforge.net/\n"
echo "================================================================================"
echo

infiles="/usr/local/work/liao/Rsubread/testsuit/SEQC2011_for_3tests_wild-coorsorted-Picard-SamtAgain.bam /usr/local/work/liao/Rsubread/testsuit/SEQC2011-A-FCpaper-SAM1.4.junc"

for inf in $infiles
do
  intag="complex"
  if [[ $inf =~ FCpaper ]]
  then
    intag="simple"
  fi
  for level in GENE EXON
  do 
    op_level=
    if [[ $level == "EXON" ]]
    then
       op_level=" -f "
    fi
    for moverlap in YES NO
    do 
      op_moverlap=
      if [[ $moverlap == "YES" ]]
      then
        op_moverlap=" -O "         
      fi
      for mmapping in YES NO
      do 
        op_mmapping=
        if [[ $mmapping == "YES" ]]
        then
          op_mmapping=" -M "         
        fi
        for fraction in YES NO
        do 
          op_fraction=
          if [[ $fraction == "YES" ]]
          then
            if [[ $mmapping == "NO" ]]
            then
		continue
	    fi
            op_fraction=" --fraction "
          fi
	  ora=data/commonusage-$intag-$level-$moverlap-$mmapping-$fraction.FC
	  echo $ora
          ../../bin/featureCounts -a /usr/local/work/liao/Rsubread/testsuit/ensembl-for_3tests-shuf.GTF -o data/del4.FC $op_level $op_moverlap $op_mmapping $op_fraction  -T7 -p $inf  &>/dev/null
          cat $ora|grep -v ^#|grep -vi Geneid |md5sum
	  cat data/del4.FC |grep -v ^#|grep -vi Geneid |md5sum

          cat $ora.summary|grep -v ^#|grep -vi Status |md5sum
	  cat data/del4.FC.summary |grep -v ^#|grep -v Status |md5sum
	  rm -f data/del4*
        done
      done
    done
  done
done

