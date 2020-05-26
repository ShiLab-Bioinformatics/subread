SUBREAD_HOME=../../bin/
PYTHON_EXEC=python

rm test-tmp.log
mkdir  -p result

$SUBREAD_HOME/subread-buildindex -B -F -o ../small1 -M100 ../chr901.fa
md5sum ../small1.00.b.array >> test-tmp.log
md5sum ../small1.00.b.tab >> test-tmp.log	

echo "*************************************************" >> test-tmp.log
echo "*** SINGLE-END READS NO ERROR              ******" >> test-tmp.log
echo "*************************************************" >> test-tmp.log
echo >>test-tmp.log

$SUBREAD_HOME/subread-align --SAMoutput -t0 -P6 -i ../small1 -r data/test-noerror-r1.fq -o result/test-tmp.sam -H -J
cat result/test-tmp.sam | $PYTHON_EXEC readname_ora_match.py >>test-tmp.log


echo "*************************************************" >> test-tmp.log
echo "*** SINGLE-END READS NO ERROR NO DUP       ******" >> test-tmp.log
echo "*************************************************" >> test-tmp.log
echo >>test-tmp.log

$SUBREAD_HOME/subread-align --SAMoutput  -t0 -P6 -i ../small1 -r data/test-noerror-r1.fq -o result/test-tmp.sam -H -J
cat result/test-tmp.sam | $PYTHON_EXEC readname_ora_match.py >>test-tmp.log



echo >>test-tmp.log
echo "*************************************************" >> test-tmp.log
echo "***  READS WITH NO ERROR                   ******" >> test-tmp.log
echo "*************************************************" >> test-tmp.log
echo >>test-tmp.log

$SUBREAD_HOME/subread-align  --SAMoutput  -t0 -P6 -i ../small1 -r data/test-noerror-r1.fq -R data/test-noerror-r2.fq -o result/test-tmp.sam -H -J
cat result/test-tmp.sam | $PYTHON_EXEC readname_ora_match.py >> test-tmp.log
echo >>test-tmp.log

echo "*************************************************" >> test-tmp.log
echo "***  READS NO ERROR, NO DUPLICATED REPORT  ******" >> test-tmp.log
echo "*************************************************" >> test-tmp.log
echo >>test-tmp.log

$SUBREAD_HOME/subread-align  --SAMoutput  -t0 -P6 -i ../small1 -r data/test-noerror-r1.fq -R data/test-noerror-r2.fq -o result/test-tmp.sam  -Q -J
cat result/test-tmp.sam | $PYTHON_EXEC readname_ora_match.py >>test-tmp.log


echo >>test-tmp.log

echo "*************************************************" >> test-tmp.log
echo "***  READS WITH ONLY SEQUENCING ERROR      ******" >> test-tmp.log
echo "*************************************************" >> test-tmp.log
echo >>test-tmp.log

$SUBREAD_HOME/subread-align  --SAMoutput  -t0 -P6 -i ../small1 -r data/test-error-r1.fq -R data/test-error-r2.fq -o result/test-tmp.sam -H -J
cat result/test-tmp.sam | $PYTHON_EXEC readname_ora_match.py >>test-tmp.log

echo >>test-tmp.log
echo "*************************************************" >> test-tmp.log
echo "***  READS WITH SEQUENCING ERROR AND MUTATION ***" >> test-tmp.log
echo "***  SUBREAD IS RUN WITH LONG INDEL DETECTION ***" >> test-tmp.log
echo "*************************************************" >> test-tmp.log
echo >>test-tmp.log

$SUBREAD_HOME/subread-align  --SAMoutput  -t0 -P6 -i ../small1 --gzFASTQinput -r data/test-err-mut-r1.fq.gz -R data/test-err-mut-r2.fq.gz -o result/test-tmp.sam -H -J --rg-id MyTestGroup --rg SM:sample1 --rg TP:1 --rg XX:YY 
cat result/test-tmp.sam | $PYTHON_EXEC readname_ora_match.py >>test-tmp.log

cat test-tmp.log
