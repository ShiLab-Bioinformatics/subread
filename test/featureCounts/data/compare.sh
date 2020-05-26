SAM_FILE=$1
ORA_FILE=$2
ANNO_FILE=$3
PARAMETERS=$4
TEST_NAME=$5
IS_FEATURE_LEVEL=$6


TMPF=data/DEL4-`date '+%s'`
printf "Testing %-60s [" "$TEST_NAME ... "
#echo nohup ../../bin/featureCounts $PARAMETERS -o $TMPF.FC -a $ANNO_FILE $SAM_FILE 
 nohup ../../bin/featureCounts $PARAMETERS -o $TMPF.FC -a $ANNO_FILE $SAM_FILE  >/dev/null 2>&1


if [ "$IS_FEATURE_LEVEL"  == "" ]
then
	cat $ORA_FILE $TMPF.FC |grep -v ^# |grep -iv Geneid | awk 'BEGIN{is_faild = 0; nl2=0; nl3=0} NF==2{ora[$1]=$2; is_faild++; nl2++} NF>3{if(ora[$1]==$7){is_faild -- } nl3++} END{if(is_faild || nl2!=nl3)printf("%c[31mFAILED%c[0m", 27,27);else printf("%c[32mPASS%c[0m", 27,27)}'
elif [ "$IS_FEATURE_LEVEL"  == "FL" ]
then
	cat $ORA_FILE $TMPF.FC |grep -v ^# |grep -iv Geneid | awk 'BEGIN{is_faild = 0; nl2=0; nl3=0} NF==5{ora[$1 $2 $3]=$5; is_faild++; nl2++} NF>6{if(ora[$1 $2 $3]==$7){is_faild -- } nl3++} END{if(is_faild||nl2!=nl3)printf("%c[31mFAILED%c[0m", 27,27);else printf("%c[32mPASS%c[0m", 27,27)}'
else
	cat $ORA_FILE $TMPF.FC |grep -v ^# |grep -iv Geneid | awk 'BEGIN{is_faild = 0; nl2=0; nl3=0} NF==2{ora[$1]=$2; is_faild++; nl2++} NF>3{if(ora[$1]==$7){is_faild -- } nl3++} END{if(is_faild || nl2!=nl3)printf("%c[31mFAILED%c[0m", 27,27);else printf("%c[32mPASS%c[0m", 27,27)}'
	lines_res=`cat $ORA_FILE.jcounts $TMPF.FC.jcounts |grep -v ^# |grep -v PrimaryGene|sort |uniq -c|awk '$1!=2' |wc -l`
	if [[ $lines_res -gt 0 ]]
	then
		echo |awk '{ printf(",%c[31mFAILED%c[0m", 27,27) }'
	else
		echo |awk '{ printf(",%c[32mPASS%c[0m", 27,27) }'
	fi
fi

echo "]"

rm -f $TMPF.FC*
