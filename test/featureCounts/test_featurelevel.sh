SH_CMD=bash


for test in intron_between across_genes across_intron
do
    for ams in data/$test*am 
    do
	$SH_CMD data/compare.sh $ams $ams.ora data/$test*gtf ' -p -f -s 2 ' $test FL
    done
done

