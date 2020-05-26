#!bash
echo |awk '{printf("%c[2J%c[0;0H%c[30;47m", 27,27,27)}'
echo
echo
echo
echo
echo "  **************************************************   "
echo "  **************************************************   "
echo "  ***                                            ***   "
echo "  *** This script will test the major functions  ***   "
echo "  *** in our package, including the index build- ***   "
echo "  *** er, subread-align, subjunc, featureCounts  ***   "
echo "  *** and exactSNP.                              ***   "
echo "  ***                                            ***   "
echo "  *** Test will start in   seconds.              ***   "
echo "  ***                                            ***   "
echo "  **************************************************   "
echo "  **************************************************   "
echo

echo |awk -v secs=9 '{printf("%c[13;26H%c[33;44m%s",27,27,secs)}'
for secs in {0..8}
do
	sleep 1
	sec0=` echo 8-$secs | bc `
	echo |awk -v secs=$sec0 '{printf("%c[13;26H%s",27,secs)}'
done

echo |awk '{printf("%c[0m%c[2J%c[0,0H",27,27, 27)}'
#test subread-align
cd subread-align
sh subread-align-test.sh

#test subjunc
cd ../subjunc
sh subjunc-test.sh

#test featureCounts
cd ../featureCounts
sh featureCounts-test.sh

#test exactSNP
cd ../exactSNP
sh exactSNP-test.sh

echo |awk '{printf("%c[30;47m", 27)}'
echo
echo
echo "  **************************************************   "
echo "  **************************************************   "
echo "  ***                                            ***   "
echo "  *** Test finished.                             ***   "
echo "  ***                                            ***   "
echo "  *** Should there be any error, please visit    ***   "
echo | awk '{printf("  *** %c[34mhttp://subread.sourceforge.net/%c[0;30;47m for more   ***   \n", 27,27);}'
echo "  *** information.                               ***   "
echo "  ***                                            ***   "
echo "  **************************************************   "
echo "  **************************************************   "
echo
echo |awk '{printf("%c[0m", 27)}'
echo
echo
