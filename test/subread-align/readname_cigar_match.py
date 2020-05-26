from sys import stdin

def cigar_method(cigar, method):
	ret = 0;
	tmp_ret = 0;
	for cc in cigar:
		if cc.isdigit():
			tmp_ret = tmp_ret*10+int(cc)
		else:
			if cc==method:
				ret += tmp_ret
			tmp_ret=0
	return ret

all_r = 0
unmapped = 0
error_pos = 0
error_cig = 0
correct_indel = 0
while True:
	fl = stdin.readline()
	if not fl: break
	fl=fl.strip()
	fls=fl.split('\t')
	if len(fls)<5: continue

	read_name = fls[0]
	chro = fls[2]
	chro_pos = int(fls[3])
	cigar = fls[5]

	# R_chr12_110788490_49M1I50M

	oracle_info = read_name.split('_')
	oracle_cigar = oracle_info[-1]
	oracle_pos = int(oracle_info[-2])
	oracle_chro = '_'.join(oracle_info[1:-2])
	oracle_is_reverse = oracle_info[0]=='R'

	is_correct = True
	is_cigar_problem = False

	all_r +=1

	if len(chro)<2:
		unmapped+=1
		continue

	if chro != oracle_chro : is_correct = False

	if is_correct:
		if abs(oracle_pos  - chro_pos) > 5: is_correct = False

	if is_correct:
		if cigar_method(oracle_cigar,'I') != cigar_method(cigar,'I'):
			is_correct = False
			is_cigar_problem = True
		if cigar_method(oracle_cigar,'D') != cigar_method(cigar,'D'):
			is_correct = False
			is_cigar_problem = True
		if cigar_method(oracle_cigar,'N') != cigar_method(cigar,'N'):
			is_correct = False
			is_cigar_problem = True

	if cigar_method(oracle_cigar,'I')>0 or cigar_method(oracle_cigar,'D')>0:
		if is_correct:
			correct_indel+=1

	if not is_correct:
		if is_cigar_problem:
			print "%s\t%s"%(cigar, oracle_cigar)
			error_cig +=1
		else:
			error_pos +=1
			print "%s,%d\t%s,%d\t\t%s"%(chro, chro_pos, oracle_chro, oracle_pos, read_name)

print "RES: %d reads\t%d unmapped\t%d error_pos\t%d error_cigar\t%d correct_indel"%(all_r,unmapped, error_pos, error_cig, correct_indel)
