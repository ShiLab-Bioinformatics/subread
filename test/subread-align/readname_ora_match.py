from sys import stdin

def main():
	unmatched = 0
	unmapped=  0
	matched = 0
	paired_match =0
	NN = 0

	line = 0
	old_read_name = None
	old_pos = 0

	while True:
		fl = stdin.readline()
		if not fl: break
		if fl[0]=='@':continue
		linfo = fl.split('\t')
		read_name = linfo[0]
		if len(linfo[2])>30:
			sam_chro = linfo[7]
			sam_offset = int(linfo[8])
		elif linfo[3]=='+' or linfo[3]=='-':
			sam_chro = linfo[1]
			sam_offset = int(linfo[2])
		else:
			sam_chro = linfo[2]
			sam_offset = int(linfo[3])
		line +=1

		if line % 100000 == 0:
			print "L=",line

		if fl.find("NNNNNNNN")>0:
			NN+=1
			continue

		if len(sam_chro)<2:
			unmapped +=1
			continue

		if read_name.find ('/')>0:
			name_info = read_name.split('/')
			pair_number = int(name_info[1])

		if read_name.find ('/')>0:
			name_info = read_name.rsplit('_',5)
		else:
			name_info = read_name.rsplit('_',5)

		if len(name_info)<3: continue

		if name_info[2].count(':')>0:
			ora_chro = sam_chro 
			ora_pos1 = int(name_info[1]) 
			ora_pos2 = int(name_info[1])
		else:
			ora_chro = name_info[0]
			ora_pos1 = int(name_info[1])
			ora_pos2 = int(name_info[2])

		t_name = ora_chro+name_info[1]+name_info[2]
		if t_name != old_read_name:
			old_pos = sam_offset

		if (not sam_chro.endswith( ora_chro)) or not (abs(ora_pos1 - sam_offset) < 1200 or abs(ora_pos2 - sam_offset)<1200):
			unmatched +=1
#			if sam_chro != ora_chro:
#			print read_name,"\t", sam_chro, ora_chro, ora_pos1 , sam_offset, ora_pos2, fl.strip()
		else:
			matched +=1
			if  t_name  == old_read_name and  (abs(old_pos - ora_pos1) < 1200):
				paired_match +=2
#			else:
#				print "Unpired_match: old=", old_pos," new1", sam_chro, sam_offset, "ora1", ora_chro, ora_pos1 

		old_read_name = t_name

	print "unmatched=",unmatched, ";   matched=",matched, ";     unmapped=", unmapped, ";     reads=", line, "    ;NN=", NN
	print "accuracy=",(matched*1.)/(matched+unmatched)," ;    sensitivity=",(matched+unmatched)*1./(line-NN)
	print "paired_match=", paired_match,"  ;     paired=", (paired_match*1./(matched))

main()
