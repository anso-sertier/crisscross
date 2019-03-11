#!/home/anneso/myroot/bin/python
# -*- encoding: utf8 -*-

## CREATED 22/08/2017 by anso
## MODIFIED 24/04/2018
## - take ref genome b38 into account (chromosome name are chr1, chr2, etc ...)

import argparse
from itertools import dropwhile
import re
from pyfaidx import Fasta


def findTag(tag,lst):
	for i in lst:
		if i.startswith(tag):
			return i.split(":")[2]
	return ""

def getStrand(flag):
	if flag == 0:
		return "F"
	elif (flag & 16)== 16:
		return "R"
	elif (flag & 4) == 4:
		return ""
	else:
		print flag

def getEndAli(start,cigar):
	nbs = re.findall("[0-9]+",cigar)
	ope = re.findall("[SMID]+",cigar)
	endali = start
	for i in range(len(ope)):
		o = ope[i]
		if o in ["M","D"]:
			endali += int(nbs[i])
	return endali

def getRealignBkptPos(scstrand, mapstrand,alistart,aliend):
	if scstrand == "5" and mapstrand == "F":
		return aliend,"5'"
	elif scstrand == "5" and mapstrand == "R":
		return alistart-1,"3'"
	elif scstrand == "3" and mapstrand == "F":
		return alistart-1,"3'"
	else:
		return aliend,"5'"

def readSamFile(fic,CHROMOSOMES):
	alignments   = {key :{} for key in CHROMOSOMES}
	sclst        = set()
	nb_ali = nb_chr = nb_unaligned = nb_notkept = nb_R = nb_N = nb_U = 0
	for line in dropwhile(lambda line: line.startswith('@'), fic):
		nb_ali   += 1
		flag_ali  = ""
		alistart,aliend,bkptpos,sideali,xt = [-1,-1,-1,"",""]
		l         = line.strip().split("\t")
		chrom_ali = l[2]
		strand    = getStrand(int(l[1]))
		chrom,pos,side,numCons,nb,lg = l[0].split(":")
		scname = "%s:%s:%s':%s" % (chrom,pos,side,numCons)
		if strand == "":
			nb_unaligned += 1
			flag_ali = "unaligned"
		elif chrom_ali not in CHROMOSOMES:
			nb_chr += 1
			flag_ali = "badchrom"
		else:
			cigar    = l[5]
			alistart = int(l[3])
			aliend   = getEndAli(alistart-1,cigar)
			xt  = findTag("XT",l[11:])      # type (U/R/N/M)
			x0  = int(findTag("X0",l[11:])) # number of best hits
			xm  = int(findTag("XM",l[11:])) # number of mismatch
			bkptpos,sideali = getRealignBkptPos(side, strand, alistart, aliend ) ## WARNING side without "'"
			if xt == "U" or (xt == "R" and xm < 2 and cigar == lg+"M" and x0 < 5 ):
				flag_ali = "goodali"
				sclst.add(scname)
			else:
				flag_ali = "badali"
				nb_notkept += 1
			if xt == "R":
				nb_R += 1
			elif xt == "N":
				nb_N += 1
			elif xt == "U":
				nb_U += 1
		try:
			alignments[chrom][int(pos)][scname] =   (flag_ali,chrom_ali,bkptpos,sideali)
		except KeyError:
			alignments[chrom][int(pos)] = {scname : (flag_ali,chrom_ali,bkptpos,sideali)}
	fic.close()
	print "ALIGNMENT SUMMARY\t: %s realigned" % nb_ali
	print "\tAlignment outside chr : ", nb_chr
	print "\tConsensus not found   : ", nb_unaligned
	print "\tAlignment not kept    : ", nb_notkept
	print "\tAlignment kept        : ", len(sclst)
	print "Alignment XT = R:%s ; N:%s ; U:%s" % (nb_R,nb_N,nb_U)
	return alignments

def getScposDefined(filename,CHROMOSOMES):
	scpos = {key :{} for key in CHROMOSOMES}  # chrom : pos : name : [nbr, maxsc, seq]
	fic = open(filename,"r")
	fic.next()
	for line in fic:
		chrom,pos,side,numCons,nbr,maxsc,seq = line.strip().split()[0:7]
		scname = "%s:%s:%s:%s" % (chrom,pos,side,numCons)
		try:
			scpos[chrom][int(pos)][scname] = [nbr,maxsc,seq]
		except KeyError:
			scpos[chrom][int(pos)] = {scname : [nbr,maxsc,seq]}
	return scpos

def searchSimpleRepeats(seq):
	for x in range(2, 7):
		motif = seq[0:x]
		nb = 1
		y = x
		while True:
			if motif == seq[y:(y+x)] :
				nb += 1
				y += x
			else:
				break
		if nb > 1:
			return motif,nb
	return "",0

def analyseRefSequence(chrom,pos1,pos2,fasta):
	dist = abs(pos1-pos2)
	if pos1 < pos2 :
		seq = fasta[chrom][(pos1-1):pos2].seq
	else:
		seq = fasta[chrom][(pos2-1):pos1].seq
	motif,nb = searchSimpleRepeats(seq)
	if nb == 0 or len(motif)*nb < 0.5*dist: # if total repeat length is less than 50% of deleted or inserted sequence
		return "",0
	else:
		return motif,nb

def getDistItself(side1,pos1,pos2):
	if side1 == "5'":
		return pos2 - pos1
	else:
		return pos1 - pos2

def main():
	parser = argparse.ArgumentParser(description='resolve sc position pairing')
	parser.add_argument("-i",action="store", type=file,                   required=True, dest="sam",    help="sam realignment",            metavar="inputfilename")
	parser.add_argument("-s",action="store", type=str,                    required=True, dest="cons",   help="summary with cons",          metavar="inputfilename")	
	parser.add_argument("-r",action="store", type=str,                    required=True, dest="ref",    help="fasta reference genome",     metavar="inputfilename")	
	parser.add_argument("-o",action="store", type=argparse.FileType('w'), required=True, dest="out",    help="out bkpt reali summary",     metavar="outfilename")
	parser.add_argument("-m",action="store", type=int, default = 100,     required=False,dest="maxdist",help="max distance between sc pos",metavar="int")
	args = parser.parse_args()

	fasta       = Fasta(args.ref)
	CHROMOSOMES = fasta.keys()[0:24]
	
	alignments = readSamFile(args.sam,CHROMOSOMES)
	all_scpos  = getScposDefined(args.cons,CHROMOSOMES)
	# alignements = {chrom : pos : name : [flag_ali, chromali, bkptpos, side]}
	# all_scpos   = {chrom : pos : name : [nr, lg, seq]}

	to_out = {key :{} for key in CHROMOSOMES} ## {chrom : pos : side : numCons : [nr,lg,flag,info,seq]}
	nb_trash = 0
	nb_all = 0 ; nb_intra = 0 ; nb_inter = 0 ; nb_ins = 0 ; nb_del = 0 ; nb_ident = 0
	args.out.write("scpos\tflag\tchr_ali\tpos_ali\tside_ali\tdist\tmotif\tnbMotif\n")
	outline = "%s\t"*7 + "%s\n"
	for chrom in CHROMOSOMES:
		allpos = sorted(all_scpos[chrom].keys())
		for pos in allpos:
			for scname,infos in all_scpos[chrom][pos].iteritems():
				nb_all += 1
				side,numcons = scname.split(":")[2:4]
				newflag = ""
				try:
					flag_ali, chrom_ali, pos_ali, side_ali = alignments[chrom][pos][scname]
					pos_ali = int(pos_ali)
					if flag_ali in ["unaligned","badchrom","badali"]:
						newflag = "reali_"+flag_ali
					elif flag_ali == "goodali":
						if chrom == chrom_ali:
							dist = abs(pos-pos_ali)
							if dist > args.maxdist:
								nb_intra += 1
								newflag = "reali_intra"
								args.out.write(outline % (scname,newflag,chrom_ali,pos_ali,side_ali,dist,"",0))
							else:
								if side == side_ali:
									if   (side == "3'" and pos == pos_ali -1) or (side == "5'" and pos == pos_ali):
										nb_trash += 1
										continue
									elif (side == "3'" and pos > pos_ali )    or (side == "5'" and pos < pos_ali):
										nb_ins += 1
										newflag = "reali_ins"
										motif,nb = analyseRefSequence(chrom,pos,pos_ali,fasta)
										args.out.write(outline % (scname,newflag,chrom_ali,pos_ali,side_ali,dist,motif,nb))
									elif (side == "3'" and pos <= pos_ali )   or (side == "5'" and pos > pos_ali):
										nb_del += 1
										newflag = "reali_del"
										motif,nb = analyseRefSequence(chrom,pos,pos_ali,fasta)
										args.out.write(outline % (scname,newflag,chrom_ali,pos_ali,side_ali,dist,motif,nb))
									else:
										print "Error goodali: ", chrom, pos, scname, bkpt_reali
								else: #itself
									nb_ident += 1
									dist = getDistItself(side,pos,pos_ali) # negative if overlaps
									newflag = "reali_itself"
									args.out.write(outline % (scname,newflag,chrom_ali,pos_ali,side_ali,dist,"",0))
						else:
							nb_inter += 1
							newflag = "reali_inter"
							args.out.write(outline % (scname,newflag,chrom_ali,pos_ali,side_ali,-1,"",0))
					else:
						print "Error : flag for alignement is missing : %s:%s %s %s" % (chrom,pos,side,numcons)
				except KeyError:
					newflag = "not_realigned"
				## ADD positions to results
				if to_out[chrom].has_key(pos):
					try:
						to_out[chrom][pos][side][numcons]   = [infos[0],infos[1],newflag,infos[2]]
					except KeyError:
						to_out[chrom][pos][side] = {numcons : [infos[0],infos[1],newflag,infos[2]]}
				else:
					to_out[chrom][pos] = {side :   {numcons : [infos[0],infos[1],newflag,infos[2]]}}
	args.out.close()
	
	out = open(args.cons,"w")
	out.write("chr\tpos\tside\tnum\tnr\tlg\tflag\tseq\n")
	outline = "%s\t"*4 + "%s\n"
	for chrom in CHROMOSOMES:
		try:
			for pos in sorted(to_out[chrom].keys()):
				for side in sorted(to_out[chrom][pos]):
					for numCons in sorted(to_out[chrom][pos][side]):
						dd = "\t".join(to_out[chrom][pos][side][numCons])
						out.write(outline % (chrom,pos,side,numCons,dd))
		except KeyError:
			pass
	out.close()
	
	print "ANALYSE ALL SC POS : %s" % nb_all
	print "\tEqual to reference : %s" % nb_trash
	print "\tNb intra reali bkpt added : %s " % nb_intra
	print "\tNb ident reali bkpt added : %s " % nb_ident
	print "\tNb ins   reali bkpt added : %s " % nb_ins
	print "\tNb del   reali bkpt added : %s " % nb_del
	print "\tNb inter reali bkpt added : %s " % nb_inter


if __name__=="__main__":
    main()




