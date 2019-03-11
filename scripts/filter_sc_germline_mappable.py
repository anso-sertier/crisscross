#!/home/anneso/myroot/bin/python
# -*- encoding: utf8 -*-

## CREATED 24/08/2015 by anso
## MODIFIED 9/09/2016
## 		- SPLIT into : germline / specific : normal-somatic
## 		- 2 positions found in tumor and normal but with different consensus : not germline : count how many times it appends
##		- modify ReadCons : fasta file of consensus format has changed
##		- modify Write* : idem format has changed
## MODIFIED 23/08/2017
##		- add mappable region overlap
## MODIFIED 24/04/2018
##		- change CHROMOSOMES definition to take b38 ref into account

import argparse
import re
import multiprocessing as mp
import operator


def ReadCons(filename):
	fic = open(filename,"r")
	fic.next()
	consensus = {}
	for line in fic:
		l = line.strip().split("\t")
		name = ":".join(l[0:3])
		try:
			consensus[name][l[3]] = (l[7],l[4])
		except KeyError:
			consensus[name] = {l[3] : (l[7],l[4])}
	fic.close()
	return consensus

def ReadMappableOverlap(fic):
	scpos = {}
	for line in fic:
		l = line.strip().split("\t")
		mstart = int(l[4])
		mstop  = int(l[5])
		pos    = int(l[2])
		nbp    = int(0.05*(mstop-mstart))
		if pos <  mstart + nbp or pos > mstop  - nbp:
			flag = "border"
		else:
			flag = "unmap"
		scpos[l[0] + ":" + l[2]] = flag
	return scpos

def PairwiseScore(seq1, seq2):
	common = 0
	for s1, s2 in zip(seq1, seq2):
		if s1 == s2:
			common += 1
	return common*2

def ComputeAlnCons(consensus):
	nameSC,constu,consno = consensus
	side       = nameSC.split(":")[2]
	germline   = {nameSC : {}}
	lst = []
	for numTU in constu:
		seqTU  = constu[numTU][0]
		max_scores = {}
		for n in consno:
			max_scores[n] = 2*len(min([consno[n][0], seqTU],key=len))
		if (side == "5'"): ## align from right to left
			seq = seqTU[::-1]
			scores = [(PairwiseScore(seq, consno[i][0][::-1])*100.0/max_scores[i], len(consno[i][0]) ,i,numTU) for i in consno]
		elif (side == "3'"): ## align from left to right
			scores = [(PairwiseScore(seqTU, consno[i][0])*100.0/max_scores[i], len(consno[i][0]), i,numTU) for i in consno]
		else:
			print "ERROR : bad side : ", side
		kept = [i for i in scores if i[0] >= 90]
		lst.extend(kept)
	lst.sort(key=operator.itemgetter(0,1))
	already_kept = []
	while lst != []:
		score,lg,numNO,numTU = lst.pop()
		if (not numNO in already_kept) and (not numTU in already_kept):
			germline[nameSC][numTU] = numNO
			already_kept.extend([numNO,numTU])
	return germline

def WriteSpecific(filename,outtxt,germline,unmap):
	fic_sc  = open(filename,"r")
	header  = fic_sc.readline()
	outtxt.write(header.strip("\n") + "\tmap\n")
	germline_sc  = {}
	for line in fic_sc:
		l       = line.strip("\n").split("\t")
		namepos = ":".join(l[0:2])
		namesc  = ":".join(l[0:3])
		numcons = l[3]
		## get mappable flag
		try:
			flag = unmap[namepos]
		except KeyError:
			flag = "mappable"
		## if specific write, else add 
		try:
			if numcons in germline[namesc]:
				try:
					germline_sc[namesc][numcons] = l[4:8]  # nr,lg,flag,seq
				except KeyError:
					germline_sc[namesc] = {numcons : l[4:8]}
			else:
				outtxt.write(line.strip("\n") + "\t"+ flag +"\n")
		except KeyError:
			outtxt.write(line.strip("\n") + "\t"+ flag +"\n")
	fic_sc.close()
	outtxt.close()
	return germline_sc

def WriteGermline(germline,consensusTU,consensusNO,outtxt,unmap):
	outtxt.write("chrom\tpos\tside\tnum\tnrTot\tmaxLg\tflag\tseq\tmap\tconsTU\tconsNO\n")
	outline = "%s\t"*10 + "%s\n"
	if germline.keys()[0].startswith("chr"):
		CHROMOSOMES = ["chr"+str(i) for i in range(1,23)]
		CHROMOSOMES.extend(["chrX","chrY"])
	else:
		CHROMOSOMES = [str(i) for i in range(1,23)]
		CHROMOSOMES.extend(["X","Y"])

	for chrom in [i for i in CHROMOSOMES if i in germline.keys()]:
		for pos in sorted(germline[chrom].keys()):
			namepos = chrom + ":" + str(pos)
			try:
				flagmapp = unmap[namepos]
			except KeyError:
				flagmapp = "mappable"
			for side in sorted(germline[chrom][pos].keys()):
				namesc = "%s:%s:%s" % (chrom,pos,side)
				for numTU in germline[chrom][pos][side]:
					numNO  = germline[chrom][pos][side][numTU]
					consTU = consensusTU[namesc][numTU] # nr,lg,flag,seq
					try:
						consNO = consensusNO[namesc][numNO]
					except KeyError:
						print namesc,numNO,consensusNO[namesc]
						return
					nrTot = int(consTU[0]) + int(consNO[0])
					maxlg = max(int(consTU[1]), int(consNO[1]))
					flagali = consNO[2]
					if int(consNO[1]) == maxlg:
						seq = consNO[3]
					else:
						seq = consTU[3]
					txtTU = "%s,%s,%s" % (numTU,consTU[0],consTU[2])
					txtNO = "%s,%s,%s" % (numNO,consNO[0],consNO[2])
					outtxt.write(outline % (chrom,pos,side,numNO,nrTot,maxlg,flagali,seq,flagmapp,txtTU,txtNO))
	outtxt.close()



def main():
	parser = argparse.ArgumentParser(description='First filtering of germline abnormal pairs')
	parser.add_argument("-i",action="store", type=str,  required=True, dest="tusc",    help="tab file of tumor sc pos",metavar="file")
	parser.add_argument("-j",action="store", type=str,  required=True, dest="nosc",    help="tab file of normal sc pos",metavar="dir")
	parser.add_argument("-k",action="store", type=file,  required=True, dest="maptu",   help="intersect tu/mappable",metavar="dir")
	parser.add_argument("-l",action="store", type=file,  required=True, dest="mapno",   help="intersect no/mappable",metavar="dir")
	parser.add_argument("-t",action="store", type=int, 	required=True, dest="nbtread", help="nb cpu to use",metavar="number")
	parser.add_argument("-a",action="store", type=argparse.FileType('w'), required=True, dest="outtu",  help="output directory",metavar="dir")
	parser.add_argument("-b",action="store", type=argparse.FileType('w'), required=True, dest="outno",  help="output directory",metavar="dir")
	parser.add_argument("-c",action="store", type=argparse.FileType('w'), required=True, dest="outgerm",  help="output directory",metavar="dir")
	args = parser.parse_args()

	## COMPARE CONSENSUS SEQUENCES OF NORMAL AND TUMOR COMMON SoftClip
	pool    = mp.Pool(processes=(args.nbtread))
	results = pool.map(ReadCons, [args.tusc,args.nosc])
	pool.close()
	pool.join()
	tuCons = results[0]   # chrom:pos:side : num : (seq,nr)
	noCons = results[1]
	nameTU = set(tuCons.keys())
	nameNO = set(noCons.keys())
	inputs = [(name, tuCons[name], noCons[name]) for name in nameTU.intersection(nameNO)]   # [ ( chrom:pos:side, {consTU} , {consNO} )]
	pool    = mp.Pool(processes=args.nbtread)
	results = pool.map(ComputeAlnCons, inputs)  # [{chrom:pos:side : { numConsTu : numConsNo }}]
	pool.close()
	pool.join()
	
	germlineTU = {}		 # chrom:pos:side : [numTU]
	germlineNO = {}		 # chrom:pos:side : [numNO]
	germline_sorted = {} # chrom : pos : side : numTU : [numNO]
	for i in results:
		tmp = i.copy()
		if len(tmp) == 0:
			continue
		namesc          = tmp.keys()[0]
		chrom,pos,side  = namesc.split(":")
		if germline_sorted.has_key(chrom):
			try:
				germline_sorted[chrom][pos][side] = {}
			except KeyError:
				germline_sorted[chrom][pos] = {side : {}}
		else:
			germline_sorted[chrom] = {pos : {side : {}}}
		for numTU in tmp[namesc]:
			numNO = tmp[namesc][numTU]
			germline_sorted[chrom][pos][side][numTU] = numNO
			try:
				germlineTU[namesc].append(numTU)
			except KeyError:
				germlineTU[namesc] = [numTU]
			try:
				germlineNO[namesc].append(numNO)
			except KeyError:
				germlineNO[namesc] = [numNO]
	
	## READ MAPPABILITY
	tumap = ReadMappableOverlap(args.maptu) # chrom:pos : tag
	nomap = ReadMappableOverlap(args.mapno)
	## WRITE OUTPUTS
	germ_consTU = WriteSpecific(args.tusc, args.outtu, germlineTU, tumap)
	germ_consNO = WriteSpecific(args.nosc, args.outno, germlineNO, nomap)
	WriteGermline(germline_sorted,germ_consTU,germ_consNO,args.outgerm,tumap)

	print "done."
	return 0

if __name__=="__main__":
    main()
