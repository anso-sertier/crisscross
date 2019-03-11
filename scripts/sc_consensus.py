#!/usr/bin/env python
# -*- encoding: utf8 -*-

## CREATED 3/08/2016 by anso

## MODIFICATIONS:
##		- search with custom distance no pairwise function
## MODIFICATIONS 4/01/2017:
##		- remove fasta output
##		- summary file format changes : 1 lane per consensus + nt seq
## MODIFICATIONS 24/08/2017:
##		- add partial reference alignment informations
##		- if sc pos already exist where seq is different from reference then add this read to already existing position
##		- otherwise create two positions : this one with a tag and a new one where ? or delete


import argparse
import re
import multiprocessing as mp
import itertools

def ReadFasta(fic):
	sequences = []
	line = fic.readline()
	if not line.startswith(">"):
		print "Error in fasta file format : first line should start with >"
	l = line.strip().split()
	scpos = l[1]
	side = l[2]
	seqs = []
	for line in fic:
		if line.startswith(">"):
			## treat previous
			sequences.append((scpos,side,seqs))
			## initialise new
			l = line.strip().split()
			scpos = l[1]
			side = l[2]
			seqs = []
		else:
			seqs.append(line.strip())
	## treat last
	sequences.append((scpos,side,seqs))
	fic.close()
	return sequences

def ReadPartial(fic):
	sequences = []
	scpos_current = ""
	side_current  = ""
	nbi_current   = ""
	seqs = []
	fic.next()
	for line in fic:
		chrom,pos,side,nbident,seq = line.strip().split("\t")
		if side  == "5'":
			newpos = int(pos) - int(nbident)
		else:
			newpos = int(pos) + int(nbident)
		scpos = "%s:%s" % (chrom,newpos)
		if scpos == scpos_current and side == side_current and nbident == nbi_current:
			seqs.append(seq)
		else:
			sequences.append((scpos_current,side_current,seqs,nbi_current))
			seqs = [seq]
			scpos_current = scpos
			side_current  = side
			nbi_current   = nbident
	return sequences

def TryAddSequence(seq,nr,cc,side):
	flag = 0
	maxscore = 2*len(seq)
	if side == "5'":
		aln_scores = [PairwiseScore(seq[::-1],cc[i][0][::-1])*100/maxscore for i in range(len(cc))]
	else:
		aln_scores = [PairwiseScore(seq,cc[i][0])*100/maxscore for i in range(len(cc))]
	maxi = max(aln_scores)
	if maxi > 90:
		idx = aln_scores.index(maxi)
		cc[idx] = (cc[idx][0],nr+cc[idx][1])
		flag = nr
	return flag,cc

def ComputeConsensus(seqs,longuest):
	consensus = ""
	tmp = [longuest]
	tmp.extend(seqs)
	for i in range(len(longuest)):
		nt = [s[i] for s in tmp if i < len(s)]
		if len(nt) == 0:
			continue
		consensus += max(set(nt), key=nt.count)
	return consensus

def PairwiseScore(seq1, seq2):
	common = 0
	for s1, s2 in zip(seq1, seq2):
		if s1 == s2:
			common += 1
	return common*2

def SearchConsensusPairwise(seqs,side):
	uniq_seqs = list(set(seqs))
	if len(uniq_seqs) == 1:
		return uniq_seqs[0], len(seqs), []
	sort_seqs = sorted(uniq_seqs, key=lambda s:len(s),reverse=True)
	longuestSeq = sort_seqs.pop(0)
	max_scores = [2*len(s) for s in sort_seqs]
	if (side == "5'"): ## align from right to left
		revlonguest = longuestSeq[::-1]
		aln_scores = [PairwiseScore(revlonguest,sort_seqs[i][::-1])*100/max_scores[i] for i in range(len(sort_seqs))]
		kept_seqs  = [sort_seqs[i] for i in range(len(aln_scores)) if aln_scores[i] >= 90]
		revkept    = [i[::-1] for i in kept_seqs]
		consensus   = ComputeConsensus(revkept,revlonguest)[::-1]
	elif(side == "3'"): ## align from left to right
		aln_scores = [PairwiseScore(longuestSeq,sort_seqs[i])*100/max_scores[i] for i in range(len(sort_seqs))]
		kept_seqs  = [sort_seqs[i] for i in range(len(aln_scores)) if aln_scores[i] >= 90]
		consensus = ComputeConsensus(kept_seqs,longuestSeq)
	seqs_last = seqs[:]
	seqs_last = filter(lambda a: a!= longuestSeq, seqs_last)
	for s in  kept_seqs:
		seqs_last = filter(lambda a: a!= s, seqs_last)
	nb = len(seqs) - len(seqs_last)
	return consensus,nb,seqs_last


def TreatSequences(data,maxsc):
	name = data[0]
	side = data[1]
	seqs = data[2]
	outputs = ("","",[])
	if len(seqs) == 1:
		if len(seqs[0]) >= maxsc:
			cons = seqs[0]
			outputs = (name,side,[(cons,1)])
	else:
		cons_list = []
		while len(seqs) > 0:
			cons,nb,seqs = SearchConsensusPairwise(seqs,side)
			if len(cons) >= maxsc:
				cons_list.append((cons,nb))
		outputs = (name,side,cons_list)
	if len(data) == 4:
		outputs = (outputs[0],outputs[1],outputs[2],data[3])
	return outputs

def func_star(arg):
	return TreatSequences(*arg)


def main():
	parser = argparse.ArgumentParser(description='Compute soft-clip sequence consensus on each position')
	parser.add_argument("-i",action="store", type=file,                   required=True, dest="input",   help="txt summary of mappable sc",metavar="filename")
	parser.add_argument("-j",action="store", type=file,                   required=True, dest="partial", help="txt summary  partially equal to reference",metavar="filename")
	parser.add_argument("-f",action="store", type=file,                   required=True, dest="fasta",   help="fasta file of sc sequence",metavar="filename")
	parser.add_argument("-o",action="store", type=argparse.FileType('w'), required=True, dest="outtxt",  help="new summary",metavar="filename")
	parser.add_argument("-m",action="store", type=int, 	default=0,  	  required=False, dest="maxsc",  help="lg min to consider a consensus ",metavar="number")
	parser.add_argument("-c",action="store", type=int, 					  required=True, dest="nbtread", help="nb cpu to use",metavar="number")
	args = parser.parse_args()

	## TREAT GOOD SC SEQUENCES
	sequences = ReadFasta(args.fasta)  ## [(chrom:pos, side, [seqs])]
	pool      = mp.Pool(processes=args.nbtread)
	results   = pool.map(func_star, itertools.izip(sequences,itertools.repeat(args.maxsc))) ## [(chrom:pos, side, [(cons,nr)])]
	pool.close()
	pool.join()
	del sequences
	## TREAT PARTIALLY SC SEQUENCES
	sequences    = ReadPartial(args.partial)  ## [(chrom:newpos, side, [seqs], nbident)]
	pool         = mp.Pool(processes=args.nbtread)
	resu_partial = pool.map(func_star, itertools.izip(sequences,itertools.repeat(args.maxsc))) ## [(chrom:newpos, side, [(cons,nr)], nbident)]
	pool.close()
	pool.join()
	del sequences
	
	## Sort results of good sc
	consensus = {}  ## chrom:pos : side : [(cons,nr)]
	for name,side,cons_list in results:
		if consensus.has_key(name):
			consensus[name][side] = cons_list
		else:
			consensus[name] = {side : cons_list}
	del results
	## ADD partially aligned on reference sc pos (add only those that complete a consensus already defined)
	nbAdded = 0
	nbAll   = 0
	for name,side,cons_list,nbident in resu_partial:
		nbAll += sum([nr for seq,nr in cons_list])
		try:
			cc = consensus[name][side]
			all_flags = []
			for seq,nr in cons_list:
				nbi = int(nbident)
				if side == "5'":
					newseq = seq[0:-nbi]
				else:
					newseq = seq[nbi:]
				flag,cc = TryAddSequence(newseq,nr,cc,side) # try to add to existing consensus
				if flag == 0 and nr > 1: # add new consensus if more than 1 read
					flag = nr
					cc.append((seq,nr))
				nbAdded += flag
			consensus[name][side] = cc
		except KeyError: # add new position with consensus supported by at least 2 reads
			tosave = [i for i in cons_list if i[1] > 1]
			if len(tosave) > 0:
				try: # other side already save for example
					consensus[name][side] = tosave
				except KeyError:
					consensus[name] = {side : tosave}
				nbAdded += sum([nr for seq,nr in tosave])
				
	print "soft-clip sequences to rescue : %s (rescued : %s)" % (nbAll,nbAdded)

	header = args.input.readline()
	args.outtxt.write("chrom\tpos\tside\tnum\tnbr\tmaxsc\tseq\n")
	nb_wo = 0
	for line in args.input:
		l = line.strip().split()
		scpos = l[0] + ":" + l[1]
		side = l[2]
		try:
			cons_list = consensus[scpos][side]
			num = 1
			nbr_all = 0
			for cons,nb in cons_list:
				args.outtxt.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (l[0],l[1],side,num,nb,len(cons),cons))
				num += 1
				nbr_all += nb
		except KeyError:
			nb_wo += 1
			print "No consensus found",scpos,side
			continue
	
	args.input.close()
	args.outtxt.close()


if __name__=="__main__":
    main()
