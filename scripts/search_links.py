#!/home/anneso/myroot/bin/python
# -*- encoding: utf8 -*-

## CREATED 12/08/2015 by anso

## search links between sc pos and abnormal clusters (based on snooper_v2/search_links_2.py)

## MODIFICATIONS
	## ParseCluster() : remove header line skip and add only overlapping clusters
	## Modify algorithm:
		## search each sc pos independantly (even if candidates exist on both side)
		## then search best pair if many solutions (by using start of alignment)
		## remove ADD reference sequence around sc pos
## MODIFICATIONS 9/10/2017
	## PLUG algo 2.3 on 2.4


import argparse
import re
from math import log
from Bio import cpairwise2
from Bio import pairwise2
from pyfaidx import Fasta
import multiprocessing as mp
from itertools import izip, repeat

pairwise2.MAX_ALIGNMENTS=1

DELTA = 0
REFERENCE = ""

######################
## CLASS DEFINITION
######################

class SCPOS(object):
	def __init__(self,c,p,s):
		self.chrom = c
		self.pos   = p
		self.side  = s
		self.clusters = []
		self.consensus = {}
		self.flagali = {}
		self.reali = {}
		self.reference = ""
		self.nbcons = 0
		self.mappable = ""
	def GetName(self):
		return "%s:%s:%s" % (self.chrom, self.pos,self.side)
	def AddClus(self,name):
		self.clusters.append(name)
	def AddConsensus(self,numC, nb,lgC,seq,mapp,reali):
		self.nbcons += 1
		self.consensus[numC] = (seq,nb)
		self.flagali[numC] = reali
		self.mappable = mapp
	def AddReali(self,num,chrom,pos,side,dist):
		self.reali[num] = [chrom,pos,side,dist]
	def AddReference(self,ref):
		self.reference = ref
	def __repr__(self):
		return "chr %s : %s %s (nbclus = %s)\n" % (self.chrom,self.pos,self.side,len(self.clusters))
	def __str__(self):
		return "chr %s : %s %s (nbclus = %s)\n" % (self.chrom,self.pos,self.side,len(self.clusters))


class CLUSTER(object):
	def __init__(self, l):
		self.type_ab = l[0]
		self.num_ab  = int(l[1])
		self.chr1    = l[2]
		self.chr2    = l[5]
		self.start1  = int(l[3])
		self.start2  = int(l[6])
		self.end1    = int(l[4])
		self.end2    = int(l[7])
		self.nbp     = int(l[8]) + int(l[9])
		self.nbo     = int(l[9])
		self.insert  = int(l[10])
		self.unmap    = l[11]
		self.germline = l[12]
		self.strand   = self.type_ab[-2:]
		self.sc1_candidates = []
		self.sc2_candidates = []
	def GetName(self):
		return "%s%s" % (self.type_ab , self.num_ab)
	def AddCandidates(self,scpos,ab):
		if ab == "a":
			self.sc1_candidates.append(scpos)
		else:
			self.sc2_candidates.append(scpos)
	def GetClusPos(self,ab):
		if ab == "a" and self.strand in ["FF","FR"]:
			return self.end1
		elif ab == "a" and self.strand in ["RR","RF"]:
			return self.start1
		elif ab == "b" and self.strand in ["FF","RF"]:
			return self.end2
		else:
			return self.start2
	def GetStrand(self,ab):
		if (ab == "a" and self.strand in ["FF","FR"]) or (ab == "b" and self.strand in ["FF","RF"]):
			return "F"
		else:
			return "R"
	def GetChrom(self,ab):
		if ab == "a":
			return self.chr1
		else:
			return self.chr2
	def IsInverted(self):
		if self.strand in ["FF","RR"]:
			return True
		else:
			return False
	def __repr__(self):
		return "%s %s [chr %s : %s-%s] - [chr %s : %s-%s] (nbp = %s)" % (self.type_ab,self.num_ab,self.chr1,self.start1,self.end1,self.chr2,self.start2,self.end2,self.nbp)
	def __str__(self):
		return "%s %s [chr %s : %s-%s] - [chr %s : %s-%s] (nbp = %s)" % (self.type_ab,self.num_ab,self.chr1,self.start1,self.end1,self.chr2,self.start2,self.end2,self.nbp)



######################
## INPUTS PARSING - OUTPUT WRITING
######################
def ParseClusters(fic):
	clusters = {}  # name : clus object
	i = 0
	fic.next()
	for line in fic:
		i += 1
		l = line.strip().split("\t")
		name = l[0]+l[1]
		clus = CLUSTER(l)
		clusters[name] = clus
	print "\t%s\tclusters kept" % i
	fic.close()
	return clusters

def TestOrientation(side,strand):
	if (strand == "F" and side == "3'") or (strand == "R" and side == "5'"):
		return True
	else:
		return False

def ParseOverlapClusSc(fic,clusters,sc_defined):
	nb_rejected  = 0
	nb_all       = 0
	for line in fic:
		nb_all += 1
		l        = line.strip().split()
		clusname = l[3][0:-1]
		clusAB   = l[3][-1]
		side     = l[7]
		nameSC   = ":".join([l[0],l[6],side])
		clus     = clusters[clusname]
		strand   = clus.GetStrand(clusAB)
		if TestOrientation(side,strand):
			clus.AddCandidates(sc_defined[nameSC],clusAB)
		else:
			nb_rejected += 1
	fic.close()
	print "\t%s\toverlapping lines read" % nb_all
	print "\t%s\toverlaps rejected" % nb_rejected
	return

def ParseScPos(fic):
	fic.next()
	i = 0
	sc_defined = {}
	for line in fic:
		i += 1
		l      = line.strip().split()
		nameSc = ":".join(l[0:3])
		if sc_defined.has_key(nameSc):
			sc_defined[nameSc].AddConsensus(l[3],l[4],l[5],l[7],l[8],l[6]) #num nr	lg seq map flag
		else:
			scpos  = SCPOS(l[0],int(l[1]),l[2])
			scpos.AddConsensus(l[3],l[4],l[5],l[7],l[8],l[6])
			sc_defined[nameSc] = scpos		
	fic.close()
	print "\t%s\tconsensus added" % i
	return sc_defined

def AddRealiInfo(fic,sc_defined):
	fic.next()
	i = 0
	for line in fic:
		l = line.strip().split("\t")
		chrom,pos,side,numCons = l[0].split(":")
		try:
			sc_defined[":".join([chrom,pos,side])].AddReali(numCons,l[2],l[3],l[4],l[5])
			i += 1
		except KeyError:
			pass
	fic.close()
	print "\t%s\trealignment added" % i
	return


def WriteSVoutput(cluster,scpos1,scpos2,cons1,cons2,score1,score2,aln1Kept,aln2Kept,filterValue):
	if cons1 == 0:
		seq1 = ""
		nr1 = 0
	else:
		seq1 = scpos1.consensus[cons1][0]
		nr1 = scpos1.consensus[cons1][1]
	if cons2 == 0:
		seq2 = ""
		nr2 = 0
	else:
		seq2 = scpos2.consensus[cons2][0]
		nr2 = scpos2.consensus[cons2][1]
	scmap = "%s,%s" % (scpos1.mappable,scpos2.mappable)
	outline = "%s\t"*26 + "%s\n"
	return outline % (cluster.type_ab, cluster.num_ab, cluster.chr1, scpos1.pos, cluster.chr2, scpos2.pos,
			   cluster.start1,cluster.end1,cluster.start2, cluster.end2, cluster.insert, cluster.unmap, cluster.germline,
			   cluster.nbp,cluster.nbo, nr1, nr2,scmap,filterValue, score1, score2, aln1Kept, aln2Kept, len(seq1),len(seq2),seq1,seq2)

#0  type      num        chr1       sc1     chr2    sc2     start1  end1        start2  end2
#10 insert    ab_unmap   ab_germ    nbpall  nbo     nr1     nr2     sc_unmapp   filter  scoreSC1
#20 scoreSC2  scoreAli1  scoreAli2  lgcons1 lgcons2 cons1   cons2


###########################
## SEARCH FOR SV EVENTS
###########################

def ReverseTranslate(seq):
	seq = seq[::-1]
	new = ""
	tr = {"A":"T","T":"A", "C":"G","G":"C","N":"N"}
	for i in seq:
		try:
			new = new + tr[i]
		except KeyError:
			new = new + i
	return new

def ComputeScore(scpos,cluster,ab):
	dist = abs(cluster.GetClusPos(ab) - scpos.pos)
	allScore = {}
	log_text = ""
	for cons_num in scpos.consensus.keys():
		cons = scpos.consensus[cons_num][0]
		nr = int(scpos.consensus[cons_num][1])
		maxsc = len(cons)
		bonus_lg = 0
		if maxsc >= 10 and maxsc <=30:
			bonus_lg = 30
		elif maxsc > 30:
			bonus_lg = 50
		score = max(0,100-(100/DELTA)*dist) + min(100,int(round(25*log(nr)))) + bonus_lg
		log_text += "\t\tSCORE : %s:%s (%s , nbr=%s, dist=%s, maxsc=%s)\n" % (scpos.GetName(),cons_num,score,nr,dist,maxsc)
		allScore[cons_num] = score
	return allScore,log_text


def MakeAln(cons,ref,pos):
	score = 0
	start = -1
	aln  = pairwise2.align.localms(ref,cons,2,-1,-4,-2)
	if len(aln) > 0:
		ali   = aln[0]
		score = round(ali[2]*100.0/(len(cons)*2),1)
		start = pos - 150 + len(re.split("[A-Z]*",ali[1])[0])
	return score,start
	

def ComputeAlnScore(scpos,cluster,ab):
	log_text = ""
	nameSC   = scpos.GetName()
	clus_pos = cluster.GetClusPos(ab)
	chrom    = cluster.GetChrom(ab)
	if cluster.IsInverted():
		ref = ReverseTranslate(REFERENCE[chrom][(clus_pos-150):(clus_pos+150)].seq.upper())
	else:
		ref = REFERENCE[chrom][(clus_pos-150):(clus_pos+150)].seq.upper()
	consGood = {}
	for num in scpos.consensus.keys():
		flagali = scpos.flagali[num]
		cons    = scpos.consensus[num][0]
		if flagali in ["not_realigned","reali_badali","reali_badchrom","reali_unaligned"]:
			score,startAli = MakeAln(cons,ref,clus_pos)
			dist = abs(clus_pos - startAli)
			if score >= 80:
				consGood[num] = (score,startAli)
				log_text += "\t\tALN   : %s:%s (%s , dist=%s, start=%s)\n" % (nameSC, num, score, dist, startAli)
			elif startAli != -1:
				log_text += "\t\tALN   : %s:%s (%s , dist=%s, start=%s , NOT KEPT)\n" % (nameSC, num, score, dist, startAli)
			else:
				log_text += "\t\tALN   : %s:%s , (no alignment, NOT KEPT)\n" % (nameSC, num)
		elif flagali in ["reali_del","reali_ins","reali_inter","reali_intra","reali_itself"]:
			align = scpos.reali[num] # chr_ali   pos_ali   side_ali   dist
			dist_bwa  = abs(clus_pos - int(align[1]))
			if chrom == align[0] and dist_bwa < 150:
				consGood[num] = (100,int(align[1]))
				log_text += "\t\tBWA   : %s:%s (100 , dist=%s , start=%s)\n" % (nameSC, num, dist_bwa, align[1])
			else:
				score,startAli = MakeAln(cons,ref,clus_pos)
				dist_aln = abs(clus_pos - startAli)
				if score >= 80:
					consGood[num] = (score,startAli)
					log_text += "\t\tBWA   : %s:%s (dist=%s, NOT KEPT)\t\tALN   : (%s , dist=%s, start=%s)\n" % (nameSC, num, dist_bwa, score, dist_aln, startAli)
				elif startAli != -1:
					log_text += "\t\tBWA   : %s:%s (dist=%s, NOT KEPT)\t\tALN   : (%s , dist=%s, start=%s , NOT KEPT)\n" % (nameSC, num, dist_bwa, score, dist_aln, startAli)
				else:
					log_text += "\t\tBWA   : %s:%s (dist=%s, NOT KEPT)\t\tALN   : (no alignment, NOT KEPT)\n" % (nameSC, num, dist_bwa)
		else:
			log_text += "\t\t BAD FLAG ALI : %s : %s:%s" % (flagali,nameSC,num)
	return consGood,log_text



#######################
def SearchBestCons(ali_score,sc_score,scpos):
	consKept = 0
	score = -1
	if len(ali_score) == 1:
		num = ali_score.keys()[0]
		return num
	else:
		for num in ali_score:
			curr_score    = sc_score[num]
			if curr_score > score:
				score    = curr_score
				consKept = num
		return consKept

def SearchSCpos(sc1_scores,sc2_scores,sc1_aligned,sc2_aligned,sc_defined):
	gscoreMax = -1
	sc1Kept = ""
	sc2Kept = ""
	cons1Kept = 0
	cons2Kept = 0
	pairs_scores = {}
	text = ""
	for name1 in sc1_aligned:
		scpos1  = sc_defined[name1]
		for name2 in sc2_aligned:
			scpos2 = sc_defined[name2]
			text += "\n\t\tScore ( %s with %s ):  " % (name1,name2)
			for num1 in sc1_aligned[name1]:
				sc1_sc    = sc1_scores[name1][num1]
				sc1_alisc = sc1_aligned[name1][num1][0]
				sc1_start = sc1_aligned[name1][num1][1]
				for num2 in sc2_aligned[name2]:
					sc2_sc    = sc2_scores[name2][num2]
					sc2_alisc = sc2_aligned[name2][num2][0]
					sc2_start = sc2_aligned[name2][num2][1]
					dist1 = abs(scpos1.pos - sc2_start)
					dist2 = abs(scpos2.pos - sc1_start)
					g_score = sc1_sc + sc1_alisc + sc2_sc + sc2_alisc + max(0,20-dist1) + max(0,20-dist2)
					text += "(%svs%s = %s) " % (num1,num2,g_score)
					if g_score > gscoreMax:
						gscoreMax = g_score
						sc1Kept = name1
						sc2Kept = name2
						cons1Kept = num1
						cons2Kept = num2
					elif g_score == gscoreMax and (sc1_alisc + sc2_alisc > sc1_aligned[sc1Kept][cons1Kept][0] + sc2_aligned[sc2Kept][cons2Kept][0]):
						gscoreMax = g_score
						sc1Kept = name1
						sc2Kept = name2
						cons1Kept = num1
						cons2Kept = num2
	text  = text.lstrip("\n")
	text += "\n\tFINAL CHOICE : %s with %s (cons1=%s ; cons2=%s)\n" % (sc1Kept,sc2Kept,cons1Kept,cons2Kept)
	return sc1Kept,cons1Kept,sc2Kept,cons2Kept,text

def SearchSCpos1Side(sc_scores,aln_scores,sc_defined):
	bestScore = -1
	logtxt = ""
	for sc in aln_scores:
		num_cons_kept  = SearchBestCons(aln_scores[sc],sc_scores[sc],sc_defined[sc])
		scoreAlnKept   = aln_scores[sc][num_cons_kept][0]
		logtxt += "\t\t%s:\tbest is cons%s : (sc_score=%s, ali_score=%s, maxsc=%s)\n" % \
			(sc,num_cons_kept, sc_scores[sc][num_cons_kept],scoreAlnKept,len(sc_defined[sc].consensus[num_cons_kept][0]))
		if scoreAlnKept > bestScore:
			bestScore = scoreAlnKept
			scKept = sc
			consKept = num_cons_kept
	logtxt += "\tFINAL CHOICE : %s\n" % (scKept)
	return scKept,consKept,logtxt

def AbScLinkage(cluster):
	log = "\n" + str(cluster) + "\n"
	sc_NULL = SCPOS("",-1,"")
	sc_NULL.AddConsensus(0,0,0,"-","","")
	
	## COMPUTE SCORES
	flaga = 0
	flagb = 0
	if len(cluster.sc1_candidates) > 0:
		flaga = 1
	if len(cluster.sc2_candidates) > 0:
		flagb = 1
	
	sc_used = {}
	if flaga == 1:
		nbc_sc1 = sum([scpos.nbcons for scpos in cluster.sc1_candidates])
		sc1_scores = {}   # nameSC : cons : score
		sc1_aligned = {}  # nameSC : cons : (score, alistart)
		log += "\tCOMPUTE SC1 candidates score (%s scpos ; %s cons):\n" % (len(cluster.sc1_candidates), nbc_sc1)
		for scpos in cluster.sc1_candidates:
			sc_used[scpos.GetName()] = scpos
			sc1_scores[scpos.GetName()],  logtxt1 = ComputeScore(   scpos,cluster,"a")
			sc1_aligned[scpos.GetName()], logtxt2 = ComputeAlnScore(scpos,cluster,"b")
			log += logtxt1 + logtxt2
		sc1_aligned = {i:j for i,j in sc1_aligned.iteritems() if len(j) > 0}
	if flagb == 1:
		nbc_sc2 = sum([scpos.nbcons for scpos in cluster.sc2_candidates])
		sc2_scores = {}
		sc2_aligned = {}
		log += "\tCOMPUTE SC2 candidates score (%s scpos ; %s cons):\n" % (len(cluster.sc2_candidates), nbc_sc2)
		for scpos in cluster.sc2_candidates:
			sc_used[scpos.GetName()] = scpos
			sc2_scores[scpos.GetName()],  logtxt1 = ComputeScore(   scpos,cluster,"b")
			sc2_aligned[scpos.GetName()], logtxt2 = ComputeAlnScore(scpos,cluster,"a")
			log += logtxt1 + logtxt2
		sc2_aligned = {i:j for i,j in sc2_aligned.iteritems() if len(j) > 0}
	
	## SEARCH BEST SV CANDIDATE
	if flaga + flagb == 0:
		out = WriteSVoutput(cluster,sc_NULL,sc_NULL,0,0,0,0,0,0,"REJECTED_no.SC.candidates")
	
	elif flaga + flagb == 2:
		if len(sc1_aligned) > 0 and len(sc2_aligned) > 0:
			log += "\tSEARCH BEST CONS BOTH SIDE:\n"
			sc1Kept,numcons1Kept,sc2Kept,numcons2Kept,logtxt = SearchSCpos(sc1_scores,sc2_scores,sc1_aligned,sc2_aligned,sc_used)
			log += logtxt
			if sc1Kept == "" and sc2Kept == "":
				out = WriteSVoutput(cluster,sc_NULL,sc_NULL,0,0,0,0,0,0,"REJECTED_both")
			elif sc1Kept == "":
				out = WriteSVoutput(cluster,sc_NULL,sc_used[sc2Kept],0,numcons2Kept,0,sc2_scores[sc2Kept][numcons2Kept],
									0,sc2_aligned[sc2Kept][numcons2Kept][0],"PASS_both.onlySC2")
			elif sc2Kept == "":
				out = WriteSVoutput(cluster,sc_used[sc1Kept],sc_NULL,numcons1Kept,0,sc1_scores[sc1Kept][numcons1Kept],
									0,sc1_aligned[sc1Kept][numcons1Kept][0],0,"PASS_both.onlySC1")
			else:
				out = WriteSVoutput(cluster,sc_used[sc1Kept],sc_used[sc2Kept],numcons1Kept,numcons2Kept,
							  sc1_scores[sc1Kept][numcons1Kept],sc2_scores[sc2Kept][numcons2Kept],
							  sc1_aligned[sc1Kept][numcons1Kept][0],sc2_aligned[sc2Kept][numcons2Kept][0],"PASS_both")
		elif len(sc1_aligned) > 0 :
			sc1Kept,numcons1Kept,logtxt = SearchSCpos1Side(sc1_scores,sc1_aligned,sc_used)
			log += "\tSEARCH BEST CONS ONLY SIDE A:\n" + logtxt
			out = WriteSVoutput(cluster,sc_used[sc1Kept],sc_NULL,numcons1Kept,0,sc1_scores[sc1Kept][numcons1Kept],
								0,sc1_aligned[sc1Kept][numcons1Kept][0],0,"PASS_onlySC1")
		elif len(sc2_aligned) > 0 :
			sc2Kept,numcons2Kept,logtxt = SearchSCpos1Side(sc2_scores,sc2_aligned,sc_used)
			log += "\tSEARCH BEST CONS ONLY SIDE B:\n" + logtxt
			out = WriteSVoutput(cluster,sc_NULL,sc_used[sc2Kept],0,numcons2Kept,0,sc2_scores[sc2Kept][numcons2Kept],
								0,sc2_aligned[sc2Kept][numcons2Kept][0],"PASS_onlySC2")
		else:
			log += "\t\tNO ALIGNMENT ON BOTH SIDE\n"
			out = WriteSVoutput(cluster,sc_NULL,sc_NULL,0,0,0,0,0,0,"REJECTED_both")
	
	elif flaga == 1:
		if len(sc1_aligned) > 0:
			sc1Kept,numcons1Kept,logtxt = SearchSCpos1Side(sc1_scores,sc1_aligned,sc_used)
			log += "\tSEARCH BEST CONS SIDE A:\n" + logtxt
			out = WriteSVoutput(cluster,sc_used[sc1Kept],sc_NULL,numcons1Kept,0,sc1_scores[sc1Kept][numcons1Kept],
								0,sc1_aligned[sc1Kept][numcons1Kept][0],0,"PASS_onlySC1")
		else:
			log += "\t\tNO ALIGNMENT SIDE A\n"
			out = WriteSVoutput(cluster,sc_NULL,sc_NULL,0,0,0,0,0,0,"REJECTED_onlySC1")
	
	elif flagb == 1:
		if len(sc2_aligned) > 0:
			sc2Kept,numcons2Kept,logtxt = SearchSCpos1Side(sc2_scores,sc2_aligned,sc_used)
			log += "\tSEARCH BEST CONS SIDE B:\n" + logtxt
			out = WriteSVoutput(cluster,sc_NULL,sc_used[sc2Kept],0,numcons2Kept,0,sc2_scores[sc2Kept][numcons2Kept],
								0,sc2_aligned[sc2Kept][numcons2Kept][0],"PASS_onlySC2")
		else:
			log += "\t\tNO ALIGNMENT SIDE B\n"
			out = WriteSVoutput(cluster,sc_NULL,sc_NULL,0,0,0,0,0,0,"REJECTED_onlySC2")
	return {cluster.GetName() : [out,log]}
	


###########################
## MAIN PROGRAM
###########################



def main():
	parser = argparse.ArgumentParser(description='First filtering of germline abnormal pairs')
	parser.add_argument("-i",action="store", type=file, required=True, dest="overlap",  help="bed file of overlaps between sc pos and clusters of abnormal pairs",metavar="file")
	parser.add_argument("-c",action="store", type=file, required=True, dest="cluster",  help="text file describing clusters of abnormal pairs",metavar="file")
	parser.add_argument("-s",action="store", type=file, required=True, dest="consensus",help="consensus file of sc positions",metavar="dir")
	parser.add_argument("-a",action="store", type=file, required=True, dest="reali",    help="realignment file",metavar="dir")
	parser.add_argument("-r",action="store", type=str, required=True, dest="reference",  help="ref genome in fasta",metavar="dir")
	parser.add_argument("-o",action="store", type=argparse.FileType('w'), required=True, dest="out",  help="output filename",metavar="dir")
	parser.add_argument("-l",action="store", type=argparse.FileType('w'), required=True, dest="log",  help="output log filename",metavar="dir")
	parser.add_argument("-d",action="store", type=int, required=True, dest="delta", help="delta used around sc positions",metavar="dir")
	parser.add_argument("-t",action="store", type=int, required=True, dest="nbtread", help="nb parallel thread",metavar="dir")
	args = parser.parse_args()
	
	global DELTA
	global REFERENCE
	DELTA = args.delta
	REFERENCE = Fasta(args.reference)

	print "Reading input data : "
	clustersAb     = ParseClusters(args.cluster)
	sc_defined     = ParseScPos(args.consensus)
	AddRealiInfo(args.reali,sc_defined)
	ParseOverlapClusSc(args.overlap,clustersAb,sc_defined)
	# sc_defined      = nameSC  : SCPOS
	# clustersAb      = name : CLUSTER
	# nbpclus         = nbp : name
	del sc_defined
	
	print "Start linkage between abnormal clusters and soft-clipped positions"
	if args.nbtread > 1:
		pool    = mp.Pool(processes=args.nbtread)
		results = pool.map(AbScLinkage, clustersAb.values())
		pool.close()
		pool.join()
	else:
		results = []
		for cluster in clustersAb.values():
			results.append(AbScLinkage(cluster))

	args.out.write("type\tnum\tchr1\tsc1\tchr2\tsc2\tstart1\tend1\tstart2\tend2\tinsert\tab_unmap\tab_germ\
\tnbpall\tnbo\tnr1\tnr2\tsc_unmapp\tfilter\tscoreSC1\tscoreSC2\tscoreAli1\tscoreAli2\tlgcons1\tlgcons2\tcons1\tcons2\n")
	for i in results:
		tmp = i.copy()
		clusname = tmp.keys()[0]
		args.out.write(tmp[clusname][0])
		args.log.write(tmp[clusname][1])
	args.out.close()
	args.log.close()
	return 0

main()

