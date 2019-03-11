#!/data-ddn/software/Python/current/python
# -*- encoding: utf8 -*-

## CREATED 03/11/2014 by anso

## search links between sc pos and abnormal clusters (based on search_links_v2.py)
## Add better definition of SV => combine oriRR/oriFF to define inversion when possible etc ...

## MODIFICATIONS:
	## 4/09/2015 : undersrtand and treat CLUSTERING WARNING
	##			   abnormal clustering is not perfect, some pairs of blocs overlaps thus
	##			   they overlap the same oriRR pair and are associated to the same
	##				 => only one association is done

import argparse
import re
from math import log


chrom_sorted = [str(i) for i in range(1,23)]
chrom_sorted.extend(["X","Y"])

######################
## CLASS DEFINITION
######################
#0  type      num        chr1       sc1     chr2    sc2     start1  end1        start2  end2
#10 insert    ab_unmap   ab_germ    nbpall  nbo     nr1     nr2     sc_unmapp   filter  scoreSC1
#20 scoreSC2  scoreAli1  scoreAli2  lgcons1 lgcons2 cons1       cons2
class SV(object):
	def __init__(self, l):
		self.type_ab = l[0]
		self.num     = l[1]
		self.chr1    = l[2]
		self.chr2    = l[4]
		self.sc1     = int(l[3])
		self.sc2     = int(l[5])
		self.start1  = int(l[6])
		self.end1    = int(l[7])
		self.start2  = int(l[8])
		self.end2    = int(l[9])
		self.nbp     = int(l[13])
		self.nbo     = int(l[14])
		self.nr1     = int(l[15])
		self.nr2     = int(l[16])
		self.lg1     = int(l[23])
		self.lg2     = int(l[24])
		if self.sc1 == -1:
			if self.type_ab in ["intraFR","intraFF","interFF","interFR"]:
				self.sc1 = int(l[7])
			elif self.type_ab in ["intraRF","intraRR","interRR","interRF"]:
				self.sc1 = int(l[6])
		if self.sc2 == -1:
			if self.type_ab in ["intraFR","intraRR","interRR","interFR"]:
				self.sc2 = int(l[8])
			elif self.type_ab in ["intraRF","intraFF","interFF","interRF"]:
				self.sc2 = int(l[9])
		self.ABunmap = l[11]
		self.ABgerm  = l[12]
		self.SCunmap = l[17]
		self.descr = "%s;%s(%s),%s(%s);%s(%s,%s),%s-%s,%s-%s;%s,%s,%s" % \
		(self.num,self.nr1,self.lg1,self.nr2,self.lg2,self.nbp,self.nbp-self.nbo,self.nbo,self.start1,self.end1,self.start2,self.end2,self.ABunmap,self.ABgerm,self.SCunmap)
	def merge(self,sv):
		self.type_ab += "RR"
		unmapAB = set()
		unmapSC = set()
		if self.ABunmap != "none":
			unmapAB = unmapAB.union(set(list(self.ABunmap)))
		if sv.ABunmap != "none":
			unmapAB = unmapAB.union(set(list(sv.ABunmap)))
		if unmapAB == set():
			unmapAB = "none"
		else:
			unmapAB = "".join(sorted(list(unmapAB)))
		if self.SCunmap != "none":
			unmapSC = unmapSC.union(set(list(self.SCunmap)))
		if sv.SCunmap != "none":
			unmapSC = unmapSC.union(set(list(sv.SCunmap)))
		if unmapSC == set():
			unmapSC = "none"
		else:
			unmapSC = "".join(sorted(list(unmapSC)))
		germline = set()
		if self.ABgerm != "none":
			germline = germline.union(set(list(self.ABgerm)))
		if sv.ABgerm != "none":
			germline = germline.union(set(list(sv.ABgerm)))
		if germline == set():
			germline = "none"
		else:
			germline = "".join(sorted(list(germline)))
		new = "%s;%s(%s),%s(%s);%s(%s,%s),%s-%s,%s-%s;%s,%s,%s" % \
		(self.num, self.nr1+sv.nr1, max(self.lg1,sv.lg1), self.nr2+sv.nr2, max(self.lg2,sv.lg2),
		 self.nbp + sv.nbp, self.nbp + sv.nbp - self.nbo -sv.nbo ,self.nbo + sv.nbo,
		 min(self.start1,sv.start1), max(self.end1,sv.end1), min(self.start2,sv.start2), max(self.end2,sv.end2),
		 unmapAB,germline,unmapSC)
		self.descr =  new + ";FF;" + self.descr + ";RR;" + sv.descr
		return self.type_ab

######################
## INPUTS PARSING
######################

def ParseOverlap(fic,links):
	asso = {}  # nameFF : ab : (nameRR)
	nbpRR = {}
	### Read overlap file
	for line in fic:
		l = line.strip().split()
		nameFF = l[3]
		nameRR = l[8]
		abFF   = nameFF[-1]
		abRR   = nameRR[-1]
		if (abFF == abRR):
			nameFF = nameFF[0:-1]
			nameRR = nameRR[0:-1]
			nbpRR[nameRR] = int(l[9])
			if(asso.has_key(nameFF)):
				try:
					asso[nameFF][abFF].add(nameRR)
				except KeyError:
					asso[nameFF][abFF] = set([nameRR])
			else:
				asso[nameFF] = {abFF : set([nameRR])}
	fic.close()
	### Search associations
	trueLinks = []
	alreadyUsed = set()
	for nameFF in asso:
		if asso[nameFF].has_key("a") and asso[nameFF].has_key("b"):
			setA = asso[nameFF]["a"]
			setB = asso[nameFF]["b"]
			overlap = setA.intersection(setB)
			if len(overlap) == 1:
				overB = list(overlap)[0]
				if overB not in alreadyUsed:
					trueLinks.append((nameFF,overB))
					alreadyUsed.add(overB)
			elif len(overlap) > 1:
				candidat = list(overlap)
				bestRR = ""
				nbp = 0
				for currentRR in candidat:
					if currentRR in alreadyUsed:
						continue
					if nbpRR[currentRR] > nbp:
						bestRR = currentRR
						nbp = nbpRR[currentRR]
				if bestRR != "":
					trueLinks.append((nameFF,bestRR))
					alreadyUsed.add(bestRR)
	### Combine
	for nameFF, nameRR in trueLinks:
		svFF = links[nameFF]
		try:
			svRR = links[nameRR]
		except KeyError:
			print "CLUSTERING WARNING : ",nameRR
			continue
		del links[nameFF]
		del links[nameRR]
		name = svFF.merge(svRR)
		links[name + str(i)] = svFF
	return links


def ParseLinks(fic):
	data = {} # name : SV
	fic.next()
	for line in fic:
		l = line.strip("\n").split("\t")
		name = l[0]+l[1]
		data[name] = SV(l)
	fic.close()
	return data

def main():
	parser = argparse.ArgumentParser(description='First filtering of germline abnormal pairs')
	parser.add_argument("-i",action="store", type=file, required=True, dest="links",  help="file of sv reconstructed",metavar="file")
	parser.add_argument("-j",action="store", type=file, required=True, dest="overlap",  help="RR/FF overlaps",metavar="file")
	parser.add_argument("-g",action="store_true",required=False, dest="germline",help="germline case")
	parser.add_argument("-o",action="store", type=argparse.FileType('w'), required=True, dest="out",  help="output filename",metavar="dir")
	parser.add_argument("-f",action="store", type=argparse.FileType('w'), required=True, dest="filt",  help="output filename filtered",metavar="dir")
	args = parser.parse_args()

	links = ParseLinks(args.links)
	links = ParseOverlap(args.overlap,links)

	args.out.write("type\tchr1\tsc1\tchr2\tsc2\tdescription\n")
	args.filt.write("type\tchr1\tsc1\tchr2\tsc2\tdescription\n")
	for name in sorted(links):
		sv = links[name]
		outline = "%s\t%s\t%s\t%s\t%s\t%s\n" % (sv.type_ab,sv.chr1,sv.sc1,sv.chr2,sv.sc2,sv.descr)
		args.out.write(outline)
		if args.germline:
			if sv.ABunmap != "ab":
				if sv.type_ab in ["intraRF", "intraFR"] and sv.nbo == sv.nbp: ## only overlap
					if  sv.nr1 > 1 and sv.lg1 >= 15 and sv.nr2 > 1 and sv.lg2 >= 15:
						args.filt.write(outline)
				elif "intra" in sv.type_ab:
					if (sv.nr1 > 1 and sv.lg1 >= 15) or (sv.nr2 > 1 and sv.lg2 >= 15):
						args.filt.write(outline)
				else: # inter
					if sv.lg1 >= 20 and sv.lg2 >= 20 and (sv.nr1 > 2 or sv.nr2 > 2) and sv.nbp > 2:
						args.filt.write(outline)
		else:
			if sv.ABunmap != "ab" and sv.ABgerm == "none":
				if sv.type_ab in ["intraRF", "intraFR"] and sv.nbo == sv.nbp: ## only overlap
					if  sv.nr1 > 1 and sv.lg1 >= 15 and sv.nr2 > 1 and sv.lg2 >= 15:
						args.filt.write(outline)
				elif "intra" in sv.type_ab:
					if (sv.nr1 > 1 and sv.lg1 >= 15) or (sv.nr2 > 1 and sv.lg2 >= 15):
						args.filt.write(outline)
				else: # inter
					if sv.lg1 >= 20 and sv.lg2 >= 20 and (sv.nr1 > 2 or sv.nr2 > 2) and sv.nbp > 2:
						args.filt.write(outline)
	args.out.close()

	return 0

main()

