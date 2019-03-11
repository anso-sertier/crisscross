#!/usr/bin/env python
# -*- encoding: utf8 -*-

## CREATED 17/08/2015 by anso
## MODIFIED 9/09/2016

import argparse
import re


class OVERLAPS(object):
	def __init__(self,norm):
		self.map1 = set([])
		self.map2 = set([])
		if not norm :
			self.norm  = {"a" : set([]) , "b" : set([])}

def ReadInput(home,clus_obj):
	fic = open(home,"r")
	for line in fic:
		l = line.strip().split()
		tuclus = l[3][0:-1]
		tutype = l[3][-1]
		noclus = l[8][0:-1]
		notype = l[8][-1]
		if tutype == notype:
			if not clus_obj.has_key(tuclus):
				clus_obj[tuclus] = OVERLAPS(False)
			clus_obj[tuclus].norm[tutype].add(noclus)
	fic.close()
	return clus_obj

def ReadInputMappable(home,clus_obj,mapp,norm):
	fic = open(home,"r")
	for line in fic:
		l = line.strip().split()
		clus = l[3][0:-1]
		num = l[3][-1]
		if not clus_obj.has_key(clus):
			clus_obj[clus] = OVERLAPS(norm)
		if mapp == 1:
			clus_obj[clus].map1.add(num)
		else:
			clus_obj[clus].map2.add(num)
	fic.close()
	return clus_obj



def WriteOutput(list_germline, clus_obj, fic_clus, outspe, outunm, list_somatic):
	germline_clusters = {}
	header            = fic_clus.readline()
	outunm.write(header)
	outspe.write(header.strip("\n") + "\tunmappable\tgermline\n")
	for line in fic_clus:
		l = line.split()
		clusname = l[0] + l[1]
		## completely unmappable blocs
		try:
			if clus_obj[clusname].map1 == set(["a","b"]):
				outunm.write(line)
				if clusname in list_germline:
					germline_clusters[clusname] = [line.strip("\n").split(),"ab"]
				continue
		except KeyError:
			pass
		## partially unmappable blocs
		unmap = "none"
		try:
			unmap = "".join(sorted(list(clus_obj[clusname].map2)))
		except KeyError:
			pass
		## treat germline
		if clusname in list_germline:
			germline_clusters[clusname] = [line.strip("\n").split(),unmap]
		else:
			if list_somatic.has_key(clusname):
				value = list_somatic[clusname]
			else:
				value = "none"
			outspe.write("%s\t%s\t%s\n"   % (line.strip("\n"), unmap, value ))

	fic_clus.close()
	outspe.close()
	outunm.close()
	return germline_clusters


def WriteGermlines(list_germline, tuclus,noclus, tu_obj, no_obj, outgerm, outover):
	outgerm.write("type\tnumTU\tchrom1\tstart1\tend1\tchrom2\tstart2\tend2\tnbpab\tnbpov\tinsert\tunmappable\tgermline\tnumNO\tdescrTU\tdescrNO\n")
	for clusname in list_germline:
		normclus = list_germline[clusname].split(",")
		noinfos  = [noclus[i] for i in normclus]
		tuinfos  = tuclus[clusname]
		## totally unmappable blocs
		if tuinfos[1] == "ab" or "ab" in [i[1] for i in noinfos]:
			continue
		## partially unmappable blocs
		try:
			unmapTU = tu_obj[clusname].map2
		except KeyError:
			unmapTU = set([])
		unmapNO = set().union(*[no_obj[i].map2 for i in normclus if no_obj.has_key(i)])
		unmap   = "".join(sorted(list(unmapTU.union(unmapNO))))
		if unmap == "":
			unmap = "none"
		## treat germline
		nonum = ",".join([i[0][1] for i in noinfos])
		start1 = min(int(tuinfos[0][3]), min([int(i[0][3]) for i in noinfos]))
		end1   = max(int(tuinfos[0][4]), max([int(i[0][4]) for i in noinfos]))
		start2 = min(int(tuinfos[0][6]), min([int(i[0][6]) for i in noinfos]))
		end2   = max(int(tuinfos[0][7]), max([int(i[0][7]) for i in noinfos]))
		nbpab  = int(tuinfos[0][8]) + sum([int(i[0][8]) for i in noinfos])
		nbpov  = int(tuinfos[0][9]) + sum([int(i[0][9]) for i in noinfos])
		if clusname.startswith("inter"):
			insert = 0
		else:
			insert = start2 - end1
		descrTU = "%s,%s,%s-%s,%s-%s" % (tuinfos[0][8],tuinfos[0][9],tuinfos[0][3],tuinfos[0][4],tuinfos[0][6],tuinfos[0][7])
		descrNO = ";".join(["%s,%s,%s-%s,%s-%s" % (i[0][8],i[0][9],i[0][3],i[0][4],i[0][6],i[0][7]) for i in noinfos])
		## write
		line = "%s\t"*15+"%s\n"
		outgerm.write(line % (tuinfos[0][0],tuinfos[0][1],tuinfos[0][2],start1,end1,tuinfos[0][5],start2,end2,nbpab,nbpov,insert,unmap,"germline",nonum,descrTU,descrNO))
		outover.write("%s\t%s\n" % (clusname,list_germline[clusname]))

	outgerm.close()
	outover.close()
	return




def main():
	
	parser = argparse.ArgumentParser(description='First filtering of germline abnormal pairs')
	parser.add_argument("-i",action="store", type=file, required=True, dest="tuclus",  help="tab file of tumor blocs",metavar="file")
	parser.add_argument("-j",action="store", type=file, required=True, dest="noclus",  help="tab file of normal blocs",metavar="file")
	parser.add_argument("-w",action="store", type=str,  required=True, dest="dir",     help="work directory (contains overlaps files)",metavar="dir")
	parser.add_argument("-a",action="store", type=argparse.FileType('w'), required=True, dest="outtuspe",  help="output directory",metavar="dir")
	parser.add_argument("-b",action="store", type=argparse.FileType('w'), required=True, dest="outtuunm",  help="output directory",metavar="dir")
	parser.add_argument("-c",action="store", type=argparse.FileType('w'), required=True, dest="outnospe",  help="output directory",metavar="dir")
	parser.add_argument("-d",action="store", type=argparse.FileType('w'), required=True, dest="outnounm", help="output directory",metavar="dir")
	parser.add_argument("-e",action="store", type=argparse.FileType('w'), required=True, dest="outgerm",  help="output directory",metavar="dir")
	parser.add_argument("-f",action="store", type=argparse.FileType('w'), required=True, dest="outover", help="output directory",metavar="dir")
	args = parser.parse_args()

	clustu_obj = {} # clus : OVERLAPS
	clusno_obj = {} # clus : OVERLAPS
	clustu_obj = ReadInputMappable(args.dir + "/overlap_tumapp_90.bed",clustu_obj,1,False)
	clusno_obj = ReadInputMappable(args.dir + "/overlap_nomapp_90.bed",clusno_obj,1,True)
	clustu_obj = ReadInputMappable(args.dir + "/overlap_tumapp_75.bed",clustu_obj,2,False)
	clusno_obj = ReadInputMappable(args.dir + "/overlap_nomapp_75.bed",clusno_obj,2,True)
	print "input data read : \n\t",len(clustu_obj), " potential TUMOR sv candidates to remove"
	print "\t",len(clusno_obj), " potential NORMAL sv candidates to remove"

	list_tugermline = {}
	list_nogermline = set()
	list_tusomatic = {}
	clus_obj = ReadInput(args.dir + "/overlap_tu_no.bed",clustu_obj) # {clus_tu : {a/b : clus_no } }
	for tuclus,overlap in clus_obj.iteritems():
		try:
			setA = overlap.norm["a"]
			setB = overlap.norm["b"]
			intersect = setA.intersection(setB)
			if len(intersect) > 0:
				list_tugermline[tuclus] = ",".join(list(intersect))
				list_nogermline         = list_nogermline.union(intersect)
			elif len(setA) > 0 and len(setB) > 0:
				list_tusomatic[tuclus] = "AB"
			elif len(setA) > 0:
				list_tusomatic[tuclus] = "A"
			elif len(setB) > 0:
				list_tusomatic[tuclus] = "B"
		except KeyError:
			pass

	germline_tuclus = WriteOutput(list_tugermline, clustu_obj, args.tuclus, args.outtuspe,args.outtuunm,list_tusomatic)
	germline_noclus = WriteOutput(list_nogermline, clusno_obj, args.noclus, args.outnospe,args.outnounm,{})
	WriteGermlines(list_tugermline,germline_tuclus,germline_noclus,clustu_obj,clusno_obj,args.outgerm,args.outover)

	print "done."
	return 0

if __name__=="__main__":
    main()
