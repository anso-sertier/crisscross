#!/usr/bin/env python
# -*- encoding: utf8 -*-

## CREATED 17/08/2015 by anso
## MODIFIED 9/09/2016

import argparse
import re

class OVERLAPS(object):
	def __init__(self,tutype,noclus):
		if tutype == "a":
			self.a = set([noclus])
			self.b = set([])
		else:
			self.a = set([])
			self.b = set([noclus])
	def addNorm(self, tutype,noclus):
		if tutype == "a":
			self.a.add(noclus)
		else:
			self.b.add(noclus)

def ReadInput(fic):
	dd = {}
	for line in fic:
		l = line.strip().split()
		tuclus = l[3][0:-1]
		tutype = l[3][-1]
		noclus = l[8][0:-1]
		notype = l[8][-1]
		if tutype == notype:
			if not dd.has_key(tuclus):
				dd[tuclus] = OVERLAPS(tutype,noclus)
			else:
				dd[tuclus].addNorm(tutype, noclus)
	fic.close()
	return dd

def WriteSomatic(list_germline, dd_somatic, out, clusfic):
	header = clusfic.readline()
	out.write(header.strip("\n") + "\tgermline\n")
	dd_germ = {}
	for line in clusfic:
		l     = line.strip("\n").split()
		cname = l[0] + l[1]
		if cname in list_germline:
			dd_germ[cname] = l
			continue
		else:
			if dd_somatic.has_key(cname):
				germ = dd_somatic[cname]
			else:
				germ = "none"
			out.write("%s\t%s\n" % (line.strip("\n"),germ))
	clusfic.close()
	out.close()
	return dd_germ

def WriteGermline(dd_germline, tuclus, noclus, outgerm, outover):
	outgerm.write("type\tnumTU\tchrom1\tstart1\tend1\tchrom2\tstart2\tend2\tnbpab\tnbpov\tinsert\tunmappable\tnumNO\tdescrTU\tdescrNO\n")
	for clusname in dd_germline:
		nonames = dd_germline[clusname].split(",")
		noinfos = [noclus[i] for i in nonames]
		tuinfos = tuclus[clusname]
		
		start1 = min(int(tuinfos[3]), min([int(i[3]) for i in noinfos]))
		end1   = max(int(tuinfos[4]), max([int(i[4]) for i in noinfos]))
		start2 = min(int(tuinfos[6]), min([int(i[6]) for i in noinfos]))
		end2   = max(int(tuinfos[7]), max([int(i[7]) for i in noinfos]))
		nbpab  =    int(tuinfos[8]) + sum([int(i[8]) for i in noinfos])
		nbpov  =    int(tuinfos[9]) + sum([int(i[9]) for i in noinfos])
		if clusname.startswith("inter"):
			insert = 0
		else:
			insert = start2 - end1

		unmap   = "".join(list( set(list(tuinfos[11])).union([j for i in noinfos for j in list(i[11])])))
		nonum   = ",".join([i[1] for i in noinfos])
		descrTU = "%s,%s,%s-%s,%s-%s" % (tuinfos[8],tuinfos[9],tuinfos[3],tuinfos[4],tuinfos[6],tuinfos[7])
		descrNO = ";".join(["%s,%s,%s-%s,%s-%s" % (i[8],i[9],i[3],i[4],i[6],i[7]) for i in noinfos])
		b1      = "\t".join(tuinfos[0:3])

		## write
		line = "%s\t"*12+"%s\n"
		outgerm.write(line % (b1,start1,end1,tuinfos[5],start2,end2,nbpab,nbpov,insert,unmap,nonum,descrTU,descrNO))
		outover.write("%s\t%s\n" % (clusname,dd_germline[clusname]))

	outgerm.close()
	outover.close()
	return
	


def main():
	
	parser = argparse.ArgumentParser(description='First filtering of germline abnormal pairs')
	parser.add_argument("-i",action="store", type=file, required=True, dest="tuclus",  help="tab file of mappable tumor blocs",metavar="file")
	parser.add_argument("-j",action="store", type=file, required=True, dest="noclus",  help="tab file of mappable normal blocs",metavar="file")
	parser.add_argument("-k",action="store", type=file, required=True, dest="overlap", help="bed file of overlap blocs",metavar="file")
	parser.add_argument("-o",action="store", type=argparse.FileType('w'), required=True, dest="outtuspe", help="outfile of specific tumor blocs",metavar="file")
	parser.add_argument("-p",action="store", type=argparse.FileType('w'), required=True, dest="outnospe", help="outfile of specific normal blocs",metavar="file")
	parser.add_argument("-q",action="store", type=argparse.FileType('w'), required=True, dest="outgerm",  help="outfile of germline blocs",metavar="file")
	parser.add_argument("-r",action="store", type=argparse.FileType('w'), required=True, dest="outover",  help="outfile of overlaps",metavar="file")
	args = parser.parse_args()

	overlaps = ReadInput(args.overlap) # {clus_tu : {a/b : clus_no } }
	print "overlaps read"
	
	tugermline = {}	   # {tu_clusname : "noname1,noname2"}
	nogermline = set() # (noname1, noname2,...)
	tusomatic  = {}    # {tu_clusname : "AB"}
	i = 0
	for tuclusname,over in overlaps.iteritems():
		setA = over.a
		setB = over.b
		noclusnames = setA.intersection(setB)
		if len(noclusnames) > 0:
			i += 1
			tugermline[tuclusname] = ",".join(list(noclusnames))
			nogermline             = nogermline.union(noclusnames)
		elif len(setA) > 0 and len(setB) > 0:
			tusomatic[tuclusname] = "AB"
		elif len(setA) > 0:
			tusomatic[tuclusname] = "A"
		elif len(setB) > 0:
			tusomatic[tuclusname] = "B"
	print i, " intersection not null"
	print "treatment done, write now"	
	tugermblocs = WriteSomatic(set(tugermline.keys()),tusomatic,args.outtuspe, args.tuclus)
	nogermblocs = WriteSomatic(nogermline, {},       args.outnospe, args.noclus)
	WriteGermline(tugermline, tugermblocs, nogermblocs,   args.outgerm, args.outover)

	print "done."
	return 0

if __name__=="__main__":
    main()
