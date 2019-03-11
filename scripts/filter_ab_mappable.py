#!/usr/bin/env python
# -*- encoding: utf8 -*-

## CREATED 17/08/2015 by anso
## MODIFIED 9/09/2016

import argparse
import re


class OVERLAPS(object):
	def __init__(self):
		self.map1 = set([])
		self.map2 = set([])

def ReadInput(fic,clus_obj,mapp):
	for line in fic:
		l = line.strip().split()
		clus = l[3][0:-1]
		num = l[3][-1]
		if not clus_obj.has_key(clus):
			clus_obj[clus] = OVERLAPS()
		if mapp == 1:
			clus_obj[clus].map1.add(num)
		else:
			clus_obj[clus].map2.add(num)
	fic.close()
	return clus_obj




def main():
	
	parser = argparse.ArgumentParser(description='First filtering of germline abnormal pairs')
	parser.add_argument("-i",action="store", type=file, required=True, dest="clus",  help="tab file of blocs", metavar="file")
	parser.add_argument("-j",action="store", type=file, required=True, dest="over90",help="overlap file (90%)",metavar="dir")
	parser.add_argument("-k",action="store", type=file, required=True, dest="over75",help="overlap file (75%)",metavar="dir")
	parser.add_argument("-o",action="store", type=argparse.FileType('w'), required=True, dest="outunmap", help="outfile name of unmappable blocs (tab)",metavar="file")
	parser.add_argument("-t",action="store", type=argparse.FileType('w'), required=True, dest="outtab", help="outfile name of mappable blocs (tab)",  metavar="file")
	parser.add_argument("-b",action="store", type=argparse.FileType('w'), required=True, dest="outbed", help="outfile name of mappable blocs (bed)",  metavar="file")
	args = parser.parse_args()

	clus_obj = {} # clus : OVERLAPS
	clus_obj = ReadInput(args.over90,clus_obj,1)
	clus_obj = ReadInput(args.over75,clus_obj,2)
	print "input data read : \n\t",len(clus_obj), " potential blocs candidates to remove"

	header  = args.clus.readline()
	args.outunmap.write(header)
	args.outtab.write(header.strip("\n") + "\tunmappable\n")
	for line in args.clus:
		l = line.split()
		clusname = l[0] + l[1]
		if clus_obj.has_key(clusname):
			unmap1 = "".join(sorted(list(clus_obj[clusname].map1)))
			unmap2 = "".join(sorted(list(clus_obj[clusname].map2)))
		else:
			unmap1   = ""
			unmap2   = "none"
		if (len(unmap1) == 2) or (len(unmap1) == 1 and len(unmap2) == 2): # either fully unmappable or 1 is 90% and other is 75%
			args.outunmap.write(line)
		else: # either fully mappable or 1 or 2 are 75%
			nbp      = int(l[8]) + int(l[9])
			args.outbed.write("%s\t%sa\t%s\n" % ("\t".join(l[2:5]),clusname,nbp))
			args.outbed.write("%s\t%sb\t%s\n" % ("\t".join(l[5:8]),clusname,nbp))
			args.outtab.write(line.strip("\n") + "\t" + unmap2 + "\n")

	args.clus.close()
	args.outbed.close()
	args.outtab.close()
	args.outunmap.close()

	print "done."
	return 0

if __name__=="__main__":
    main()
