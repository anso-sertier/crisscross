#!/usr/bin/env python
# -*- encoding: utf8 -*-

## CREATED 15/09/2016 by anso


import argparse

def main():
	parser = argparse.ArgumentParser(description='create fastq file from consensus sequences')
	parser.add_argument("-i",action="store", type=file, required=True,               dest="sccons",  help="summary with consensus sequences",metavar="inputfilename")
	parser.add_argument("-n",action="store", type=int, required=False, default = 2,  dest="nrmin",    help="min number of supporting reads",metavar="int")
	parser.add_argument("-l",action="store", type=int, required=False, default = 15, dest="lgmin",    help="min length of clipping",metavar="int")
	parser.add_argument("-o",action="store", type=argparse.FileType('w'), required=True, dest="fastq", help="out fastq file",metavar="outfilename")
	args = parser.parse_args()

	outline = "@%s\n%s\n+\n%s\n"
	args.sccons.next()
	for line in args.sccons:
		l = line.strip().split()
		scname = "%s:%s:%s" % (l[0],l[1], l[2].rstrip("'"))
		numCons = l[3]
		nbr     = int(l[4])
		cons    = l[6]
		lg = len(cons)
		if nbr >= args.nrmin and lg >= args.lgmin:
			readname = "%s:%s:%s:%s" % (scname,numCons,nbr,lg)
			args.fastq.write(outline % (readname, cons,"I"*lg))
	args.fastq.close()
	args.sccons.close()
	return 0


if __name__=="__main__":
    main()
