#!/usr/bin/env pypy
import sys
exon = {}
single = open(sys.argv[2],'w')
multi = open(sys.argv[3],'w')

for line in open(sys.argv[1]):
	a = line.split('"')
	exon.setdefault(a[3],[]).append(a[5])

for i in exon:
	if len(exon[i])>=2:
		multi.writelines(i+'\n')
	else:
		single.writelines(i+'\n')

