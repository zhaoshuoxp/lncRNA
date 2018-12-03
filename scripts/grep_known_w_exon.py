#!/usr/bin/env pypy
import sys
trans_list = []
f = open(sys.argv[3],'w')

for line in open(sys.argv[1]):
	a = line.split()
	trans_list.append(a[0])

for line in open(sys.argv[2]):
	a = line.split('"')
	for i in a:
		if i in trans_list:
			f.writelines(line)

