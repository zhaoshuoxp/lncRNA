#!/usr/bin/env pypy
import sys
single = open(sys.argv[4],'w')
multi = open(sys.argv[5],'w')

def get_trans(i):
	trans_list = []
	for line in open(i):
		a = line.split()
		trans_list.append(a[0])
	return trans_list
s = get_trans(sys.argv[1])
m = get_trans(sys.argv[2])

for line in open(sys.argv[3]):
	a = line.split('"')
	if a[3] in s:
		single.writelines(line)
	elif a[3] in m:
		multi.writelines(line)

