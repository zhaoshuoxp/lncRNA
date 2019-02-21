#!/usr/bin/env pypy
#####################################
# Usage:                            #
# Manual:                           #
#####################################

import sys

trans = {}

for line in open(sys.argv[1]):
	a = line.split()
	b = line.split('"')
	length = int(a[4]) - int(a[3])
	trans.setdefault(b[3],[]).append(length)

for i in trans:
	for j in trans[i]:
		length+=j
	if length >= 200:
		print(i)



################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################