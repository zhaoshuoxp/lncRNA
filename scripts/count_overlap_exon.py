#!/usr/bin/env pypy
#####################################
# Usage:                            #
# Manual:                           #
#####################################

import sys

trans = {}

for line in open(sys.argv[1]):
	a = line.split()
	trans.setdefault(a[0], []).append(a[1])
	
for line in open(sys.argv[2]):
	a = line.split()
	if a[0] in trans:
		trans.setdefault(a[0], []).append(a[1])
		
for i in trans:
	if int(i[0])>(int(i[1])/2.0):
		print i



################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################