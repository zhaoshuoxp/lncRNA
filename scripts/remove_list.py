#!/usr/bin/env pypy
#####################################
# Usage:                            #
# Manual:                           #
#####################################
import sys

trans = []
f = open(sys.argv[3],'w')

for line in open(sys.argv[1]):
	a = line.split()
	trans.append(a[0])
	
for line in open(sys.argv[2]):
	a = line.split('"')
	if a[3] in trans:
		pass
	else:
		f.writelines(line)
	


################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################