#!/usr/bin/python

###################################################
# Write to TSVs reads overlapping 7076, noting A/G
###################################################

import sys
import re
import os
import pysam
import numpy as np
from collections import defaultdict

bamfile = sys.argv[1]
outpre = bamfile
maxBP = 16569
base_qual = 0

# Export Functions
def writeSparseMatrix(mid, vec):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(0,int(maxBP)):
			if(vec[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec[i])+"\n")

n = int(maxBP)
alignment_quality = 0

# initialize with a pseudo count to avoid dividing by zero
counts_A = [0] * n 
counts_G = [0] * n 

# organize reads into a dict where key is readname
bam2 = [x for x in pysam.AlignmentFile(bamfile, "rb")]
ordered_bam2 = defaultdict(list)
for read in bam2:
	ordered_bam2[read.query_name].append(read)

bam2 = pysam.AlignmentFile(bamfile, "rb")
for read in bam2:
	seq = read.seq
	reverse = read.is_reverse
	quality = read.query_qualities
	align_qual_read = read.mapping_quality
	for qpos, refpos in read.get_aligned_pairs(True):
		if qpos is not None and refpos is not None and align_qual_read > alignment_quality:
			if(refpos + 1 == 7076 and (seq[qpos] == "G" or seq[qpos] == "A")):
				print(read.reference_start, read.reference_length, seq[qpos])
			exit

