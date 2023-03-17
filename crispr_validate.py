#!/usr/env/python

# Simple script to double-check all the primers desigend by the CRISPR pipeline and make
# sure none create tabix hits (i.e. contain an SNP with frequency > 0.01)

import io
import sys
import os
from configs import *

class Primer:

	def __init__(self, chr, num, start, end):
		self.chr = chr
		self.number = num
		self.start = start
		self.end = end

	def chr(self):
		return self.chr

	def start(self):
		return self.start

	def end(self):
		return self.end

	def number(self):
		return self.number

assayfile = open(sys.argv[1],"r")
primers = []
pnum = 1

for line in assayfile:
	data = line.split(";")
	if (data[0] == 'Chrom'):
		continue
	fp = Primer(data[0], pnum, data[4], data[5])
	rp = Primer(data[0], pnum, data[11], data[10])
	primers.append(fp)
	primers.append(rp)
	pnum += 1

for primer in primers:
	tabix_query = "tabix " + TABIX_PATH + " " + primer.chr + ":" + primer.start + "-" + primer.end
	stream = os.popen(tabix_query)
	results = stream.readlines()
	for i in results:
		result_tabs = i.split('\t')
		if (float(result_tabs[3]) <= 0.01):
			continue
		else:
			print("SNP found on primer " + str(primer.number) + "!")
			print("Primer loc: " + str(primer.start) + "-" + str(primer.end) + "  |  SNP at " + result_tabs[1] + "/" + result_tabs[2] + "\n")

