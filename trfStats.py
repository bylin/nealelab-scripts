#!/usr/bin/python
#Author: Jacob Zieve
import sys
from pickler import sendToPickleJar, getFromPickleJar
from classes import trfHit
from collections import defaultdict
from operator import itemgetter, attrgetter

"""Prints non-overlapping and cleaned version of trf output for load into database"""
"""Sample line:
33493 33544 22 2.3 22 83 3 59 34 3 1 59 1.27 ATTTTTAAAACATTTTTAATTC ATTTTAAAAATATTTTTAATTTATTTTTAAAACATTTTTTATTGCATTTTTA
"""
"""Handling overlap like greedy activity scheduling"""


def parse_dat(dat_file):
	hitsDict = defaultdict(list)
	seq_name = ""
	for line in open(dat_file).read().split("\n"):	
		if line.startswith("Seq"):
			seq_name = line[len("Sequence: "):]
		if line != "" and line.split()[0].isdigit():
			fields = line.split()
			fields.insert(0,seq_name)
			hitsDict[seq_name].append((trfHit(fields)))

	return hitsDict

def greedyOverlapFiltering(hitsDict):
	tandems = []
	for seq_name in hitsDict.keys():
		seqs=sorted(hitsDict[seq_name],key=attrgetter('end'), reverse=True)
		tandems.append(seqs[0])
		f = seqs[0].end
		for i in range(1,len(seqs)):
			if seqs[i].start >= f:
				tandems.append(seqs[i])
				f = seqs[i].end
	return tandems

def getPeriodFreqs(tandems):
	periodDict = defaultdict(int)
	for t in tandems:
		periodDict[t.period]+=1
	return periodDict

def printPeriodFreqs(tandems):
	for k,v in sorted(getPeriodFreqs(tandems).items()):
		print "{:d},{:d}".format(k,v)

def getCumulativeBpByPeriodRange(tandems,prange):
	return sum([t.length for t in filter(lambda x: x.period in range(prange[0],prange[1]),tandems)])

def printCumulativeBpByPeriodRange(tandems,prange):
	print "Between {:d} and {:d}:\t{:d}".format(prange[0],prange[1],getCumulativeBpByPeriodRange(tandems,(prange[0],prange[1])))

def main():
#	tandems = greedyOverlapFiltering(parse_dat(sys.argv[1]))
#	sendToPickleJar(tandems,'Ptaeda_trf.txt')
#	sendToPickleJar(tandems,'Pglauca_trf.txt')
#	sendToPickleJar(tandems,'Pabies_trf.txt')
#	tandems = getFromPickleJar('Ptaeda_trf.txt')
	tandems = getFromPickleJar('Pglauca_trf.txt')
#	tandems = getFromPickleJar('Pabies_trf.txt')
#	printPeriodFreqs(tandems)	
	printCumulativeBpByPeriodRange(tandems,(2,3))
	printCumulativeBpByPeriodRange(tandems,(21,22))
	printCumulativeBpByPeriodRange(tandems,(102,103))
	
	#printCumulativeBpByPeriodRange(tandems,(2,9))
#	printCumulativeBpByPeriodRange(tandems,(10,101))
#	printCumulativeBpByPeriodRange(tandems,(102,501))
		
if __name__=="__main__":
	main()

