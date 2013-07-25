#!/usr/bin/python
# Author: Brian Lin
import subprocess, argparse, timestamp, re, os
from FileIO import Fasta


def parseArgs():
	parser = argparse.ArgumentParser(description="Uses probabilistic models to find likely LTR retrotransposon candidates")
	parser.add_argument('fasta', metavar='[fasta]', help="Input .fasta file, containing transposable elements")
	return parser.parse_args()

def main():
	hmmer()
	rawSeqs = {}
	for seq in Fasta(orfFile).seqRecordList():
		rawSeqs[seq.id] = ExtendedSeqRecord(seq)
	parseHits(rawSeqs)
	examineDomainsAndClassify(rawSeqs)
	os.remove(hitFile)

def hmmer():
	hmmernCmd = 'hmmern.py -hmm {} -i {} -o {}'.format(hmmProfilesForLTRs, args.fasta, hitFile)
	subprocess.call(hmmernCmd, shell=True)

def parseHits(rawSeqs):
	for line in open(hitFile):
		tabs = re.findall('\S+', line)
		if tabs[0][0] == '#': continue
		rawSeq = tabs[0]
		hmmName = tabs[3]
		orfstart = int(re.findall(':(\d+)-', tabs[24])[0])
		start = int(tabs[19]) * 3 - 2 + orfstart
		end = int(tabs[20]) * 3 - 5 + orfstart
		score = float(tabs[13])
		hmmhit = HmmHit(hmmName, start, end, score)
		try: rawSeqs[rawSeq].addHit(hmmhit)
		except KeyError:
			print '{} not in rawSeqs'.format(rawSeq)

class HmmHit(object):
	def __init__(self, name, start, end, score):
		self.start = start
		self.end = end
		self.name = name
		self.score = score
	def __str__(self):
		return '{}[{}-{}]'.format(self.name, self.start, self.end)
	def __cmp__(self, other):
		if self.start == other.start and self.end == other.end:
			return 0
		elif self.start > other.start:
			return 1
		else: return -1

class ExtendedSeqRecord(object):
	def __init__(self, seq):
		self.seq = seq
		self.hits = []

	def addHit(self, hmmhit):
		self.hits.append(hmmhit)

def examineDomainsAndClassify(rawSeqs):
	for xseq in rawSeqs:
		if len(rawSeqs[xseq].hits) == 0: continue
		string = "{}: ".format(xseq)
		string += ' '.join('{}({})'.format(hit,hit.score) for hit in sorted(rawSeqs[xseq].hits))
		#print string
		classif = classify(sorted(rawSeqs[xseq].hits), [], 0)
		print '{} : {}\n'.format(xseq, classif)

def classify(hitList, classifs, currentScore):
	for hit,i in zip(hitList, range(0,len(hitList))):
		if hit.name[0:3] == 'GAG':
			classify_GAG(hitList[i+1:], classifs, hit.score)
		elif hit.name[0:6] == 'POL-AP':
			classify_AP(hitList[i+1:], classifs, hit.score)
		elif hit.name[0:6] == 'POL-RT':
			classify_RT(hitList[i+1:], classifs, hit.score)
		elif hit.name[0:6] == 'POL-RH':
			classify_RH(hitList[i+1:], classifs, hit.score)
		elif hit.name[0:7] == 'POL-INT':
			classify_INT(hitList[i+1:], classifs, hit.score)
	return findBestClassif(classifs)

def findBestClassif(classifs):
	if len(classifs) > 0:
		classifs = sorted(classifs, key=lambda x:x[1], reverse=True)
		#print bestClassif[1]
		#print '\t' + ','.join('{}:{}'.format(c[0], c[1]) for c in classifs)
		return classifs[0][0]
	else:
		return "Unable to be categorized"

def classify_GAG(hitList, classifs, currentScore):
	for hit, i in zip(hitList, range(0, len(hitList))):
		if hit.name[0:6] == 'POL-AP':
			classify_AP(hitList[i+1:], classifs, hit.score + currentScore)
		elif hit.name[0:6] == 'POL-RT':
			classify_RT(hitList[i+1:], classifs, hit.score + currentScore)
		elif hit.name[0:6] == 'POL-RH':
			classify_RH(hitList[i+1:], classifs, hit.score + currentScore)
		elif hit.name[0:6] == 'POL-IN':
			classify_INT(hitList[i+1:], classifs, hit.score + currentScore)

def classify_AP(hitList, classifs, currentScore):
	for hit, i in zip(hitList, range(0, len(hitList))):
		if hit.name[0:6] == 'POL-RT':
			classify_RT(hitList[i+1:], classifs, hit.score + currentScore)
		elif hit.name[0:6] == 'POL-RH':
			classify_RH(hitList[i+1:], classifs, hit.score + currentScore)
		elif hit.name[0:7] == 'POL-INT':
			classify_INT(hitList[i+1:], classifs, hit.score + currentScore)

def classify_RT(hitList, classifs, currentScore):
	for hit, i in zip(hitList, range(0, len(hitList))):
		if hit.name[0:6] == 'POL-RH':
			classify_RH(hitList[i+1:], classifs, hit.score + currentScore)
		if hit.name[0:7] == 'POL-INT':
			classifs.append(('Gypsy', hit.score + currentScore))

def classify_RH(hitList, classifs, currentScore):
	for hit in hitList:
		if hit.name[0:7] == 'POL-INT':
			classifs.append(('Gypsy', hit.score + currentScore))

def classify_INT(hitList, classifs, currentScore):
	for hit, i in zip(hitList, range(0, len(hitList))):
		if hit.name[0:6] == 'POL-RT':
			classify_INT_RT(hitList[i+1:], classifs, hit.score + currentScore)
		elif hit.name[0:6] == 'POL-RH':
			classifs.append(('Copia', hit.score + currentScore))

def classify_INT_RT(hitList, classifs, currentScore):
	for hit in hitList:
		if hit.name[0:6] == 'POL-RH':
			classifs.append(('Copia', hit.score + currentScore))


if __name__ == '__main__':
	#setup global & static variables
	args = parseArgs()
	timestamp = timestamp.timestamp()
	hmmProfilesForLTRs = '/share/jyllwgrp/nealedata/databases/ltr-Pfam.hmm'
	targetSeqs = '/share/jyllwgrp/nealedata/databases/pier-1.3.fa'
	orfFile = args.fasta + '.orfs'
	hitFile = 'hmmer-output-{}.txt'.format(timestamp)
	main()
