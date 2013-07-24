#!/usr/bin/python
# Author: Brian Lin
import subprocess, argparse, timestamp, FileIO


def parseArgs():
	parser = argparse.ArgumentParser(description="Uses probabilistic models to find likely LTR retrotransposon candidates")
	parser.add_argument('fasta', metavar='[fasta]', help="Input .fasta file, containing transposable elements")
	return parser.parse_args()

def main():
	hmmer()
	parseHits()
	examineDomainsAndClassify()

def hmmer():
	hmmernCmd = 'hmmern.py -hmm {} -i {} -o {}'.format(hmmProfilesForLTRs, args.fasta, ltrFile)
	subprocess.call(hmmernCmd, shell=True)

def parseHits():

class HmmHit(object):
	def __init__(self, name, start, end):
		self.start = int(start)
		self.length = int(end) - int(start) + 1
		self

class rawTargetSeqs():
	seqs = (s for s in Fasta('/share/jyllwgrp/nealedata/databases/pier-1.3.fa').seqRecordList())
	def __init__(self):
		self.hmmhits = []
	

def examineDomainsAndClassify():



if __name__ == '__main__':
	args = parseArgs()
	timestamp = timestamp.timestamp()
	hmmProfilesForLTRs = '/share/jyllwgrp/nealedata/databases/ltr-Pfam.hmm'
	targetDatabase = '/share/jyllwgrp/nealedata/databases/pier-1.3.fa'
	ltrFile = 'LTRs-{}.txt'.format(timestamp)
	main()
