#!/usr/bin/python
# Author: Brian Lin
# Loblolly pine genome annotation project: Find total BP coverage of tandems in a sequence set, given interspersed repeats and tandem repeats, with priority going to interspersed repeats. For figure ?.
import FileIO, bitwiseSeq, argparse, re
from Bio import SeqIO

def parseArgs():
	argParser = argparse.ArgumentParser(description="Gets total bp coverage of SSRs that aren't already annotated by interspersed repeats")
	argParser.add_argument('-trf', '--trf_fasta_file', required=True, help='.fa file containing TRF annotated output')
	argParser.add_argument('-rm', '--rm_output', required=True, help='RepeatMasker output for interspersed')
	argParser.add_argument('-raw', '--raw_sequence_file', required=True, help='Raw sequence .fasta file')
	return argParser.parse_args()

args = parseArgs()

def main():
	interspersed = addInterspersedToRawSeqs()
	tandems, ntandems = addTandemsToRawSeqs()
	print "Intersprsed: {}\nTandems: {}\nNumber of tandem hits: {}\n".format(interspersed, tandems, ntandems)

def addInterspersedToRawSeqs():
	infile = open(args.rm_output)
	for x in range(0,3): infile.readline()
	sortedFile = sorted(infile, reverse=True)
	total = 0
	for line in sortedFile:
		try:
			total += parseRMLineIntoBlockSize(line)
		except Exception as e:
			print e
			print ("Unable to parse line: {}".format(line))
	return total

def addTandemsToRawSeqs():
	infile = FileIO.Fasta(args.trf_fasta_file)
	total = 0
	nhits = 0
	for seq in infile:
		try:
			bp = parseTRFLineIntoBlockSize(seq.description)
			if bp == 0: continue
			total += bp
			nhits += 1
		except Exception as e:
			print e
			print ("Unable to parse seq: {}".format(seq.description))
	return total, nhits

def parseRMLineIntoBlockSize(line):
	tracker = SeqTracker()
	try:
		tabs = re.findall('\S+', line)
		rawSeqName = tabs[4]
		block = (int(tabs[5]), int(tabs[6]) + 1)
	except:
		raise
	blockSize = tracker.addBlockReturnDifference(rawSeqName, block)
	return blockSize

def parseTRFLineIntoBlockSize(line):
	tracker = SeqTracker()
	try:
		tabs = line.split('|')
		rawSeqName = tabs[0]
		block = (int(tabs[2].split(':')[0]), int(tabs[2].split(':')[1]))
	except:
		raise
	blockSize = tracker.addBlockReturnDifference(rawSeqName, block)
	return blockSize

def storeRawSeqs(rawSeqFile):
	seqs = {}
	for seq in SeqIO.parse(rawSeqFile, 'fasta'):
		seqName = re.findall('\S+', seq.description)[0]
		seqs[seqName] = bitwiseSeq.BitwiseSeq(len(seq))
	return seqs

class SeqTracker(object):
	rawSeqs = storeRawSeqs(args.raw_sequence_file)
	def addBlockReturnDifference(self, name, block):
		diff = self.rawSeqs[name].addBlockReturnDifference(block)
		return diff

if __name__ == '__main__':
	main()
