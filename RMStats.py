#!/usr/bin/python
# Author: Brian Lin
# Get repeat element stats from redundant RepeatMasker .out output. Assume annotations include a mixture of Wicker annotations, Repbase annotations, and custom annotations.
# Need to sort .out file by S-W score before running!
import sys, re, glob, classes, argparse, subprocess, pickle, bitwiseSeq
from pickler import getFromPickleJar
from Bio import SeqIO

pickledRepbaseFile = 'repBaseDict.pkl'
rawSeqFile = 'ptaeda.v0.9_bin001_of100.fsa'

def main():
	args = parseArgs()
	fileList = getFileList(args)
	stats = classes.RepeatStats()
	for i, currentFile in zip(range(1, len(fileList)+1), fileList):
		sys.stderr.write('Parsing {} ...\n'.format(currentFile))
		addStatsFromFile(stats, currentFile)
		sys.stderr.write('{} finished, {} files remaining\n'.format(currentFile, str(len(fileList) - i)))
	sys.stderr.write(str(stats))

def getFileList(args):
	if args.file:
		return [args.file]
	elif args.directory:
		fileList = glob.glob(args.directory + '/*')
		return fileList
	else: exit("Could not get file list")
	
def parseArgs():
	argParser = argparse.ArgumentParser(description='Get repeat element stats from RepeatMasker GFF output. Assume annotations include a mixture of Wicker annotations, Repbase annotations, and custom annotations.')
	group = argParser.add_mutually_exclusive_group(required=True)
	group.add_argument('-f', '--file', help='Input file')
	group.add_argument('-d', '--directory', help='Input directory')
	argParser.add_argument('-o', '--output_file', help='Output file')
	argParser.add_argument('-e', '--error_file', help='Error file, status messages will also be posted in this file')
	args = argParser.parse_args()
	if args.output_file:
		sys.stdout = open(args.output_file, 'w')
	if args.error_file:
		sys.stderr = open(args.error_file, 'w')
	return args
	

def addStatsFromFile(stats, fileHandle):
	infile = open(fileHandle)
	# need to skip first 3 lines of GFF files
	for x in range(0,3): infile.readline()
	sortedFile = sorted(infile, reverse=True)
	counter = 0
	for line in sortedFile:
		counter += 1
		if x % 100 == 0:
			sys.stdout.write('.')
			sys.stdout.flush()
		try:
			repeatTuple = parseLineIntoRepeatTuple(line)
			print '{} {}'.format(repeatTuple[0], repeatTuple[1])
			stats += repeatTuple
		except:
			sys.stderr.write("Unable to parse line: {}".format(line))

def parseLineIntoRepeatTuple(line):
	try:
		tabs = re.findall('\S+', line)
		repeatName = tabs[9]
		rawSeqName = tabs[4]
		block = (int(tabs[5]), int(tabs[6]) + 1)
	except:
		raise
	tracker = RawSeqTracker()
	blockSize = tracker.addBlockReturnDifference(rawSeqName, block)
	repeat = buildRepeatFromName(repeatName)
	return (repeat, blockSize)

# check if re.match is faster than array matching
def buildRepeatFromName(repeatName):
	# if not a Wicker annotation: use PierRepeat().
	isWicker = repeatName[0:2] == 'Pt' and (repeatName[5] == '_' or (len(repeatName) > 5 and (repeatName[2:6] in ['NoCa', 'Pote', 'RDNA'] or repeatName[2:5] == 'SSR')))
	matcher = RepbaseMatcher()
	if isWicker:
		code = repeatName[2:].split('_')[0]
		repeat = classes.WickerRepeat(code, repeatName)
	elif matcher.isRepbaseRepeat(repeatName):
		rawRepeat = matcher.match(repeatName)
		repeat = classes.conversion(rawRepeat)
	else:
		repeat = classes.PierRepeat(repeatName)
	return repeat

class RepbaseMatcher(object):
	repbase = getFromPickleJar(pickledRepbaseFile) # dict containing repBaseRepeat objects
	
	def isRepbaseRepeat(self, repeatName):
		if repeatName in self.repbase:
			return True
		else:
			return False

	def match(self, repeatName):
		return self.repbase[repeatName]

def storeRawSeqs(rawSeqFile):
	seqs = {}
	for seq in SeqIO.parse(rawSeqFile, 'fasta'):
		seqs[seq.description] = bitwiseSeq.BitwiseSeq(len(seq))
	return seqs

class RawSeqTracker(object):
	rawSeqs = storeRawSeqs(rawSeqFile)
	def addBlockReturnDifference(self, name, block):
		diff = self.rawSeqs[name].addBlockReturnDifference(block)
		print name, block, diff
		return diff

if __name__ == '__main__':
	sys.stderr.write('Raw sequences loaded!\n')
	main()
