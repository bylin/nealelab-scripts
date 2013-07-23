#!/usr/bin/python
# Author: Brian Lin
# Get repeat element stats from redundant RepeatMasker .out output. Assume annotations include a mixture of Wicker annotations, Repbase annotations, and custom annotations.
# Need to sort .out file by S-W score before running!
import sys, re, glob, classes, argparse, subprocess, pickle, bitwiseSeq
from pickler import *
from Bio import SeqIO

pickledRepbaseFile = 'repbase-17.07.pkl'

def parseArgs():
	argParser = argparse.ArgumentParser(description='Get repeat element stats from RepeatMasker GFF output. Assume annotations include a mixture of Wicker annotations, Repbase annotations, and custom annotations.')
	group = argParser.add_mutually_exclusive_group(required=True)
	group.add_argument('-f', '--input_file', help='Input file')
	group.add_argument('-d', '--directory', help='Input directory')
	argParser.add_argument('-o', '--output_file', help='Output file')
	argParser.add_argument('-raw', '--raw_sequence_file', required=True, help='Raw sequence .fasta file')
	argParser.add_argument('-lib', '--library', help='Reference library used (need this flag if looking for full-length hits only)')
	argParser.add_argument('-full', '--full_length', action='store_true', help='Only calculate stats for full length hits')
	argParser.add_argument('-pickle', '--pickle_file', help='Write stats object to pickle')
	args = argParser.parse_args()
	if args.output_file:
		sys.stdout = open(args.output_file, 'w')
	return args

args = parseArgs()

def main():
	fileList = getFileList()
	stats = classes.RepeatStats()
	sys.stdout.write('Raw seqs loaded, args parsed, ready to go\n')
	for i, currentFile in zip(range(1, len(fileList)+1), fileList):
		sys.stdout.write('Parsing {}'.format(currentFile))
		addStatsFromFile(stats, currentFile)
		sys.stdout.write('{} finished, {} files remaining\n'.format(currentFile, str(len(fileList) - i)))
	if args.pickle_file:
		sendToPickleJar(stats, args.pickle_file)
	sys.stdout.write(str(stats))

def getFileList():
	if args.input_file:
		return [args.input_file]
	elif args.directory:
		fileList = glob.glob(args.directory + '/*')
		return fileList
	else: exit("Could not get file list")
	

def addStatsFromFile(stats, fileHandle):
	infile = open(fileHandle)
	# need to skip first 3 lines of GFF files
	for x in range(0,3): infile.readline()
	sortedFile = sorted(infile, reverse=True)
	i = 0
	for line in sortedFile:
		i += 1
		if i % 100 == 0:
			sys.stdout.write('.')
			sys.stdout.flush()
		try:
			repeatTuple = parseLineIntoRepeatTuple(line)
			if repeatTuple[1] == 0:
				continue
			stats += repeatTuple
		except Exception as e:
			print e
			print ("Unable to parse line: {}".format(line))

def parseLineIntoRepeatTuple(line):
	tracker = SeqTracker()
	try:
		tabs = re.findall('\S+', line)
		repeatName = tabs[9]
		rawSeqName = tabs[4]
		block = (int(tabs[5]), int(tabs[6]) + 1)
		if args.full_length:
			familyLength = tracker.libSeqs[repeatName]
			percIdentity = 100 - float(tabs[1])
	except:
		raise
	blockSize = tracker.addBlockReturnDifference(rawSeqName, block)
	if args.full_length:
		if blockSize/float(familyLength) < 0.7 or percIdentity < 80:
			blockSize = 0
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
		seqName = re.findall('\S+', seq.description)[0]
		seqs[seqName] = bitwiseSeq.BitwiseSeq(len(seq))
	return seqs

def storeLibSeqs(libSeqFile):
	seqs = {}
	for seq in SeqIO.parse(libSeqFile, 'fasta'):
		seqs[seq.id] = len(seq)
	return seqs

class SeqTracker(object):
	rawSeqs = storeRawSeqs(args.raw_sequence_file)
	if args.full_length:
		libSeqs = storeLibSeqs(args.library)
	def addBlockReturnDifference(self, name, block):
		diff = self.rawSeqs[name].addBlockReturnDifference(block)
		return diff

if __name__ == '__main__':
	main()
