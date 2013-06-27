#!/usr/bin/python
# Author: Brian Lin
# Get repeat element stats from redundant RepeatMasker .out output. Assume annotations include a mixture of Wicker annotations, Repbase annotations, and custom annotations.
# Need to sort .out file by S-W score before running!
import sys, re, glob, classes, argparse, subprocess, pickle, bitwiseSeq
from pickler import getFromPickleJar

pickledRepbaseFile = 'repBaseDict.pkl'
pickledRawSequenceFile = 'bitwiseRawSeqs.pkl'

def main():
	fileList = getFileList()
	stats = classes.RepeatStats()
	for i, currentFile in zip(range(1, len(fileList)+1), fileList):
		print 'Parsing ' + currentFile + '...'
		addStatsFromFile(stats, currentFile)
		print currentFile + ' finished, ' + str(len(fileList) - i) + ' files remaining'
	print stats

def getFileList():
	args = parseArgs()
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
	return argParser.parse_args()
	

def addStatsFromFile(stats, fileHandle):
	gffFile = open(fileHandle)
	# need to skip first 3 lines of GFF files
	for x in range(0,3): gffFile.readline()
	for line in gffFile:
		stats += parseLineIntoRepeatTuple(line)

def parseLineIntoRepeatTuple(line):
	tabs = re.findall('\W+', line)
	repeatName = tabs[9]
	rawSeqName = tabs[4]
	block = (int(tabs[5]), int(tabs[6]))
	tracker = RawSeqTracker()
	difference = tracker.addBlockReturnDifference(rawSeqName, block)
	blockSize = block[1] - block[0] - difference
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

class RawSeqTracker(object):
	rawSeqs = getFromPickleJar(pickledRawSequenceFile) # dict containing bitwise seqs
	
	def addBlockReturnDifference(self, name, block):
		rawSeqs[name] += block
		return rawSeqs[name].findBlockDifference(block)

if __name__ == '__main__':
	main()
