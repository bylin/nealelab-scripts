#!/usr/bin/python
import sys
from classes import PierRepeat
from Bio import SeqIO
import TupleMergeGenerator

def checkAndReturnArguments():
	if (len(sys.argv) != 2):
		print "Usage: getPslFileRepeatStats.py [input psl file]"
	return sys.argv[1]

def convertCSVsToIntList(csvString):
	stringList = csvString.split(',')
	intList = []
	for value in stringList[:-1]: # the last value is an empty string; raw string has a trailing comma
		intList.append(int(value))
	return intList

class RepeatStats(object):

	def __init__(self):
		self.BpsByRepeatClassifs = {}

	def addRepeatCopy(self, repeat, repeatLength): # repeat is of type PierRepeat
		if repeat.CLASS in self.BpsByRepeatClassifs:
			self.BpsByRepeatClassifs[repeat.CLASS] += blockSize
		else: 
			self.BpsByRepeatClassifs[repeat.CLASS] = 0
		if repeat.SUPER in self.BpsByRepeatClassifs:
			self.BpsByRepeatClassifs[repeat.SUPER] += blockSize
		else: 
			self.BpsByRepeatClassifs[repeat.SUPER] = 0
		if repeat.ORDER in self.BpsByRepeatClassifs:
			self.BpsByRepeatClassifs[repeat.ORDER] += blockSize
		else: 
			self.BpsByRepeatClassifs[repeat.ORDER] = 0

	def __str__(self):
		outputString = ''
		for repeat in self.BpsByRepeatClassifs:
			repeatBps = str(self.BpsByRepeatClassifs[repeat])
			outputString += repeat + ': ' + repeatBps + ', '
		return outputString[:-2] # truncate the trailing ', '

def storePIERAnnotationsAsDict():
	annotationDict = {}
	for seq in SeqIO.parse('/share/jyllwgrp/nealedata/databases/pier-1.3.fa', 'fasta'):
		(seqID, annotation) = getSeqIdAndAnnotation(seq.description)
		annotationDict[seqID] = annotation
	return annotationDict

def getSeqIdAndAnnotation(header):
	seqID = header.split('\t')[0]
	annotation = header.split('\t')[1]
	return (seqID, annotation)

class NonredundantBlockStorage(object):
	def __init__(self):
		self.targetBlocks = []

	def addAndFilterBlock(self, seqID, newBlock):
		if seqID not in self.targetSeqs:
			newTargetSeq = TargetSequenceSubseqs(seqID)
			targetSeqs.append(newTargetSeq)
		self.targetBlocks[seqID].mergeNewSeq(newBlock)


class TargetSequenceSubseqs(object):
	def __init__(self, seqID):
		self.seqID = seqID
		self.subseqs = []
	
	def __eq__(self, other):
		return self.seqID == other.seqID

	def mergeNewSeq(self, newSubseq):
		newSubseqList = self.subseqs.append(newSubseq)
		self.subseqs = list(TupleMergeGenerator.merge(newSubseqList))

if __name__ == '__main__':
	filename = checkAndReturnArguments()
	psl_file = open(filename)
	annotationDict = storePIERAnnotationsAsDict()
	timer = 0
	stats = RepeatStats()
	NonredundantBlocks = NonredundantBlockStorage()
	for line in psl_file:
		tabs = line.split('\t')
		seqName = tabs[9]
		targetSeqName = tabs[12]
		seqAnnotation = annotationDict[seqName]
		currentHit = PierRepeat(seqAnnotation)

		#TODO: convert to generator
		for block in convertCSVsToIntList(tabs[18], tabs[19], tabs[20]) # tabs[18] is a CSV string
			filteredBlockSize = NonredundantBlocks.addAndFilterBlock(targetSeqName, block)
			stats.addRepeatCopy(currentHit, blockSize)
		if timer % 400000 == 0:
			print stats
		timer += 1
	print stats
