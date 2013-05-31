#!/usr/bin/python
import sys
from classes import PierRepeat
from Bio import SeqIO

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

	def addRepeatCopy(self, repeat, blockSize): # repeat is of type PierRepeat
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

if __name__ == '__main__':
	filename = checkAndReturnArguments()
	psl_file = open(filename)
	annotationDict = storePIERAnnotationsAsDict()
	timer = 0
	stats = RepeatStats()
	for line in psl_file:
		tabs = line.split('\t')
		seqName = tabs[9]
		seqAnnotation = annotationDict[seqName]
		currentHit = PierRepeat(seqAnnotation)
		blockSizes = convertCSVsToIntList(tabs[18]) # tabs[18] is a CSV string
		for blockSize in blockSizes:
			stats.addRepeatCopy(currentHit, blockSize)
		if timer % 400000 == 0:
			print stats
		timer += 1
	print stats
