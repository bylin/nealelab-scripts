#!/usr/bin/python
from classes import PierRepeatElement

def checkAndReturnArguments():
	if (len(sys.argv) != 1):
		print "Usage: getPslFileRepeatStats.py [input psl file]"
	return sys.argv[1]

def convertCSVsToIntList(csvString):
	stringList = csvString.split(',')
	intList = []
	for value in stringList:
		intList.append(int(value))
	return intList



class RepeatStats(object):

	def __init__(self):
	self.BpsByRepeatClassifs = {}

	def addRepeatCopy(repeat, blockSize): # repeat is of type PierRepeatElement
		if repeat.CLASS in BpsByRepeatClassifs:
			BpsByRepeatClassifs[repeat.CLASS] += blockSize
		if repeat.SUPER in BpsByRepeatClassifs:
			BpsByRepeatClassifs[repeat.SUPER] += blockSize
		if repeat.ORDER in BpsByRepeatClassifs:
			BpsByRepeatClassifs[repeat.ORDER] += blockSize

	def __str__(self):
		outputString = ''
		for repeat in BpsByRepeatClassifs:
			outputString += repeat + ': ' + BpsByRepeatClassifs[repeat] + '\n'
		return outputString

def storePIERAnnotationsAsDict():
	annotationDict = {}
	for seq in SeqIO.parse('pier-1.3.fa', 'fasta'):
		(seqID, annotation) = getSeqIdAndAnnotation(seq.description)
		annotationDict[seqID] = annotation
	return annotationDict

def getSeqIdAndAnnotation(header):
	seqID = header.id
	annotation = header.description.split('\t')[1]
	return (seqID, annotation)

if __name__ == '__main__':
	filename = checkAndReturnArguments()
	annotationDict = storePIERAnnotationsAsDict()
	for line in filename:
		tabs = filename.split('\t')
		seqName = tabs[9]
		seqAnnotation = annotationDict[seqName]
		currentHit = PierRepeatElement(seqAnnotation)
		blockSizes = convertCSVsToIntList(tabs[18]) # tabs[18] is a CSV string
		stats = RepeatStats()
		for blockSize in blockSizes:
			stats.addRepeatCopy(currentHit)
	print stats