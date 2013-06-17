#!/usr/bin/python
# Author: Brian Lin
# Get repeat element stats from RepeatMasker GFF output. Assume annotations include a mixture of Wicker annotations, Repbase annotations, and custom annotations.
import sys, re, classes, argparse, subprocess

def main():
	fileList = getFileList()
	stats = classes.RepeatStats()
	for currentFile in fileList:
		addStatsFromFile(stats, currentFile)
	print stats

def getFileList():
	args = parseArgs()
	if args.file:
		return [args.input]
	elif args.directory:
		fileList = subprocess.check_output('ls' + args.input).split('\n')
		return fileList[:-1]
	else: exit("Could not get file list")
	
def parseArgs():
	argParser = argparse.ArgumentParser(description='Get repeat element stats from RepeatMasker GFF output. Assume annotations include a mixture of Wicker annotations, Repbase annotations, and custom annotations.')
	group = parser.add_mutually_exclusive_group()
	group.add_argument('-f', '--file', help='Input file')
	group.add_argument('-d', '--directory', help='Input directory')
	group.add_argument('input', help='Input file/directory name')
	return argParser.parse_args()
	

def addStatsFromFile(stats, fileHandle):
	myfile = open(fileHandle)
	for line in myfile:
		stats += parseLineIntoRepeatTuple(line)

def parseLineIntoRepeatTuple(line):
	parsedLine = line.strip().split('\t')[8][14:].split(' ')
	repeatName = parsedLine[0][:-1]
	blockSize = int(parsedLine[2]) - int(parsedLine[1]) + 1
	repeat = buildRepeatFromName(repeatName)
	(repeat, blockSize) = trivialParse

# check if re.match is faster than array matching
def buildRepeatFromName(repeatName):
	if repeatName[0:2] == 'Pt': # either Wicker or custom
		if repeatName[5] == '_':
			repeat = classes.WickerRepeat(repeatName[2:5], repeatName)
		elif repeatName[7] == '_':
			repeat = classes.WickerRepeat('Uncategorized', repeatName)
#		elif repeatName[2:11] == 'Potential':
#			repeat = classes.Wi
		else:
			if repeatName in ['PtOuachita', 'PtBastrop', 'PtOzark', 'PtAppalachian', 'PtAngelina', 'PtTalladega']:
				repeat = classes.WickerRepeat('RLG', repeatName)
			elif repeatName in ['PtCumberland', 'PtPineywoods', 'PtConagree']:
				repeat = classes.WickerRepeat('RLC', repeatName)
			elif repeatName == 'PtPiedmont':
				repeat = classes.WickerREepat('RLX', repeatName)
#	else: # Repbase annotation
#		repeat = classes.PierRepeat(repeatName)
	return repeat