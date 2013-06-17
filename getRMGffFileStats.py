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
		return [args.file]
	elif args.directory:
		fileList = subprocess.check_output('ls' + args.directory).split('\n')
		return fileList[:-1]
	else: exit("Could not get file list")
	
def parseArgs():
	argParser = argparse.ArgumentParser(description='Get repeat element stats from RepeatMasker GFF output. Assume annotations include a mixture of Wicker annotations, Repbase annotations, and custom annotations.')
	group = argParser.add_mutually_exclusive_group(required=True)
	group.add_argument('-f', '--file', help='Input file')
	group.add_argument('-d', '--directory', help='Input directory')
	return argParser.parse_args()
	

def addStatsFromFile(stats, fileHandle):
	myfile = open(fileHandle)
	# need to skip first 3 lines of GFF files
	for x in range(0,3): myfile.readline()
	for line in myfile:
		stats += parseLineIntoRepeatTuple(line)

def parseLineIntoRepeatTuple(line):
	parsedLine = line.strip().split('\t')[8][14:].split(' ')
	repeatName = parsedLine[0][:-1]
	blockSize = int(parsedLine[2]) - int(parsedLine[1]) + 1
	repeat = buildRepeatFromName(repeatName)
	return (repeat, blockSize)

# check if re.match is faster than array matching
def buildRepeatFromName(repeatName):
	if repeatName[0:2] == 'Pt': # either REPET/Wicker, or custom
		code = repeatName[2:].split('_')[0]
			repeat = classes.WickerRepeat(code)
	else: # Repbase annotation
		repeat = classes.PierRepeat(repeatName)
	print repeatName
	return repeat

if __name__ == '__main__':
	main()
