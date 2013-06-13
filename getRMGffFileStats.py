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
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument('-f', '--file', help='Input file')
	group.add_argument('-d', '--directory', help='Input directory')
	group.add_argument('input', help='Input file/directory name')
	return argParser.parse_args()
	

def addStatsFromFile(stats, fileHandle):
	myfile = open(fileHandle)
	for line in myfile:
		stats += parseLineIntoRepeatTuple(line)

def parseLineIntoRepeatTuple(line):
	tabs = line.split('\t')
	trivialParse = lambda x: return 
	# 'Target "Motif:PtRLX_772" 1 3453\n'
	(repeat, blockSize) = trivialParse