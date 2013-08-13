#!/usr/bin/python
# Author: Brian Lin

import sys
import re
import argparse
import errno

argParser = argparse.ArgumentParser()
argParser.add_argument("inputFile", help="Variable length space file to be turned into a tab delimited file")
args = argParser.parse_args()

handle = open(args.inputFile)
for line in handle:
	tabs = re.findall('\S+', line)
	try:
		print '\t'.join(tabs)
		sys.stdout.flush()
	except IOError as e:
		if e.errno == errno.EPIPE: # error 32, broken pipe
			exit()

