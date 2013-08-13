#!/usr/bin/python
# Author: Brian Lin

import re, argparse

argParser = argparse.ArgumentParser()
argParser.add_argument("input", help="Variable length space file to be turned into a tab delimited file")
args = argParser.parse_args()

for line in args.input:
	tabs = re.findall('\S+', line)
	print '\t'.join(tabs)

