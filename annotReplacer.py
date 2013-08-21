#!/usr/bin/python
import re
import sys
import pickler
if len(sys.argv) < 3:
	print 'usage: annotReplacer.py [infile] [outfile]'
	exit()

infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')
newClassifs = pickler.getFromPickleJar('newclassifs')
for line in infile:
	tabs = re.findall('\S+', line)
	try:
		name = tabs[9]
		family = '_'.join(name.split('_')[:2])
		classif = name.split('_')[0]
		if family in newClassifs:
			#print family, newClassifs[family]
			line = line.replace(family, newClassifs[family])
		elif classif in newClassifs:
			line = line.replace(classif, newClassifs[classif])
	except:
		print '1'
		if name in newClassifs:
			line = line.replace(name, newClassifs[name])
	outfile.write(line)
	
