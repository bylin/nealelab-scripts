import re
import sys
import pickler
if len(sys.argv) < 3:
	print 'usage: annotReplacer.py [infile] [outfile]'
infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')
newClassifs = pickler.getFromPickleJar('newclassifs')
i=0
for line in infile:
	i += 1
	if i <= 3:
		outfile.write(line)
		continue
	tabs = re.findall('\S+', line)
	name = tabs[9]
	family = '_'.join(name.split('\t')[:2])
	if family in newClassifs:
		line = line.replace(family, newClassifs[family])
	outfile.write(line)
	
