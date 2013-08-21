#!/usr/bin/python
import re
import str
import sys
a = open(sys.argv[1], 'r')

class rep(object):
	def __init__(self, n, bp):
		self.nCopies = n
		self.length = bp
	def __str__(self):
		return '{}\t{}'.format(self.nCopies, self.length)
	def __gt__(self, other):
		return self.length > other.length

def parse(a, classif):
	nextline = a.readline()
	tabs = nextline.split('\t')
	while len(tabs) != 1:
		#sys.stdout.write('.')
		#sys.stdout.flush()
		if tabs[2] == '0':
			nextline = a.readline()
			continue
		elif tabs[0] in classif:
			classif[tabs[0]].nCopies += int(tabs[1])
			classif[tabs[0]].length += int(tabs[2])
		else:
			classif[tabs[0]] = rep(int(tabs[1]), int(tabs[2]))
		nextline = a.readline()
		tabs = nextline.split('\t')
	return nextline


classes = {}
orders = {}
supers = {}
families = {}
total = rep(0, 0)
linecount = 0
while True:
	line = a.readline()
	if not line: break
	if line == '===classes===\n':
		line = parse(a, classes)
	if line == '===orders===\n':
		line = parse(a, orders)
	if line == '===supers===\n':
		line = parse(a, supers)
	if line == '===families===\n':
		line = parse(a, families)
	if line[0:6] == 'TOTAL\t':
		tabs = line.split('\t')
		total.nCopies += int(tabs[1])
		total.length += int(tabs[2])
print 'TOTAL\t{}'.format(total)
for classif, dicts in zip(['classes', 'orders', 'supers', 'families'], [classes, orders, supers, families]):
	print str.upper(classif)
	current = sorted(dicts, key=lambda key:dicts[key], reverse=True)
	for classif in current:
		print '{}\t{}'.format(classif, dicts[classif])

