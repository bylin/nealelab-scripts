#!/bin/usr/python
#Author: Jacob Zieve, Brian Lin

"""Defines different "types" of data inputs that are handled
int the stats scritps. Many of the classes include overloaded equality, hashing
and comparison operators to be used in the script"""

import re
from collections import defaultdict
from Bio import SeqIO
#used as main data structure for sequence
	
class trfHit():
	def __init__(self, Data=[]):
		self.seq_name = Data[0]
		self.start = int(Data[1]) #Intrices of the repeat relative to the start of the sequence 
		self.end = int(Data[2])
		self.period = int(Data[3])#Period size of the repeat
		self.copy_num = float(Data[4])#Number of copies aligned with the consensus pattern
		self.size_of_copy = int(Data[5])#Size of consensus pattern (may differ slightly from the period size)
		self.percent_matches = float(Data[6])#Percent of matches between adjacent copies overall
		self.percent_indels = float(Data[7])#Percent of indels between adjacent copies overall
		self.score = int(Data[8])#Alignment score
		self.a = float(Data[9])#Percent composition for each of the four nucleotides
		self.c = float(Data[10])
		self.g = float(Data[11])
		self.t = float(Data[12])
		self.entropy = float(Data[13])#Entropy measure based on percent composition
		self.sequence = Data[14]#The repeat motif itself
		self.length = self.end - self.start#Alignment length
					
	def __str__(self): 
		#for printing in mysql format
		return "'" + str(self.seq_name) + "'," + str(self.start) + "," + str(self.end) + ","	+ str(self.period) + "," + str(self.copy_num)+ ",'" + str(self.sequence)+ "'," + str(self.length)
	
	def __cmp__(self,other):
		return self.score > other.score #blast score, metric used to compare redundant bp
	def __eq__(self,other):
		return self.seq_name == other.seq_name

	def __hash__(self):
		return hash(self.sequence)

class repBaseRepeat:
	def __init__(self,name, accession, family, species, taxonomy, length):
		self.name = name
		self.accession = accession
		self.family = family
		self.species = species
		self.taxonomy = taxonomy
		self.length = length
	def __str__(self):
		return "name: "+self.name+"\nac: "+self.accession+"\nfam: "+self.family+"\n\
		sp: "+self.species+"\ntax: "+self.taxonomy+"len: "\
		+self.length+"\n"

class censorRepeat(object):# (object) allows for type checking, will return censorRepeat rather than instance
	def __init__(self,Data = []):
		self.seq = Data[0]
		self.sstart = int(Data[1]) #Indices of the repeat relative to the start of the sequence. 
		self.send = int(Data[2])
		self.name = Data[3]
		self.rstart = int(Data[4])
		self.rend = int(Data[5])
		self.orientation = Data[6]
		self.sim = float(Data[7])
		self.pos = float(Data[8])
		self.score = int(Data[9])
		self.per_query = float(Data[10])
		self.per_repeat = float(Data[11])	
		self.length = self.rend-self.rstart #length of alignment is defined from indices in the repeat
		self.isDnvo = False

	def __cmp__(self,other):
		return self.score > other.score
		
	def __eq__(self,other):
		return self.name == other.name
	
	def __hash__(self):
		return hash(self.name)
		
	def __str__(self):
		return 'name: {:s}\nscore: {:d}\n'.format(self.name,self.score)

class blastHit:
	def __init__(self, data = []):
		self.qid = data[0]#query id
		self.sid = data[1]#subject/database id
		self.pid = float(data[2])#percent identity
		self.alen = int(data[3])#alignment length
		self.m = int(data[4])#mismatches
		self.g = int(data[5])#gaps
		self.qstart = int(data[6])#indice start of query
		self.qend = int(data[7])#indice end of query
		self.sstart = int(data[8])#" " of subject
		self.send = int(data[9])#" " of subject
		self.ev = float(data[10])#e-value
		self.bs = float(data[11])#bit score
	def __str__(self):
		return '{:s} {:s} {:f} {:d} {:d} {:d} {:d}\
		 {:d} {:d} {:d} {:f} {:f}\n'.format(self.qid,\
		 self.sid, self.pid, self.alen, self.m, self.g, self.qstart,\
		 self.qend, self.sstart, self.send, self.ev, self.bs)

class Repeat(object):

	def __init__(self, code, name):
		self.NAME = name
		self.CLASS = self.determineClass(code)
		self.ORDER = self.determineOrder(code)
		self.SUPER = self.determineSuper(code)
		self.FAMILY = self.determineFamily(name)
	
	def __str__(self):#for debug
		return ', '.join([self.NAME, self.CLASS, self.ORDER, self.SUPER, self.FAMILY])

	def __eq__(self, other):
		return self.name == other.name

	def __hash__(self): #used in set() function later
		return hash(self.name)

class WickerRepeat(Repeat):
	# Example: myRepeat = WickerRepeat('RLX')

	def determineClass(self, code):
		myClass = code[0]
		if myClass == 'R': return 'I'
		elif myClass == 'D': return 'II'
		elif myClass == 'X': return 'Confused'
		elif code in ['PotentialHostGene', 'NoCat', 'rRNA', 'SSR']:
			return code
	
	def determineOrder(self, code):
		myClass = self.CLASS
		order = code[1]
		if myClass == 'I':
			if order == 'X': return 'Unknown'
			elif order == 'L': return 'LTR'
			elif order == 'I': return 'LINE'
			elif order == 'Y': return 'DIRS'
			elif order == 'S': return 'SINE'
			elif order == 'P': return 'Penelope'
		elif myClass == 'II':
			if order == 'X': return 'Unknown'
			elif order == 'T': return 'TIR'
			elif order == 'Y': return 'Crypton'
			elif order == 'H': return 'Helitron'
			elif order == 'M': return 'Maverick'
		elif myClass in ['Confused', 'PotentialHostGene', 'NoCat', 'rRNA', 'SSR']:
			return 'Unknown'
		else:
			print self
			exit("Unhandled case: " + code)

		
	def determineSuper(self, code):
		super = code[2]
		if super == 'G': return 'Gypsy'
		elif super == 'C': return 'Copia'
		else: return 'Unknown'

	def determineFamily(self, name):
		if re.search('_[A-Z]', name):
			return '_'.join(name.split('_')[:-1])
		return name

# TODO: flesh out class, order, super functions. Repbase can be tricky.
class PierRepeat(Repeat):
	# Example: myRepeat = WickerRepeat('PtPiedmont')
	def __init__(self, name):
		self.NAME = name
		self.FAMILY = self.determineFamily(name)
		self.CLASS = self.determineClass(self.FAMILY)
		self.ORDER = self.determineOrder(self.FAMILY)
		self.SUPER = self.determineSuper(self.FAMILY)
		
	def determineFamily(self, name): #truncate subfamily
		if re.search('_[A-Z]', name) and name[0:2] == 'Pt':
			return '_'.join(name.replace('-', '_').split('_')[:-1])
		return name

	classIRE = re.compile('ifg7|tpe1|gymny|copia|corky|ouachita|piedmont|bastrop|ozark|appalachian|angelina|talladega|cumberland|pineywoods|conagree|gypsy|ltr|class I|penelope|line|sine|non-ltr|retrotransposon|endogenous retrovirus|rte|alisei|shacop|rn107|ty1_pe|piln|sc-10', re.I)
	classIIRE = re.compile('II|helitron|dirs|mudr|maverick|enspm|harbinger|helmet', re.I)
	ltrRE = re.compile('ifg7|ouachita|bastrop|piedmont|ozark|appalachian|angelina|talladega|cumberland|pineywoods|conagree|gypsy|copia|tpe1|corky|gymny|ltr|endogenous retrovirus|alisei|shacop|rn107|sc-10', re.I)
	gypsyRE = re.compile('ifg7|gymny|gypsy|ouachita|piedmont|bastrop|ozark|appalachian|angelina|talladega|corky|alisei', re.I)
	copiaRE = re.compile('cumberland|pineywoods|conagree|copia|tpe1|shacop|rn107|sc-10', re.I)


	def determineClass(self, annotation):
		if PierRepeat.classIRE.search(annotation) or annotation in ['L1', 'I', 'Sm1'] or re.search('COP', annotation):
			return 'I'
		elif PierRepeat.classIIRE.search(annotation) or annotation in ['II', 'DIR', 'P', 'TIR'] or re.search('hAT', annotation):
			return 'II'
		elif annotation in ['CENSATC4_ZM', 'TREP60', 'ATREP18']:
			return 'SSR'
		elif annotation == 'PtRDNA1':
			return 'rDNA'
		else:
			return 'Unknown'

	def determineOrder(self, annotation):
		if self.CLASS == 'Unknown':
			print '!!! Unknown! ' + annotation
		if PierRepeat.ltrRE.search(annotation) or re.search('COP', annotation):
			return 'LTR'
		elif annotation in ['LINE', 'L1', 'RTE', 'VLINE6_VV', 'PILN1_PT']:
			return 'LINE'
		elif annotation in ['SINE', 'SINE1/7SL', 'SINE2/tRNA']:
			return 'SINE'
		elif annotation in ['P', 'Mariner/Tc1', 'Harbinger', 'TIR'] or re.search('mudr|enspm', annotation, re.I) or re.search('hAT', annotation):
			return 'TIR'
		elif annotation in ['DIR', 'DIRS']:
			return 'DIRS'
		elif annotation in ['Penelope', 'Sm1']:
			return 'Penelope'
		elif annotation in ['Helitron', 'HELMET3']:
			return 'Helitron'
		elif annotation in ['SINE', 'Maverick']:
			return annotation
		else:
			return 'Unknown'

	def determineSuper(self, annotation):
		if PierRepeat.gypsyRE.search(annotation):
			return 'Gypsy'
		elif PierRepeat.copiaRE.search(annotation) or re.search("COP", annotation, re.I):
			return 'Copia'
		elif annotation in ['DIRS', 'RTE', 'L1', 'Mariner/Tc1', 'hAT', 'Maverick', 'Helitron', 'Harbinger', 'P']:
			return annotation
		elif annotation == 'SINE2/tRNA':
			return 'tRNA'
		elif re.search('enspm', annotation, re.I):
			return 'CACTA'
		elif annotation == 'SINE1/7SL':
			return '7SL'
		elif annotation in ['Penelope', 'Sm1']:
			return 'Penelope'
		elif re.search('mudr', annotation, re.I):
			return 'Mutator'
		else:
			if self.CLASS == 'Unknown' and self.ORDER == 'Unknown':
				print '!!! Unknown! ' + annotation
			return 'Unknown'

class RepeatGroupStats(object):

	def __init__(self, classification, blockSize):
		self.classification = classification
		self.length = blockSize
		self.nElements = 1
	
	def __iadd__(self, blockSize):
		self.length += blockSize
		self.nElements += 1
		return self
	
	def getAvgLength(self):
		self.avgLength = float(self.length)/float(self.nElements)

	def __hash__(self):
		return hash(self.classification)
	
	def __eq__(self, other):
		return self.classification == other.classification
	
	def __str__(self):
		return '{}\t{:d}\t{:d}'.format(self.classification, self.nElements, self.length)

class RepeatStats(object):

	def __init__(self):
		self.classes = {}
		self.orders = {}
		self.supers = {}
		self.families = {}

	def addRepeatCopy(self, repeat, blockSize): # repeat is of type Repeat
		if repeat.CLASS in self.classes and repeat.CLASS is not None:
			self.classes[repeat.CLASS] += blockSize
		elif repeat.CLASS is not None:
			self.classes[repeat.CLASS] = RepeatGroupStats(repeat.CLASS, blockSize)
		if repeat.ORDER in self.orders and repeat.ORDER is not None:
			self.orders[repeat.ORDER] += blockSize
		elif repeat.ORDER is not None:
			self.orders[repeat.ORDER] = RepeatGroupStats(repeat.ORDER, blockSize)
		if repeat.SUPER in self.supers and repeat.SUPER is not None:
			self.supers[repeat.SUPER] += blockSize
		elif repeat.SUPER is not None:
			self.orders[repeat.SUPER] = RepeatGroupStats(repeat.SUPER, blockSize)
		if repeat.FAMILY in self.families and repeat.FAMILY is not None:
			self.families[repeat.FAMILY] += blockSize
		elif repeat.FAMILY is not None:
			self.families[repeat.FAMILY] = RepeatGroupStats(repeat.FAMILY, blockSize)
			
	def __iadd__(self, repeatTuple): #overloads stats += (repeat, blockSize)
		self.addRepeatCopy(repeatTuple[0], repeatTuple[1])
		return self
	
	def __str__(self):
		outputString = 'Group\tnElements\tLength\n===CLASSES===\n'
		
		for myClass in self.classes:
			outputString += '{}\n'.format(self.classes[myClass])
		
		outputString = outputString[:-2] + '\n===ORDERS===\n'
	
		for order in self.orders:
			outputString += '{}\n'.format(self.orders[order])
		
		outputString = outputString[:-2] + '\n===SUPERFAMILIES===\n'
		
		for mySuper in self.supers:
			outputString += '{}\n'.format(self.supers[mySuper])
		
#		outputString = outputString[:-2] + '\n===FAMILIES===\n'
		
#		for family in self.families:
#			outputString += '{}\n'.format(self.families[family])
		
		return outputString
