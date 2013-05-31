#!/bin/usr/python
#Author: Jacob Zieve, Brian Lin

"""Defines different "types" of data inputs that are handled
int the stats scritps. Many of the classes include overloaded equality, hashing
and comparison operators to be used in the script"""

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

class denovoRepeat(object):
	def __init__(self,num,cat):
		self.num = num #unique id for each consensus sequence
		self.cat = cat#categorie (e.g. RLX = LTRs)
	#	self.seq = seq #where found on the original sequence
	#	self.status = status #complete or incomplete
	
	def __str__(self):#for debug
		return 'num: {:s}\ncategory: {:s}\n'.format(self.num,self.cat)
	def __eq__(self,other):
		return self.cat == other
	def __hash__(self): #used in set() function later
		return hash(self.cat)

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

class PierRepeat(object):
	# Annotation examples: "Gypsy", "LTR", "Class II", "L1"
	def __init__(self, annotation):
		self.CLASS = determineClass(annotation)
		self.ORDER = determineOrder(annotation)
		self.SUPER = determineSuper(annotation)

	def determineClass(annotation):
		if annotation in ['Gypsy', 'Copia', 'LTR', 'Class I', 'Penelope', 'LINE', 'SINE', 'SINE1/7SL', 'SINE2/tRNA', 'Gymny', 'L1', 'Non-LTR Retrotransposon', 'LTR Retrotransposon', 'Endogenous Retrovirus', 'RTE']:
			return 'I'
		elif annotation in ['Class II', 'Helitron', 'DIR', 'DIRS', 'P', 'Maverick', 'hAT', 'TIR', 'Mariner/Tc1', 'Harbinger', 'MuDR', 'EnSpm']:
			return 'II'
		else:
			return 'Uncategorized'

	def determineOrder(annotation):
		if annotation in ['Gypsy', 'Copia', 'LTR', 'LTR Retrotransposon', 'Endogenous Retrovirus', 'Gymny']:
			return 'LTR'
		elif annotation in ['LINE', 'L1', 'RTE']:
			return 'LINE'
		elif annotation in ['SINE', 'SINE1/7SL', 'SINE2/tRNA']:
			return 'SINE'
		elif annotation in ['hAT', 'P', 'Mariner/Tc1', 'Harbinger', 'MuDR', 'TIR']:
			return 'TIR'
		elif annotation in ['DIR', 'DIRS']:
			return 'DIRS'
		elif annotation == 'EnSpm':
			return 'EnSpm'
		elif annotation in ['Penelope', 'SINE', 'Helitron', 'Maverick']:
			return annotation
		elif annotation in ['Class I', 'Class II']:
			return 'Unknown'

	def determineSuper(annotation):
		if annotation in ['Gymny', 'Gypsy']:
			return 'Gypsy'
		elif annotation in ['Copia', 'DIRS', 'Penelope', 'RTE', 'L1', 'Mariner/Tc1', 'hAT', 'Maverick', 'Helitron', 'Harbinger', 'P']:
			return annotation
		elif annotation == 'SINE2/tRNA':
			return 'tRNA'
		elif annotation == 'SINE1/7SL':
			return '7SL'
		elif annotation == 'MuDR':
			return 'Mutator'
		else:
			return 'Unknown'

