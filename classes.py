#!/bin/usr/python
#Author: Jacob Zieve, Brian Lin

"""Defines different "types" of data inputs that are handled
int the stats scritps. Many of the classes include overloaded equality, hashing
and comparison operators to be used in the script"""

from repeat import *
from collections import defaultdict
from Bio import SeqIO
import re
#used as main data structure for sequence
	
class trfHit:
	def __init__(self, Data=[]):
		self.seq_name = Data[0]
		self.start = int(Data[1]) #Intrices of the repeat relative to the start of the sequence 
		self.end = int(Data[2])
		self.period = int(Data[3])#Period size of the repeat
		self.copy_num = float(Data[4])#Number of copies aligned with the consensus pattern
		self.size_consensus = int(Data[5])#Size of consensus pattern (may differ slightly from the period size)
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
		return "'"+str(self.seq_name) + "'," + str(self.start) + "," + str(self.end) + "," + str(self.period) + "," + str(self.copy_num)+ "," + \
		str(self.size_consensus)+ "," + str(self.percent_matches) + "," + str(self.percent_indels) + "," + str(self.score) + "," + \
		str(self.a) + ","+ str(self.c) + "," + str(self.t) + ","+ str(self.g) + "," + str(self.entropy) + "," + "'" + str(self.motif)+ "'" + \
		"," + str(self.length)
	
	def __cmp__(self,other):
		return self.score > other.score #blast score, metric used to compare redundant bp

	def __eq__(self,other):
		return self.seq_name == other.seq_name

	def __hash__(self):
		return hash(self.sequence)
	
	def test(self):
		return str(self.start)+","+str(self.end)+","+str(self.score)

class BlastHit(object):
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

class HmmHit(object):
	def __init__(self, name, start, end, score):
		self.start = start
		self.end = end
		self.name = name
		self.score = score
	def __cmp__(self, other): # for efficient sorting of hits within a raw sequence
		if self.start == other.start:
			if self.end == other.end:
				return 0
			elif self.end > other.end:
				return 1
			elif self.end < other.end:
				return -1
		elif self.start > other.start:
			return 1
		else: return -1
	def __str__(self):
		return '{}[{}:{}]'.format(self.name, self.start, self.end)
	
