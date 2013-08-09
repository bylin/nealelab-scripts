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
		self.start = int(Data[1]) #Indices of the repeat relative to the start of the sequence 
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
		self.motif = Data[14]#The repeat motif itself
		self.sequence = Data[15]#The entire sequence
		self.length = self.end - self.start#Alignment length
					
	def __str__(self): 
		#for printing in mysql format
		return str(self.seq_name) + "\t" + str(self.start) + "\t" + str(self.end) + "\t" + str(self.period) + "\t" + str(self.copy_num)+ "\t" + \
		str(self.size_consensus)+ "\t" + str(self.percent_matches) + "\t" + str(self.percent_indels) + "\t" + str(self.score) + "\t" + \
		str(self.a) + "\t"+ str(self.c) + "\t" + str(self.t) + "\t"+ str(self.g) + "\t" + str(self.entropy) + "\t" + str(self.length)+ \
		"\t" + self.motif
	
	def __cmp__(self,other):
		return self.score > other.score #blast score, metric used to compare redundant bp

	def __eq__(self,other):
		return self.seq_name == other.seq_name

	def __hash__(self):
		return hash(self.sequence)
	
	def test(self):
		return str(self.start)+"\t"+str(self.end)+"\t"+str(self.score)

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

