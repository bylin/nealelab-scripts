#!/bin/usr/python
#Author: Jacob Zieve, Brian Lin

"""Defines different "types" of data inputs that are handled
int the stats scritps. Many of the classes include overloaded equality, hashing
and comparison operators to be used in the script"""

from collections import defaultdict
from Bio import SeqIO
import re
#used as main data structure for sequence

class Repeat(object):

	def __init__(self, code, name):
		self.NAME = name
		self.CLASS = self.determineClass(code)
		self.ORDER = self.determineOrder(code)
		self.SUPER = self.determineSuper(code)
		self.FAMILY = self.determineFamily(name)
	
	def __str__(self):#for debug
		return '{:s}: {:s}, {:s}, {:s}, {:s}\n'.format(self.NAME, self.CLASS, self.ORDER, self.SUPER, self.FAMILY)

	def __eq__(self, other):
		return self.name == other

	def __hash__(self): #used in set() function later
		return hash(self.name)

class repBaseRepeat: #DEPRECATED

	def __init__(self,name, accession, superfamily, order, _class, species, length):
		self.name = name
		self.accession = accession
		self.superfamily = superfamily
		self.order = order
		self._class = _class
		self.species = species
		self.length = length
		
	def __str__(self):
		return "name: "+self.name+"\nac: "+self.accession+"\nfam: "+self.superfamily+"\n\
		sp: "+self.species+"\nlen: "+str(self.length)+"\n"

class RepbaseRepeat(Repeat):
	
	def __init__(self, name, accession, superfamily, order, _class, species, length):
		self.NAME = name
		self.ORDER = order
		self.SUPER = superfamily
		self.CLASS = _class
		self.FAMILY = name
		self.accession = accession
		self.species = species
		self.length = length

def conversion(r1): #TEMPORARY, SHOULD BE DEPRECATED SOON
	r2 = RepbaseRepeat(r1.name, r1.accession, r1.superfamily, r1.order, r1._class, r1.species, r1.length)
	return r2

class WickerRepeat(Repeat):
	# Example: myRepeat = WickerRepeat('RLX')
	
	def determineClass(self, code):
		myClass = code[0]
		if myClass == 'R': return 'I'
		elif myClass == 'D': return 'II'
		elif myClass == 'X': return 'Confused'
		elif code in ['PotentialHostGene', 'NoCat', 'SSR', 'rDNA']: return code
	
	def determineOrder(self, code):
		order = code[1]
		if self.CLASS == 'I':
			if order == 'X': return 'Unknown'
			elif order == 'L': return 'LTR'
			elif order == 'I': return 'LINE'
			elif order == 'Y': return 'DIRS'
			elif order == 'S': return 'SINE'
			elif order == 'P': return 'Penelope'
		elif self.CLASS == 'II':
			if order == 'X': return 'Unknown'
			elif order == 'T': return 'TIR'
			elif order == 'Y': return 'Crypton'
			elif order == 'H': return 'Helitron'
			elif order == 'M': return 'Maverick'
		elif self.CLASS in ['Confused', 'PotentialHostGene', 'NoCat', 'rRNA', 'SSR', 'rDNA']:
			return 'Unknown'
		else:
			print self.NAME, self.CLASS
			exit("Unhandled case: " + code)
		
	def determineSuper(self, code):
		super = code[2]
		if super == 'G': return 'Gypsy'
		elif super == 'C': return 'Copia'

	def determineFamily(self, name):
		if re.search('_[A-Z]', name):
			return '_'.join(name.split('_')[:-1])
		return name

# TODO: flesh out class, order, super functions. Repbase can be tricky.
class PierRepeat(Repeat):
	
	def __init__(self, name):
		self.NAME = self.determineFamily(name)
		self.CLASS = 'I'
		self.ORDER = 'LTR'
		self.SUPER = self.determineSuper(name)
		self.FAMILY = self.NAME

	def determineSuper(self, annotation):
		if re.search('ouachita|bastrop|ozark|appalachian|angelina|talladega', annotation, re.I):
			return 'Gypsy'
		if re.search('cumberland|pineywoods|conagree', annotation, re.I):
			return 'Copia'
		return 'Unknown'

	def determineFamily(self, name): #truncate subfamily
		return '_'.join(name.split('_')[:-1])

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
		self.totalLength = 0
		self.nSeqs = 0

	def addRepeatCopy(self, repeat, blockSize): # repeat is of type Repeat
		self.totalLength += blockSize
		self.nSeqs += 1
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
			self.supers[repeat.SUPER] = RepeatGroupStats(repeat.SUPER, blockSize)

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
		outputString += '\n===TOTAL===\n\t{}\t{}\n'.format(self.nSeqs, self.totalLength)
		return outputString
