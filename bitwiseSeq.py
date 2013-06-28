#!/usr/bin/python
# Author: Brian Lin
# Encodes marks on each individual base pair within a sequence. Meant for comparison and elimination of redundant regions in a sequence.

import math

class BitwiseSeq(object):

	def __init__(self, length):
		self.length = length
		self.bitint = 0
		self.bitmax = pow(2, length) - 1

	def verifyBlockRange(self, block):
		(startPosition, endPosition) = block
		if startPosition <= 0 or endPosition > self.length or endPosition < startPosition:
			return False
		return True

	def __iadd__(self, block):
		assert self.verifyBlockRange(block) is True
		blockBitint = convertBlockToBitint(block)
		self.bitint = self.bitint|blockBitint
		return self

	def getBitAtPosition(self, position):
		assert position <= self.length - 1
		bitint = self.bitint
		position = self.length - 1
		while bitint > 0:
			bitint -= pow(2, position)
		return 1 if bitint == 0 else 0

	def findBlockDifference(self, block):
		selfBitint = self.bitint
		blockBitint = convertBlockToBitint(block)
		diffBitint = (~selfBitint) & blockBitint
		return countFlippedBits(diffBitint)

	def __str__(self):
		return '{:d}:{:d}'.format(self.length, self.bitint)

def countFlippedBits(bitint): #counts number of flipped bits in O(n), n = # of flipped bits
		nflipped = 0
		while bitint != 0:
			bitint -= 2**int(math.log(bitint, 2))
			nflipped += 1
		return nflipped

def convertIntToBitint(number):
	return (pow(2, number - 1))


def convertBlockToBitint(block):
	(start, end) = block
	bitstring = ''
	for i in range(start - 1, end):
		bitstring += '1'
		i += 1
	for i in range(0, start):
		bitstring += '0'
	return int(bitstring, 2)

def convertBlockToBitint2(block):
	(start, end) = block
	return sum(convertIntToBitint(base) for base in range(start, end))
