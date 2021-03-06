#!/usr/bin/python
# Author: Brian Lin
# Encodes marks on each individual base pair within a sequence. Meant for comparison and elimination of redundant regions in a sequence.

import math, sys

class BitwiseSeq(object):

	def __init__(self, length):
		self.length = length
		self.bitint = 0
		self.bitmax = pow(2, length) - 1

	def verifyBlockRange(self, block):
		(startPosition, endPosition) = block
		if startPosition <= 0 or endPosition > (self.length + 1) or endPosition < startPosition:
			raise ValueError
		return True

	def __iadd__(self, blockBitint):
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

	def addBlockReturnDifference(self, block):
		assert self.verifyBlockRange(block) is True
		blockBitint = convertBlockToBitint(block)
		oldSelfBitint = self.bitint
		self += blockBitint
		diffBitint = oldSelfBitint^self.bitint
		start = block[0]
		return countFlippedBits(diffBitint >> (start - 1))

	def __str__(self):
		return '{:d}:{:d}'.format(self.length, self.bitint)

# Kernighan popcount
def countFlippedBits(bitint):
	count = 0 
	while bitint:
		bitint &= bitint - 1
		count += 1
	return count

def convertIntToBitint(number):
	return (pow(2, number - 1))

def convertBlockToBitint(block):
	(start, end) = block
	left = pow(2, end-1)
	right = pow(2, start-1)
	return left-right
