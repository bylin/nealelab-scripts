#!/usr/bin/python
# Author: Brian Lin
# Encodes marks on each individual base pair within a sequence. Meant for comparison and elimination of redundant regions in a sequence.

class BitwiseSeq(object):

	def __init__(self, length):
		self.bits = convertIntToBits(length)
