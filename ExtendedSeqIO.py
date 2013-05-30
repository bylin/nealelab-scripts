from Bio import SeqIO
from Bio import SeqRecord
import re

# ExtendedSeqIO.py
# Author: Brian Lin
# Description: Extended versions of a file object. Supports simpler, abstracted sequence handling, matching, and extraction using 

class ExtendedFastaFile(file):

	def seqRecordList(self):
		return SeqIO.parse(self, 'fasta')
	
	def resetFilePointer(self):
		self.seek(0)

	def getSeqUsingDescription(self, description):
		for seq in self.seqRecordList():
			if seq.description == description:
				mySeq = seq
		self.resetFilePointer()
		if (mySeq is None): print "Sequence %s not found" % description
		else: return mySeq
	
	def getSeqsUsingDescriptionList(self, descriptionList):
		seqList = []
		numDescriptionsLeft = len(descriptionList)
		for seq in self.seqRecordList():
			if seq.description in descriptionList:
				seqList.append(seq)
				numDescriptionsLeft -= 1
		self.resetFilePointer()
		if numDescriptionsLeft != 0: print "%s sequences not found" % (numDescriptionsLeft)
		return seqList

	def getSeqsUsingRegex(self, regex):
		seqList = []
		for seq in self.seqRecordList():
			if re.search(regex, seq.description, re.IGNORECASE): 
				seqList.append(seq)
		self.resetFilePointer()
		return seqList

class ExtendedBlastTabFile(file):

	def findResultsUsingQuery(query):
		return
	def findResultsUsingSubject(subject):
		return

