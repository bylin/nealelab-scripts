#!/usr/bin/python
from Bio import SeqIO
from Bio import SeqRecord
import re

# ExtendedSeqIO.py
# Author: Brian Lin
# Description: Extended versions of a file object. Supports simpler, abstracted sequence handling, matching, and extraction using Bio::SeqIO

class ExtendedFastaFile(file):
	
	def seqRecordList(self):
		return SeqIO.parse(self, 'fasta')
	
	def resetFilePointer(self):
		self.seek(0)


	def getSeqUsingId(self, id):
		mySeq = None
		for seq in self.seqRecordList():
			if seq.id == id:
				mySeq = seq
		self.resetFilePointer()
		if (mySeq is None): print "Sequence %s not found" % id
		else: return mySeq
	
	def getSeqUsingName(self, name):
		mySeq = None
		for seq in self.seqRecordList():
			if seq.name == name:
				mySeq = seq
		self.resetFilePointer()
		if (mySeq is None): print "Sequence %s not found" % name
		else: return mySeq
	
	def getSeqUsingDescription(self, description):
		mySeq = None
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
	
	def createFastaFromRecord(self, seq):
		print seq.format("fasta")
	
	def createMultiFastaFromRecords(self, seqs):
		for seq in seqs: self.createFastaFromRecord(seq)

class ExtendedBlastTabFile(file):

	def findResultsUsingQuery(query):
		return
	def findResultsUsingSubject(subject):
		return


########
def stub():
	pier = ExtendedFastaFile('/home/jjzieve/Pita_Genome-0.9_Repeats/databases/pier-1.3.fa')
	gypsySeqs = pier.getSeqsUsingRegex('Gypsy')
	#for seq in gypsySeqs: print seq.id
	pier.createMultiFastaFromRecords(gypsySeqs)

if __name__=="__main__":
	stub()
