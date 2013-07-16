#!/usr/bin/python
from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqUtils import GC
from classes import repBaseRepeat
from subprocess import PIPE,Popen
from pickler import sendToPickleJar,getFromPickleJar
import re,pickler,sys

# FileIO.py
# Author: Brian Lin, Jake Zieve
# Description: Extended versions of a file object. Supports simpler, abstracted sequence handling, matching, and extraction using Bio::SeqIO
class ExtendedFile(file):
	
	def __init__(self, filename):
		file.__init__(self, filename)
		self.filename = filename
	
	def resetFilePointer(self):
		self.seek(0)

class Fasta(ExtendedFile):

	def seqRecordList(self):
		return SeqIO.parse(self, 'fasta')
	
	@staticmethod
	def writeToFile(seq, filename):
<<<<<<< .merge_file_ZULxLs
		handle = open(filename, 'a')
		SeqIO.write(seq, handle, 'fasta')
=======
		SeqIO.write(seq, filename, 'fasta')
>>>>>>> .merge_file_6UHb8o

	def printSeqsStats(self):
		self.printNumSeqs()
		self.printAvgSeqLen()
		self.printMedSeqLen()
		self.printN50()
		self.printShortestSeqLen()
		self.printLongestSeqLen()
		self.printTotalBp()
		self.printAvgGCcontent()
	
	def printNumSeqs(self):
		print 'Number of sequences:\t{:,d}'.format(self.getNumSeqs())
	
	def getNumSeqs(self):
		num = sum(1 for seq in self.seqRecordList())
		self.resetFilePointer()
		return num
	
	def printAvgSeqLen(self):
		print 'Average sequence length:\t{:,d}'.format(self.getAvgSeqLen())
	
	def getSeqLens(self):
		lenList = [len(seq) for seq in self.seqRecordList()]
		self.resetFilePointer()
		return lenList

	def getAvgSeqLen(self):
		total = sum(self.getSeqLens())
		num = self.getNumSeqs()
		return  total/num
	
	def printMedSeqLen(self):
		print 'Median sequence length:\t{:,d}'.format(self.getMedSeqLen())

	def getMedSeqLen(self):
		return sorted(self.getSeqLens())[self.getNumSeqs()/2]

	def printN50(self):
		print 'N50:\t{:d}'.format(self.getN50())
			
	def getN50(self):
		cumulative = 0
		lenList = self.getSeqLens()
		for s in sorted(lenList):
			cumulative += s
			if (cumulative >= (sum(lenList)/2.0)):
				return s
			
	def printShortestSeqLen(self):
		print 'Shortest sequence length:\t{:,d}'.format(self.getShortestSeqLen())
	
	def getShortestSeqLen(self):
		return sorted(self.getSeqLens())[0]

	def printLongestSeqLen(self):
		print 'Longest sequence length:\t{:,d}'.format(self.getLongestSeqLen())
	
	def getLongestSeqLen(self):
		return sorted(self.getSeqLens(),reverse=True)[0]

	def printTotalBp(self):
		print 'Total bp:\t{:,d}'.format(self.getTotalBp())

	def getTotalBp(self):
		return sum(self.getSeqLens())	
			
	def printAvgGCcontent(self):
		print 'Average GC content:\t{:.2f}%'.format(self.getAvgGCcontent())
		
	def getAvgGCcontent(self):
		gc = sum(GC(seq.seq) for seq in self.seqRecordList())
		self.resetFilePointer()
		return gc/self.getNumSeqs()
		
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
	
	def removeRedundantSeqs(self):
		outfile = '{}.nonredundant'.format(self.filename)
		seqList = []
		for seq in self.seqRecordList():
<<<<<<< .merge_file_ZULxLs
			if seq.description not in seqList:
				seqList.append(seq.description)
=======
			if seq not in seqList:
				seqList.append(seq)
>>>>>>> .merge_file_6UHb8o
				self.writeToFile(seq, outfile)

class BlastOutTab(ExtendedFile):

	def findResultsUsingQuery(query):
		return
	def findResultsUsingSubject(subject):
		return

class EMBL(ExtendedFile):
	
	def parse(self):
		emblDict = {}
		for record in self.read().split('//'):
			key_value = self.parseLines(record)
			emblDict[key_value[0]]=key_value[1]
		self.resetFilePointer()
		return emblDict

	def parseLines(self,record):
		name = ""
		accession = ""
		superfamily = ""
		species = ""
		length = 0	
		for line in record.split('\n'):
			fields = line.split()
			if len(fields) < 1: continue
			if fields[0] == 'ID':
				name = fields[1]
			if fields[0] == 'AC':
				accession = fields[1].strip(';.')
			if fields[0] == 'KW' and superfamily == "":
				superfamily = fields[1].strip(';')
			if fields[0] == 'OS':
				species = ' '.join(fields[1:])
			if fields[0] == 'SQ':
				length = int(fields[2])
		order = self.getOrder(superfamily)
		_class = self.getClass(order)
		return (name,repBaseRepeat(name,accession,superfamily,order,_class,species,length))
	
	def getOrder(self,superfamily):
		LINE = ['CRE','NeSL','R4','R2','L1','RTE','I','Jockey','CR1','Rex1','RandI','Tx1','RTEX','Crack','Nimb','Proto1','Proto2','RTETP','Hero','L2','Tad1','Loa','Ingi','Outcast','R1','Daphne','L2A','L2B','Ambal','Vingi','Kiri']
		TIR = ['Mariner/Tc1','hAT','MuDR','EnSpm','piggyBac','P','Merlin','Harbinger','Transib','Novosib','Kolobok','ISL2EU','Sola','Zator','Ginger1','Ginger2/TDD','Academ','Zisupton','IS3EU','Dada']
		if superfamily in ['Copia','Gypsy','BEL','DIRS']:
			return 'LTR'
		elif superfamily in ['ERV1','ERV2','ERV2','ERV3','ERV4','Lentivirus']:
			return 'Endogenous Retrovirus'
		elif superfamily in LINE:
			return 'LINE'
		elif superfamily in ['SINE','SINE1/7SL','SINE2/tRNA','SINE3/5S','SINE4']:
			return 'SINE'
		elif superfamily == 'Penelope':
			return 'Penelope'
		elif superfamily in TIR:
			return 'TIR'
		elif superfamily in ['Crypton','CryptonA','CryptonF','CryptonI','CryptonS','CryptonV']:
			return 'Crypton'
		elif superfamily == 'Helitron':
			return 'Helitron'
		elif superfamily == 'Polinton':
			return 'Maverick/Polinton'
		elif superfamily in ['SAT','MSAT']:
			return 'Satellite'
		elif superfamily in ['rRNA','tRNA','snRNA']:
			return 'Pseudogene/RNA'
		elif superfamily in ['DNA Virus', 'Caulimoviridae']:
			return 'Integrated Virus'
		else:
			return 'Unknown'

	def getClass(self,order):
		if order in ['LTR','Endogenous Retrovirus','LINE','SINE','Penelope']:
			return 'I'
		elif order in ['TIR', 'Crypton', 'Helitron', 'Maverick/Pointon']:
			return 'II'
		elif order in ['Satellite','Pseudogene/RNA','Integrated Virus']:
			return 'Other'
		else:
			return 'Unknown'
		
			
####Temporary###
def main():
	#repBase = ExtendedEMBLFile('/home/jjzieve/Pita_Genome-0.9_Repeats/databases/RepBase18.03.embl/all.ref').parse()
	#pickler.sendToPickleJar(repBase,'repBaseDict.pkl')
	#repBase = pickler.getFromPickleJar('repBaseDict.pkl')
	#print len(repBase) 
	#fasta = ExtendedFastaFile(sys.argv[1])
	#fasta.printSeqsStats()
	return

if __name__=="__main__":
	main()
