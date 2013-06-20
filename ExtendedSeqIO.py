#!/usr/bin/python
from Bio import SeqIO
from Bio import SeqRecord
from classes import repBaseRepeat
import re,pickler

# ExtendedSeqIO.py
# Author: Brian Lin
# Description: Extended versions of a file object. Supports simpler, abstracted sequence handling, matching, and extraction using Bio::SeqIO
class ExtendedFile(file):
	
	def resetFilePointer(self):
		self.seek(0)

class ExtendedFastaFile(ExtendedFile):
##def n50(len_list):
#	cumulative = 0
#        for s in sorted(len_list):
#                cumulative += s
#                if (cumulative >= (sum(len_list)/2.0)):
#                        return s
#
#def len_sequences(fasta_file):
#        sanity_check(fasta_file)
#        len_dict = {}
#        for seq_record in SeqIO.parse(fasta_file, "fasta"):
##               if len(seq_record)-seq_record.seq.count("N")==75791 or \
##               len(seq_record)-seq_record.seq.count("N")==67475:
##                       print ">"+seq_record.id
##                       print seq_record.seq
##                       print "\n"
#                len_dict[seq_record.id]=len(seq_record)-seq_record.seq.upper().count("N")#omits unknown bases
##       print sorted(len_dict.values())[len(len_dict.values())-1]
##       print sorted(len_dict.values())[len(len_dict.values())-2]
#        return len_dict
#	
#def nucl_lengths(fasta_file):
#        sanity_check(fasta_file)
#        A_count = []
#        C_count = []
#        T_count = []
#        G_count = []
#        N_count = []
#        seq_count = 0
#
#        for seq_record in SeqIO.parse(fasta_file, "fasta"):
#                seq_count += 1
#                seq = seq_record.seq.upper()
#                #print seq_record.id#the header (e.g. >fasta_name...)
#                #print repr(seq_record.seq)#how its represented with biopython
#                A_count.append(seq.count("A"))
#                C_count.append(seq.count("C"))
#                T_count.append(seq.count("T"))
#                G_count.append(seq.count("G"))
#
#        return(A_count,C_count,T_count,G_count)
#
	def seqRecordList(self):
		return SeqIO.parse(self, 'fasta')
	
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

class ExtendedBlastTabFile(ExtendedFile):

	def findResultsUsingQuery(query):
		return
	def findResultsUsingSubject(subject):
		return

class ExtendedEMBLFile(ExtendedFile):
	
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
#		print name
#		print accession
#		print superfamily
#		print order
#		print _class
#		print species
#		print length
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
	repBase = pickler.getFromPickleJar('repBaseDict.pkl')
	for k,v in repBase.items(): print k 

if __name__=="__main__":
	main()
