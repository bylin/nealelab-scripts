# findORFsFromFasta.py
# Author: Brian Lin

# Start with a FASTA file containing all sequences in question.
# Use USEARCH to find all ORFs
# For each ORF in each sequence, run interpro
# Get lowest e-value set of ORFs
# Extract DNA sequence from FASTA, translate to Protein

import ExtendedSeqIO, re, os
from Bio import SeqIO
from subprocess import Popen,PIPE

def checkAndReturnArguments():
	if (len(sys.argv) != 1):
		print "Usage: findORFsFromFasta.py [input fasta file]"
	return sys.argv[1]

def loadModules():
	subprocess.call('module load bigmem')
	subprocess.call('module load usearch')
	subprocess.call('module load jre')
	subprocess.call('module load interproscan/5-RC6')

def getAllORFsOfSeqAsFile(seq):
	# write sequence to temporary file so we can run findorfs on it
	seqFileName = seq.name+'.fa'
	seqFileHandle = open(seqFileName, 'w')
	SeqIO.write(seq, seqFileHandle, 'fasta')
	findorfs_command = 'usearch -findorfs ' + seqFileName ' -output ' + seqFileName+'.orfs -xlat'
	subprocess.call(findorfs_command)
	return seqFileName+'.orfs'

def runInterproAndGetGFFFile(orfFile):
	orfFilePath = '/share/nealedata/blin/'+orfFile
	interpro_command = 'interproscan.sh -appl PfamA -f GFF3 -iprlookup -goterms -i '+orfFilePath
	subprocess.call(interpro_command);
	return orfFile+'.gff3'

def filterTrueORFsFromInterproGFFFile(gffFile):

	return gffFile+'.filtered'

def cleanFiles(seqID):
	os.remove(seqID+'.fa')
	os.remove(seqID+'.fa.orfs')
	os.remove(seqID+'.fa.orfs.gff3')


if __name__ == '__main__':
	filename = checkAndReturnArguments()
	loadModules()
	for seq in SeqIO.parse(filename, 'fasta'):
		orfFile = getAllORFsOfSeqAsFile(seq)
		gffFile = runInterproAndGetGFFFile(orfFile)
		trueORFs = filterTrueORFsFromInterproGFFFile(gffFile)
		writeORFSeqsToFile(seq, trueORFs)
		cleanFiles(seq.name)