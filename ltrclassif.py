#!/usr/bin/python
# Author: Brian Lin
import subprocess, argparse, timestamp, re, os, classes
from FileIO import Fasta

def parseArgs():
	parser = argparse.ArgumentParser(description="Uses probabilistic models to find likely LTR retrotransposon candidates")
	parser.add_argument('fasta', metavar='[fasta]', help="Input .fasta file, containing transposable elements")
	parser.add_argument('-blast', '--use_blast_evidence', help="Perform a BLAST search for additional evidence. Must load ncbi-blast module before executing program with this flag!")
	parser.add_argument('-noclean', '--keep_temporary_files', action="store_true", help="Don't clean up hmm and blast output files")
	parser.add_argument('-o', '--output_fasta_file', default='classified-{}.out'.format(timestamp), help="Tab-delimited output file containing new classifications")
	return parser.parse_args()

def main():
	if args.use_blast_evidence:
		blastHits = blastAndParse()
	runHmmer()
	rawSeqs = {}
	for seq in Fasta(orfFile):
		rawSeqs[seq.id] = ExtendedSeqRecord(seq)
	parseHmmHits(rawSeqs)
	examineDomainsAndClassify(rawSeqs)
	clean()

def blastAndParse():
	blastDBLens = {}
	for seq in Fasta(blastDB): dbLens[seq.id] = len(seq)
	runBlast()
	for line in open(blastHitFile, 'r'):
		blastHit = BlastHit(line.split('\t'))
		if blastHit.pid >= 80 and blastHit.alen > blastDBLens[blastHit.id]:
			continue

def runBlast():
	checkBlastVersion()
	makeDBCmd = "makeblastdb -dbtype nucl -in {0} -title {0} -out {0}".format(blastDB)
	#blastCmd = "blastn -subject {} -query {} -outfmt '6 qseqid sseqid length pident evalue bitscore' -out {}".format(args.fasta, blastDB, blastHitFile)
	blastCmd = "blastn -db {} -query {} -outfmt 6 -out {}".format(args.fasta, blastDB, blastHitFile)
	subprocess.call(blastCmd, shell=True)

def checkBlastVersion():
	proc = subprocess.Popen(['blastn', '-version'], stdout=subprocess.PIPE)
	stdout = proc.communicate()[0]
	if not re.search('2.2.27', stdout):
		print 'Did you remember to load the ncbi-blast module?'

def runHmmer():
	hmmernCmd = 'hmmern.py -hmm {} -i {} -o {}'.format(hmmProfilesForLTRs, args.fasta, hmmHitFile)
	subprocess.call(hmmernCmd, shell=True)

def parseHmmHits(rawSeqs):
	for line in open(hmmHitFile):
		tabs = re.findall('\S+', line)
		if tabs[0][0] == '#': continue
		rawSeq = tabs[0]
		hmmName = tabs[3]
		try: orfstart = int(re.findall(':(\d+)-', tabs[-1])[0])
		except:
			print line
			continue
		start = int(tabs[19]) * 3 - 2 + orfstart
		end = int(tabs[20]) * 3 - 5 + orfstart
		score = float(tabs[13])
		hmmhit = classes.HmmHit(hmmName, start, end, score)
		try: rawSeqs[rawSeq].addHit(hmmhit)
		except KeyError:
			print '{} not in rawSeqs'.format(rawSeq)

class ExtendedSeqRecord(object):
	def __init__(self, seq):
		self.seq = seq
		self.hits = []

	def addHit(self, hmmhit):
		self.hits.append(hmmhit)

def examineDomainsAndClassify(rawSeqs):
	outputFile.write('Sequence\tOld classification\tNew classification\n')
	for xseq in rawSeqs:
		if len(rawSeqs[xseq].hits) == 0: continue
		string = "{}: ".format(xseq)
		string += ' '.join('{}({})'.format(hit,hit.score) for hit in sorted(rawSeqs[xseq].hits))
		classif, score = classify(sorted(rawSeqs[xseq].hits), [], 0)
		old = rawSeqs[xseq].seq.description.split('\t')[1]
		if classif and old != classif:
			print xseq, string
			outputFile.write('{}\t{}\t{}\t{}\n'.format(xseq, old, classif, score))

def classify(hitList, classifs, currentScore):
	for hit,i in zip(hitList, range(0,len(hitList))):
		if hit.name[0:6] == 'POL-GA':
			classify_GAG(hitList[i:], classifs, hit.score)
		if hit.name[0:6] == 'POL-AP':
			classify_AP(hitList[i:], classifs, hit.score)
		elif hit.name[0:6] == 'POL-RT':
			classify_RT(hitList[i:], classifs, hit.score)
		elif hit.name[0:6] == 'POL-RH':
			classify_RH(hitList[i:], classifs, hit.score)
		elif hit.name[0:7] == 'POL-INT':
			classify_INT(hitList[i:], classifs, hit.score)
	return findBestClassif(classifs)

def findBestClassif(classifs):
	if len(classifs) > 0:
		classifs = sorted(classifs, key=lambda x:x[1], reverse=True)
		#print bestClassif[1]
		#print '\t' + ','.join('{}:{}'.format(c[0], c[1]) for c in classifs)
		return classifs[0][0], classifs[0][1] # return classif and score
	else:
		return '', -1

# Classify_: represents states in a state diagram. 
# This is coded as a sort of DFA, but all paths are taken to maximize the cumlative HMM hit score of a series of linked hits.
def classify_GAG(hitList, classifs, currentScore):
	for hit, i in zip(hitList[i:], range(1, len(hitList))):
		if hit.name[0:6] == 'POL-AP':
			classify_AP(hitList[i:], classifs, hit.score + currentScore)
		elif hit.name[0:6] == 'POL-RT':
			classify_RT(hitList[i:], classifs, hit.score + currentScore)
		elif hit.name[0:6] == 'POL-RH':
			classify_RH(hitList[i:], classifs, hit.score + currentScore)
		elif hit.name[0:7] == 'POL-RTRH':
			classify_RTRH(hitList[i:], classifs, hit.score + currentScore)
		elif hit.name[0:6] == 'POL-IN':
			classify_INT(hitList[i:], classifs, hit.score + currentScore)

def classify_AP(hitList, classifs, currentScore):
	for hit, i in zip(hitList[1:], range(1, len(hitList))):
		if hit.start-hitList[0].end > 2000:
			continue
		if hit.name[0:6] == 'POL-RT':
			classify_RT(hitList[i:], classifs, hit.score + currentScore)
		elif hit.name[0:6] == 'POL-RH':
			classify_RH(hitList[i:], classifs, hit.score + currentScore)
		elif hit.name[0:7] == 'POL-INT':
			classify_INT(hitList[i:], classifs, hit.score + currentScore)

def classify_RT(hitList, classifs, currentScore):
	for hit, i in zip(hitList[1:], range(1, len(hitList))):
		if hit.start-hitList[0].end > 2000:
			continue
		if hit.name[0:6] == 'POL-RH':
			classify_RH(hitList[i:], classifs, hit.score + currentScore)
		elif hit.name[0:7] == 'POL-INT':
			classifs.append(('Gypsy', hit.score + currentScore))

def classify_RH(hitList, classifs, currentScore):
	for hit in hitList[1:]:
		if hit.start-hitList[0].end > 2000:
			continue
		elif hit.name[0:7] == 'POL-INT':
			classifs.append(('Gypsy', hit.score + currentScore))

def classify_INT(hitList, classifs, currentScore):
	for hit, i in zip(hitList[1:], range(1, len(hitList))):
		if hit.start-hitList[0].end > 2000:
			continue
		elif hit.name[0:6] == 'POL-RT':
			classify_INT_RT(hitList[i:], classifs, hit.score + currentScore)
		elif hit.name[0:6] == 'POL-RH':
			classifs.append(('Copia', hit.score + currentScore))

def classify_INT_RT(hitList, classifs, currentScore):
	for hit in hitList[1:]:
		if hit.start-hitList[0].end > 2000:
			continue
		elif hit.name[0:6] == 'POL-RH':
			classifs.append(('Copia', hit.score + currentScore))
	classifs.append(('Copia', currentScore))

def clean():
	os.remove(orfFile)
	os.remove(hmmHitFile)
	if args.use_blast_evidence:
		os.remove(blastHitFile)

if __name__ == '__main__':
	#setup global & static variables
	timestamp = timestamp.timestamp()
	args = parseArgs()
	blastDB = '/share/jyllwgrp/nealedata/databases/cprd.fa'
	hmmProfilesForLTRs = '/share/jyllwgrp/nealedata/databases/ltr-Pfam.hmm'
	targetSeqs = '/share/jyllwgrp/nealedata/databases/pier-1.3.fa'
	orfFile = args.fasta + '.orfs'
	hmmHitFile = 'hmmer-output-{}.txt'.format(timestamp)
	blastHitFile = 'blast-output-{}.txt'.format(timestamp)
	outputFile = open(args.output_fasta_file, 'w')
	main()
	
