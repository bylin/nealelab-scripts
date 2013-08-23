#!/usr/bin/python
#Author: Jacob Zieve
import sys,argparse,textwrap,os,re,pdb
from pickler import sendToPickleJar, getFromPickleJar
from classes import trfHit
from collections import defaultdict
from operator import itemgetter, attrgetter

def main():
	args = parseArgs()
	if args.filter_monomers:
		outFileBase = args.trf_output+".multimers"
	elif args.filter_multimers:
		outFileBase = args.trf_output+".monomers"
	else:
		outFileBase = args.trf_output
		
	trfDict = getTRFDict(args)
	getStats(args,trfDict,outFileBase)
		
def parseArgs():
	myDescription = '''
	Processes TRF output. It does nothing if no parameters are supplied.

	Example usage (if pathed):
		trfStats.py -s trf_output.dat
	Sample line from TRF .dat file:
	start end period_length copy_num consensus_length %matches %indels score %a %c %g %t entropy period sequence 
	33493 33544 22 2.3 22 83 3 59 34 3 1 59 1.27 ATTTTTAAAACATTTTTAATTC ATTTTAAAAATATTTTTAATTTATTTTTAAAACATTTTTTATTGCATTTTTA
	'''

	argParser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(myDescription))
	argParser.add_argument('trf_output')
	argParser.add_argument('-ngs', '--next_gen_seq', help='Input file was generated with TRF(v4.07+) option -ngs, a .ngs file', action='store_true')
	argParser.add_argument('-csv', '--csv', help='Prints trf_output as a csv', action='store_true')
	argParser.add_argument('-sql', '--sql', help='Prints trf_output as a MySQL dump for further querying', action='store_true')
	###Deprecated###
	#argParser.add_argument('-pis', '--send_to_pickle', nargs =1, help='Send the main dictionary to a pickle file')
	#argParser.add_argument('-pil', '--load_from_pickle', nargs=1, help='Load the main dictionary to a pickle file')
	#argParser.add_argument('-g', '--greedy_filter', help='Filter selects monomers over multimers if overlapping (i.e. http://en.wikipedia.org/wiki/Activity_selection_problem)', action='store_true')
	argParser.add_argument('-fmo', '--filter_monomers', help='Filter selects multimers over monomers if overlapping (i.e. selects the highest scoring loci), WARNING: still some overlap' , action='store_true')
	argParser.add_argument('-fmu', '--filter_multimers', help='Should be default to remover redundancy, it filters selects monomers over multimers if overlapping (i.e. selects the lowest scoring loci)' , action='store_true')
	argParser.add_argument('-pf', '--period_frequencies', help='Outputs a tab-delimited .period_freqs file with number of occurences of a period as a tab-delimeted file with the first field being the period and the second being the frequecy', action='store_true')
	argParser.add_argument('-pl', '--period_length', help='Outputs a tab-delimited .period_lengths file with the first field being the period and the second being the cumulative bp', action='store_true')
	argParser.add_argument('-pt', '--period_table', help='Outputs a tab-delimited .period_tbl file with the period and total number of loci,copies,variants and cumulative bp', action='store_true')
	argParser.add_argument('-mf', '--motif_frequencies', help='Outputs tab-delimited .motif_freqs file with the first field being the motif and the second being the frequency', action='store_true')
	argParser.add_argument('-ml', '--motif_length', help='Outputs tab-delimited .motif_lengths file with the first field being the motif and the second being the cumulative bp', action='store_true')
	argParser.add_argument('-fa', '--fasta', help='Outputs a fasta file with all the tandem repeats for further processing', action='store_true')
	argParser.add_argument('-tel', '--telomeres', help='Asseses potential telomeres and interstitial telomeric sequences', action='store_true')
	argParser.add_argument('-s', '--stats', help = 'Performs all the above output operations', action='store_true')
	
	return argParser.parse_args()

def getTRFDict(args):
	if args.next_gen_seq:
		return parseNGS(args)
	else:
		return parseDat(args)


def getStats(args,trfDict,outFileBase):
	trfList = flattenTRFDict(trfDict)
	if args.csv:
		writeCSV(trfList,outFileBase)
	if args.sql:
		writeSQL(trfList,outFileBase)
	if args.period_frequencies:
		 writePeriodFreqs(trfList,outFileBase)
	if args.period_length:
		writePeriodLengths(trfList,outFileBase)
	if args.period_table:
		writePeriodTable(trfList,outFileBase)
	if args.motif_frequencies:
		 writeMotifFreqs(trfList,outFileBase)
	if args.motif_length:
		writeMotifLengths(trfList,outFileBase)
	if args.fasta:
		writeFasta(trfList,outFileBase)	 	
	if args.telomeres:
		writeTelomeres(trfList,outFileBase)
	if args.stats:
		writeCSV(trfList,outFileBase)
		writeSQL(trfList,outFileBase)
		writePeriodFreqs(trfList,outFileBase)
		writePeriodLengths(trfList,outFileBase)
		writePeriodTable(trfList,outFileBase)
		writeMotifFreqs(trfList,outFileBase)
		writeMotifLengths(trfList,outFileBase)
		writeFasta(trfList,outFileBase)	 	
		writeTelomeres(trfList,outFileBase)
	
def parseNGS(args):
	trfDict = defaultdict(list)
	seq_name = ''
	count = 0
	for line in open(args.trf_output).read().split('\n'):	
		if line.startswith('@'):
			seq_name = line[1:]
		elif line != '':
			fields = line.split()
			fields.insert(0,seq_name)
			newTandem = trfHit(fields)
			if args.filter_monomers:
				filterTRFDict(newTandem,trfDict,"monomers")
			elif args.filter_multimers:
				filterTRFDict(newTandem,trfDict,"multimers")
			else:
				trfDict[seq_name].append(newTandem)
	return trfDict

def overlap(t1,t2):
	if ( (t1.start >= t2.end) or (t2.start >= t1.end) ):
                return False
        else:   
                return True
	
def filterTRFDict(newTandem,trfDict,discard):
	key = newTandem.seq_name
	overlapFlag = False
	if trfDict.has_key(key):
		for i,tandem in enumerate(trfDict[key][:]):
			if overlap(newTandem,tandem):
				overlapFlag = True
				if (discard == "monomers") and (newTandem.score >= tandem.score):
					trfDict[key].pop(i)
					trfDict[key].append(newTandem)
			 	if (discard == "multimers") and (newTandem.score <= tandem.score):
					trfDict[key].pop(i)
					trfDict[key].append(newTandem)

	if(not overlapFlag):
		trfDict[key].append(newTandem)

def parseDat(args):
	trfDict = defaultdict(list)
	seq_name = ''
	for line in open(args.trf_output).read().split('\n'):	
		if line.startswith('Seq'):
			seq_name = line[len('Sequence: '):]
		if line != '' and line.split()[0].isdigit():
			fields = line.split()
			fields.insert(0,seq_name)
			newTandem = trfHit(fields)
			if args.filter_monomers:
				filterTRFDict(newTandem,trfDict,"monomers")
			elif args.filter_multimers:
				filterTRFDict(newTandem,trfDict,"multimers")
			else:
				trfDict[seq_name].append(newTandem)
	return trfDict

def flattenTRFDict(trfDict):
	trfList = []
	for v in trfDict.values():
		trfList.extend(v)
	return trfList

#greedy algorithm, can be faster
#def selectMonomers(trfDict):#greedy
#	trf = []
#	for seq_name in trfDict.keys():
#		seqs=sorted(trfDict[seq_name],key=attrgetter('end'))
#		trf.append(seqs[0])
#		f = seqs[0].end
#		for i in range(1,len(seqs)):
#			if seqs[i].start >= f:
#				trf.append(seqs[i])
#				f = seqs[i].end
#	return trf
#

def writeCSV(trfList,outFileBase):#mostly for visualization in R script
	out = open(outFileBase+'.csv','w')
	motifFreqs = getMotifFreqs(trfList)
	periodFreqs = getPeriodFreqs(trfList)
	header = "seq_name,start,end,period,copies,consensus,%matches,%indels,score,%a,%c,%g,%t,entropy,motif,length,motif_freq,period_freq\n"
	out.write(header)
	for t in trfList:
		motifFreq = motifFreqs[t.motif]
		periodFreq = periodFreqs[t.period]
		out.write(str(t)+','+str(motifFreq)+','+str(periodFreq)+'\n');

def writeSQL(trfList,outFileBase):
	table = outFileBase[:outFileBase.find('.')]+'_trf'#tables can't have . in them, quick hack
	out = open(table+'.dmp','w')
	motifFreqs = getMotifFreqs(trfList)
	periodFreqs = getPeriodFreqs(trfList)
	out.write('BEGIN;\n')
	out.write('DROP TABLE IF EXISTS `'+table+'`;\n')
        out.write('CREATE TABLE '+table+'(id INT NOT NULL AUTO_INCREMENT PRIMARY KEY, seq_name VARCHAR(200), start INT, end INT, period INT, copies FLOAT, consensus INT, matches FLOAT,\
		indels FLOAT, score INT, A FLOAT, C FLOAT, G FLOAT, T FLOAT, entropy FLOAT, motif VARCHAR(2001),len INT, motif_freq INT, period_freq INT);\n')
	for t in trfList:
		motifFreq = motifFreqs[t.motif]
		periodFreq = periodFreqs[t.period]
		out.write('INSERT INTO '+table+'(seq_name,start,end,period,copies,consensus,matches,indels,score,A,C,G,T,entropy,motif,len,motif_freq,period_freq) VALUES('+str(t)+','+str(motifFreq)+','+str(periodFreq)+');\n')
	out.write('COMMIT;')

def getPeriodFreqs(trfList):
	periodDict = defaultdict(int)
	for t in trfList:
		periodDict[t.period]+=1
	return periodDict

def writePeriodFreqs(trfList,outFileBase):
	out = open(outFileBase+'.period_freqs','w')
	for k,v in sorted(getPeriodFreqs(trfList).items()):
		out.write('{:d}\t{:d}\n'.format(k,v))
	out.close()

def getCumulativeBpByPeriodRange(trfList,prange):
	return sum([t.length for t in filter(lambda x: x.period in range(prange[0],prange[1]),trfList)])


def getCumulativePeriodLengths(trfList):
	periods = defaultdict(int)
	for t in trfList:
		periods[t.period] += t.length
	return periods

def writePeriodLengths(trfList,outFileBase):
	out = open(outFileBase+'.period_lengths','w')
	for k,v in sorted(getCumulativePeriodLengths(trfList).items()):
		out.write('{:d}\t{:d}\n'.format(k,v))
	out.close()

def getMaxPeriod(trfList):
	return max([t.period for t in trfList])

'''def writePeriodRange(trfList,outFileBase):
	with open(outFileBase+'.period_range.txt','w') as out:
		out = writeCumulativeBpByPeriodRange(trfList,out,(2,3))
		out = writeCumulativeBpByPeriodRange(trfList,out,(3,4))
		out = writeCumulativeBpByPeriodRange(trfList,out,(4,5))
		out = writeCumulativeBpByPeriodRange(trfList,out,(5,6))
		out = writeCumulativeBpByPeriodRange(trfList,out,(6,7))
		out = writeCumulativeBpByPeriodRange(trfList,out,(7,8))
		out = writeCumulativeBpByPeriodRange(trfList,out,(8,9))#micros
		out = writeCumulativeBpByPeriodRange(trfList,out,(9,101))#minis
		out = writeCumulativeBpByPeriodRange(trfList,out,(101,getMaxPeriod(trfList)))#sats
'''

def getLociByPeriodRange(trfList,prange):
	return sum([1 for t in filter(lambda x:x.period in range(prange[0],prange[1]),trfList)])

def getCopiesByPeriodRange(trfList,prange):
	return sum([t.copy_num for t in filter(lambda x:x.period in range(prange[0],prange[1]),trfList)])

def getVariantsByPeriodRange(trfList,prange):
	return len(set([t.motif for t in filter(lambda x:x.period in range(prange[0],prange[1]),trfList)]))

def getCumulativeBpByPeriodRange(trfList,prange):
	return sum([t.length for t in filter(lambda x:x.period in range(prange[0],prange[1]),trfList)])

def writePeriodTable(trfList,outFileBase):
	with open(outFileBase+'.period_tbl','w') as out:
		maxPeriod = getMaxPeriod(trfList)
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('Dinucleotide',getLociByPeriodRange(trfList,(2,3)),getCopiesByPeriodRange(trfList,(2,3)),getVariantsByPeriodRange(trfList,(2,3)),getCumulativeBpByPeriodRange(trfList,(2,3))))
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('Trinucleotide',getLociByPeriodRange(trfList,(3,4)),getCopiesByPeriodRange(trfList,(3,4)),getVariantsByPeriodRange(trfList,(3,4)),getCumulativeBpByPeriodRange(trfList,(3,4))))
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('Tetranucleotide',getLociByPeriodRange(trfList,(4,5)),getCopiesByPeriodRange(trfList,(4,5)),getVariantsByPeriodRange(trfList,(4,5)),getCumulativeBpByPeriodRange(trfList,(4,5))))
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('Pentanucleotide',getLociByPeriodRange(trfList,(5,6)),getCopiesByPeriodRange(trfList,(5,6)),getVariantsByPeriodRange(trfList,(5,6)),getCumulativeBpByPeriodRange(trfList,(5,6))))
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('Hexanucleotide',getLociByPeriodRange(trfList,(6,7)),getCopiesByPeriodRange(trfList,(6,7)),getVariantsByPeriodRange(trfList,(6,7)),getCumulativeBpByPeriodRange(trfList,(6,7))))
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('Heptanucleotide',getLociByPeriodRange(trfList,(7,8)),getCopiesByPeriodRange(trfList,(7,8)),getVariantsByPeriodRange(trfList,(7,8)),getCumulativeBpByPeriodRange(trfList,(7,8))))
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('Octanucleotide',getLociByPeriodRange(trfList,(8,9)),getCopiesByPeriodRange(trfList,(8,9)),getVariantsByPeriodRange(trfList,(8,9)),getCumulativeBpByPeriodRange(trfList,(8,9))))
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('9-30 bp',getLociByPeriodRange(trfList,(9,31)),getCopiesByPeriodRange(trfList,(9,31)),getVariantsByPeriodRange(trfList,(9,31)),getCumulativeBpByPeriodRange(trfList,(9,31))))
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('31-50 bp',getLociByPeriodRange(trfList,(31,51)),getCopiesByPeriodRange(trfList,(31,51)),getVariantsByPeriodRange(trfList,(31,51)),getCumulativeBpByPeriodRange(trfList,(31,51))))
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('51-70 bp',getLociByPeriodRange(trfList,(51,71)),getCopiesByPeriodRange(trfList,(51,71)),getVariantsByPeriodRange(trfList,(51,71)),getCumulativeBpByPeriodRange(trfList,(51,71))))
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('71-100 bp',getLociByPeriodRange(trfList,(71,101)),getCopiesByPeriodRange(trfList,(71,101)),getVariantsByPeriodRange(trfList,(71,101)),getCumulativeBpByPeriodRange(trfList,(71,101))))
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('101-200 bp',getLociByPeriodRange(trfList,(101,201)),getCopiesByPeriodRange(trfList,(101,201)),getVariantsByPeriodRange(trfList,(101,201)),getCumulativeBpByPeriodRange(trfList,(101,201))))
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('201-300 bp',getLociByPeriodRange(trfList,(201,301)),getCopiesByPeriodRange(trfList,(201,301)),getVariantsByPeriodRange(trfList,(201,301)),getCumulativeBpByPeriodRange(trfList,(201,301))))
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('301-400 bp',getLociByPeriodRange(trfList,(301,401)),getCopiesByPeriodRange(trfList,(301,401)),getVariantsByPeriodRange(trfList,(301,401)),getCumulativeBpByPeriodRange(trfList,(301,401))))
		out.write('{:s}\t{:d}\t{:.1f}\t{:d}\t{:d}\n'.format('>400 bp',getLociByPeriodRange(trfList,(401,maxPeriod)),getCopiesByPeriodRange(trfList,(401,maxPeriod)),getVariantsByPeriodRange(trfList,(401,maxPeriod)),getCumulativeBpByPeriodRange(trfList,(401,maxPeriod))))
		
def getMotifFreqs(trfList):
	motifDict = defaultdict(int)
	for t in trfList:
		motifDict[t.motif]+=1
	return motifDict
	
def writeMotifFreqs(trfList,outFileBase):
	out = open(outFileBase+'.motif_freqs','w')
	for k,v in sorted(getMotifFreqs(trfList).items(),key=itemgetter(1),reverse=True):
		out.write('{:s}\t{:d}\n'.format(k,v))
	out.close()

def getCumulativeMotifLengths(trfList):
	motifs = defaultdict(int)
	for t in trfList:
		motifs[t.motif] += t.length
	return motifs

def writeMotifLengths(trfList,outFileBase):
	out = open(outFileBase+'.motif_lengths','w')
	for k,v in sorted(getCumulativeMotifLengths(trfList).items(),key=itemgetter(1),reverse=True):
		out.write('{:s}\t{:d}\n'.format(k,v))
	out.close()
	
 
def writeFasta(trfList,outFileBase):
	out = open(outFileBase+'.fa','w')
	for t in trfList:
		header = '>'+t.seq_name+'|tandem repeat|'+str(t.start)+':'+str(t.end)+'|length='+str(t.length)+'|period='+str(t.period)
		out.write(header+'\n'+t.sequence+'\n')
	out.close()

def isTelomere(tandem):
	telomere = re.search('(TTTAGGG){3,}',tandem.sequence)
	if telomere is not None and tandem.length > 1000:
		return True
	else:
		return False
	
def writeTelomeres(trfList,outFileBase):
	out = open(outFileBase+'.telomeres','w')
	for t in trfList:
		if isTelomere(t):
			out.write(str(t)+"\n")

if __name__=='__main__':
	main()
