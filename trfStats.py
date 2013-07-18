#!/usr/bin/python
#Author: Jacob Zieve
import sys,argparse,textwrap,os
from pickler import sendToPickleJar, getFromPickleJar
from classes import trfHit
from collections import defaultdict
from operator import itemgetter, attrgetter

def main():
	args = parseArgs()
	outFileBase = args.trf_output
	trfDict = getTRFDict(args)
	pickleTRF(args,trfDict)
	getStats(args,trfDict,outFileBase)
		
def parseArgs():
	myDescription = '''
	Processes TRF output. It does nothing if no parameters are supplied.

	Example usage (if pathed):
		trfStats.py trf_output.dat -f -m
	Sample line from TRF .dat file:
	start end period_length copy_num consensus_length %matches %indels score %a %c %g %t entropy period sequence 
	33493 33544 22 2.3 22 83 3 59 34 3 1 59 1.27 ATTTTTAAAACATTTTTAATTC ATTTTAAAAATATTTTTAATTTATTTTTAAAACATTTTTTATTGCATTTTTA
	'''

	argParser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(myDescription))
	argParser.add_argument('trf_output')
	argParser.add_argument('-t', '--tab_delim', help='Prints trf_output as a tab delimited file', action='store_true')
	argParser.add_argument('-ngs', '--next_gen_seq', help='Input file was generated with TRF(v4.07+) option -ngs, a .ngs file', action='store_true')
	argParser.add_argument('-pis', '--send_to_pickle', nargs =1, help='Send the main dictionary to a pickle file')
	argParser.add_argument('-pil', '--load_from_pickle', nargs=1, help='Load the main dictionary to a pickle file')
	argParser.add_argument('-g', '--greedy_filter', help='Filter selects monomers over multimers if overlapping (i.e. http://en.wikipedia.org/wiki/Activity_selection_problem)', action='store_true')
	argParser.add_argument('-f', '--filter', help='Filter selects multimers over monomers if overlapping (i.e. maximizes the highest scoring loci' , action='store_true')
	argParser.add_argument('-p', '--period_frequencies', help='Outputs a tab-delimited .period_stats file with number of occurences of a period as a tab-delimeted file with the first field being the period and the second being the frequecy', action='store_true')
	argParser.add_argument('-c', '--period_cumulative_length', help='Outputs a tab-delimited .cumulative_stats file with the first field being the period and the second being the cumulative bp', action='store_true')
	argParser.add_argument('-m', '--motif_frequencies', help='Outputs tab-delimited .motif_stats file with the first field being the motif and the second being the frequency', action='store_true')
	argParser.add_argument('-mc', '--motif_cumulative_length', help='Outputs tab-delimited .motif_bp file with the first field being the motif and the second being the frequency', action='store_true')
	argParser.add_argument('-fa', '--fasta', help='Outputs a fasta file with all the tandem repeats for further processing', action='store_true')
	
	return argParser.parse_args()

def getTRFDict(args):
	if args.load_from_pickle:
		return getFromPickleJar(args.load_from_pickle[0])
	elif args.next_gen_seq:
		return parseNGS(args)
	else:
		return parseDat(args)


def getStats(args,trfDict,outFileBase):
	trfList = flattenTRFDict(trfDict)
	if args.tab_delim:
		writeTRFDict(trfList,outFileBase)
	if args.period_frequencies:
		 writePeriodFreqs(trfList,outFileBase)
	if args.period_cumulative_length:
		 writePeriodBp(trfList,outFileBase)
	if args.motif_frequencies:
		 writeMotifFreqs(trfList,outFileBase)
	if args.motif_cumulative_length:
		writeMotifLengths(trfList,outFileBase)
	if args.fasta:
		 writeFasta(trfList,outFileBase)	 	


def pickleTRF(args,trfDict):
	if args.send_to_pickle:
		sendToPickleJar(trfDict,args.send_to_pickle[0])
	
def parseNGS(args):
	trfDict = defaultdict(list)
	seq_name = ''
	for line in open(args.trf_output).read().split('\n'):	
		if line.startswith('@'):
			seq_name = line[1:]
		elif line != '':
			fields = line.split()
			fields.insert(0,seq_name)
			newTandem = trfHit(fields)
			print newTandem._dict_.keys()
			if args.filter:
				trfDict = selectMultimer(newTandem,trfDict)
			else:
				trfDict[seq_name].append(newTandem)
		else:
			continue
	return trfDict

def overlap(t1,t2):
	if ( (t1.start >= t2.end) or (t2.start >= t1.end) ):
                return False
        else:   
                return True
	

def selectMultimer(newTandem,trfDict):
	key = newTandem.seq_name
	if trfDict.has_key(key):
		for tandem in trfDict[key]:
			if overlap(newTandem,tandem):
				if newTandem.score > tandem.score:
					trfDict[key].remove(tandem)
					trfDict[key].append(newTandem)
	else:
		trfDict[key].append(newTandem)
	return trfDict

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
			if args.filter:
				trfDict = selectMultimer(newTandem,trfDict)
			else:
				trfDict[seq_name].append(newTandem)
	return trfDict

def flattenTRFDict(trfDict):
	trfList = []
	for v in trfDict.values():
		trfList.extend(v)
	return trfList

#def selectMonomers(trfDict):
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
def writeTRFDict(trfList,outFileBase):
	out = open(outFileBase+'.txt','w')
	header = "seq_name\tstart\tend\tperiod\tcopies\tconsensus\t%matches\t%indels\tscore\t%a\t%c\t%g\t%t\tentropy\tlength\tmotif\n"
	out.write(header)
	for t in trfList:
		out.write(str(t));

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

def writeCumulativeBpByPeriodRange(trfList,outFileBase,prange):
	out = open(outFileBase+'.cumulative_stats','w')
	out.write('Between {:d} and {:d}:\t{:d}\n'.format(prange[0],prange[1],getCumulativeBpByPeriodRange(trfList,(prange[0],prange[1]))))
	out.close()

def writePeriodBp(trfList,outFileBase):
	writeCumulativeBpByPeriodRange(trfList,outFileBase,(2,3))
	writeCumulativeBpByPeriodRange(trfList,outFileBase,(3,4))
	writeCumulativeBpByPeriodRange(trfList,outFileBase,(5,6))
	writeCumulativeBpByPeriodRange(trfList,outFileBase,(7,8))
	writeCumulativeBpByPeriodRange(trfList,outFileBase,(8,9))
	writeCumulativeBpByPeriodRange(trfList,outFileBase,(10,31))
	writeCumulativeBpByPeriodRange(trfList,outFileBase,(30,51))
	writeCumulativeBpByPeriodRange(trfList,outFileBase,(51,71))
	writeCumulativeBpByPeriodRange(trfList,outFileBase,(70,101))
	writeCumulativeBpByPeriodRange(trfList,outFileBase,(100,201))
	writeCumulativeBpByPeriodRange(trfList,outFileBase,(200,301))
	writeCumulativeBpByPeriodRange(trfList,outFileBase,(300,401))
	writeCumulativeBpByPeriodRange(trfList,outFileBase,(400,501))
	writeCumulativeBpByPeriodRange(trfList,outFileBase,(1,501))

def getMotifFreqs(trfList):
	motifDict = defaultdict(int)
	for t in trfList:
		motifDict[t.motif]+=1
	return motifDict
	
def writeMotifFreqs(trfList,outFileBase):
	out = open(outFileBase+'.motif_freqs','w')
	for k,v in sorted(getMotifFreqs(trfList).items()):
		out.write('{:s}\t{:d}\n'.format(k,v))
	out.close()

def getCumulativeMotifLengths(trfList):
	motifs = defaultdict(int)
	for t in trfList:
		motifs[t.motif] += t.length
	return motifs

def writeMotifLengths(trfList,outFileBase):
	out = open(outFileBase+'.motif_lengths','w')
	for k,v in sorted(getCumulativeMotifLengths(trfList).items()):
		out.write('{:s}\t{:d}\n'.format(k,v))
	out.close()
	
def getCumulativeBpByPeriodRange(trfList,prange):
	return sum([t.length for t in filter(lambda x: x.period in range(prange[0],prange[1]),trfList)])
	
 
def writeFasta(trfList,outFileBase):
	out = open(outFileBase+'.fa','w')
	for t in trfList:
		header = '>'+t.seq_name+'|tandem repeat|'+str(t.start)+':'+str(t.end)+'|length='+str(t.length)+'|period='+str(t.period)
		out.write(header+'\n'+t.sequence+'\n')
	out.close()


if __name__=='__main__':
	main()
