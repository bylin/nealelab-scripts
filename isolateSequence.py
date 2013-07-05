#!/usr/bin/python

from Bio import SeqIO
import argparse

def parseArgs():
	argParser = argparse.ArgumentParser(description='Isolate sequences from a fasta file')
	argParser.add_argument('-in', '--input_fasta_file', required=True, help='Input fasta file')
	argParser.add_argument('-list', '--sequence_list', required=True, help='File containing sequence headers for each sequence to be isolated')
	argParser.add_argument('-one', '--one_sequence_per_file', help='Create new file for each sequence, NOT FUNCTIONAL')
	argParser.add_argument('-out', '--output_fasta_file', required=True, help='Output fasta file')
	args = argParser.parse_args()
	args.sequence_list = [seq for seq in open(args.sequence_list)]
	args.output_fasta_file = open(args.output_fasta_file, 'w')
	
def main():
	isolateIntoOneFile()
	if seqsAreMissing():
		print 'Unisolated sequences:' + ', '.join(args.sequence_list)

def isolateIntoOneFile():
	for seq in SeqIO.parse(args.input_fasta_file, 'fasta'):
		if seq.description in args.sequence_list:
			SeqIO.write(seq, args.output_fasta_file, 'fasta')
			args.sequence_list.remove(seq.description)

def seqsAreMissing():
	if len(args.sequence_list) > 0:
		return True
	return False

if __name__ == '__main__':
	args = parseArgs()
	main()
