#!/usr/bin/python
# Author: Brian Lin

import argparse, subprocess, timestamp, os

def parseArgs():
	parser = argparse.ArgumentParser(description="Translated nucleotide HMMER search, default settings. Uses the findorfs utility from USEARCH to find six-frame translations from nucleotide FASTA input")
	parser.add_argument('-hmm', '--hmm_profiles', help='Database HMM-formatted binary file')
	parser.add_argument('-i', '--input_nucl_fasta', help='Input .fasta file')
	parser.add_argument('-o', '--output_file', help='Output hmmer results')
	return parser.parse_args()

def main():
	script = generateScript()
	cmd = 'bash {}'.format(script)
	subprocess.call(cmd, shell=True)
	os.remove(script)

def generateScript():
	filename = timestamp.timestamp() + '.hmmern.sh'
	moduleScript = open(filename, 'w')
	loadModules(moduleScript)
	generateORFs(moduleScript)
	hmmer(moduleScript)
	return filename

def loadModules(moduleScript):
	moduleScript.write('module load usearch hmmer\n')

def generateORFs(moduleScript):
	orfs = args.input_nucl_fasta + '.orfs'
	if not os.path.exists(orfs):
		moduleScript.write('echo Generating ORFs ...\n')
		moduleScript.write('usearch -findorfs {} -output {} -xlat &> /dev/null\n'.format(args.input_nucl_fasta, orfs))
	moduleScript.write('echo Got ORFs\n')

def hmmer(moduleScript):
	orfs = args.input_nucl_fasta + '.orfs'
	moduleScript.write('echo Running hmmsearch ...\n')
	moduleScript.write('hmmsearch -o /dev/null --domtblout {} {} {}\n'.format(args.output_file, args.hmm_profiles, orfs))
	moduleScript.write('echo Finished hmmsearch\n')
	moduleScript.close()

if __name__ == "__main__":
	args = parseArgs()
	main()
