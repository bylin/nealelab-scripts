#!/usr/bin/python
# Author: Brian Lin
# SGE wrapper. Generates a helper file and runs it.

import subprocess, argparse, os, textwrap, datetime
global timestamp
timestamp = datetime.datetime.today().strftime("%y%m%d%H%m%S")

def main():
	args = parseArgs()
	os.mkdir(timestamp)
	generateQsubScript(args)
	cmd = 'qsub ' + timestamp + '/SGEWrapper.sh'
	subprocess.call(cmd, shell=True)

def parseArgs():
	mydesc = '''
	SGE wrapper. Contains options for qsub's built-in job arrays.
	
	Example usage: 
		SGEJobWrapper.py -m gcc repeatmasker -cmd 'RepeatMasker input_file.fa -pa 10 -lib library_file.fa -q -nolow -gff'
	or, if you want to run a command on multiple files in an input directory using SGE job arrays, 15 slots:
		SGEJobWrapper.py -m gcc repeatmasker -j -s 15 -d [input_directory] -cmd 'RepeatMasker $INPUT_FILE -pa 10 -lib library_file.fa -q -nolow -gff '''
	argParser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(mydesc))
	argParser.add_argument('-m', '--modules', nargs='+', required=True, help='Modules to load')
	argParser.add_argument('-cmd', '--command', nargs=1, required=True, help='Command to qsub')
	argParser.add_argument('-o', '--stdout', nargs=1, help='Standard output file', default='SGEWrapper.stdout')
	argParser.add_argument('-e', '--stderr', nargs=1, help='Standard error file', default='SGEWrapper.stderr')
	argParser.add_argument('-j', '--use_job_array', help='Use a SGE job array', action='store_true')
	argParser.add_argument('-s', '--slots', nargs=1, type=int, help='Number of slots to use with the job array')
	argParser.add_argument('-d', '--input_directory', nargs=1, help='Directory containing all input files to be scheduled with the command in the job array')
	return argParser.parse_args()

def generateQsubScript(args):
	if args.use_job_array:
		generateQsubArrayScript(args)
	else:
		generateQsubSingleJobScript(args)

def generateQsubArrayScript(args):
	outfile = file(timestamp + '/SGEWrapper.sh', 'w')
	inputDirectory = args.input_directory[0]
	nJobParser = lambda d: str(len(os.listdir(inputDirectory)))
	nJobs = nJobParser(inputDirectory)
	script = '#!/bin/bash\n# Qsub script generated from SGEWrapper.py\n#$ -S /bin/bash\n#$ -cwd\n#$ -N SGEWrapper\n#$ -q bigmem1.q\n'
	script += '#$ -t 1-' + nJobs + '\n'
	script += '#$ -o ' + timestamp + '/' + args.stdout + '\n'
	script += '#$ -e ' + timestamp + '/' + args.stderr + '\n\n'
	script += 'module load ' + ' '.join(args.modules) + '\n\n'
	script += 'FILES=(`ls ' + inputDirectory +'/*`)\n'
	script += 'INPUT_FILE=${FILES[SGE_TASK_ID-1]}\n'
	script += args.command[0] + '\n'
	outfile.write(script)

def generateQsubSingleJobScript(args):
	outfile = file('SGEJobWrapper.sh', 'w')
	script = '''#!/bin/bash
	#$ -S /bin/bash
	#$ -cwd
	#$ -N SGEJobWrapper
	#$ -q bigmem1.q
	#$ -o ''' + args.stdout + '''
	#$ -e ''' + args.stderr + '\n'
	for module in args.modules:
		script += 'module load ' + module + '\n'
	script += args.command[0] + '\n'
	outfile.write(script)

if __name__ == '__main__':
	main()

#RepeatMasker [input] -pa 10 -lib pier-1.3.fa (-q) -nolow -gff
