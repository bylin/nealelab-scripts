#!/usr/bin/python
# Author: Brian Lin
# SGE wrapper. Generates a helper file and runs it.

import subprocess, argparse, os, textwrap, datetime, shutil
global timestamp
timestamp = 'run-'+datetime.datetime.today().strftime("%d-%m-%y-%H%m")

def main():
	args = parseArgs()
	setUpEnv(args)
	generateQsubScript(args)
	cmd = 'qsub SGEWrapper.sh'
	subprocess.call(cmd, shell=True)

def parseArgs():
	mydesc = '''
	SGE wrapper. Contains options for qsub's built-in job arrays.
	
	Example usage: 
		SGEWrapper.py -m gcc repeatmasker -cmd 'RepeatMasker input_file.fa -pa 10 -lib library_file.fa -q -nolow -gff'
	
	or, if you want to run a command on multiple files in an input directory using SGE job arrays, 5 slots:
		SGEWrapper.py -j -tc 5 -d [input_directory] -m repeatmasker -f library_file.fa -cmd 'RepeatMasker $INPUT_FILE -lib library_file.fa -q -nolow -gff' '''
	argParser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(mydesc))
	argParser.add_argument('-j', '--use_job_array', help='Use a SGE job array', action='store_true')
	argParser.add_argument('-tc', '--max_running_tasks', nargs=1, help='Number of concurrent jobs to run at a time')
	argParser.add_argument('-m', '--modules', nargs='+', required=True, help='Modules to load')
	argParser.add_argument('-d', '--input_directory', nargs=1, help='Directory containing all input files to be scheduled with the command in the job array')
	argParser.add_argument('-f', '--extra_files', nargs='+', help='Extra files that need to be used by the program')
	argParser.add_argument('-o', '--stdout', nargs=1, help='Standard output file', default=['SGEWrapper.stdout'])
	argParser.add_argument('-e', '--stderr', nargs=1, help='Standard error file', default=['SGEWrapper.stderr'])
	argParser.add_argument('-M', '--email', nargs=1, help='Email account user wishes to send job start and end notifications to')
	argParser.add_argument('-cmd', '--command', nargs=1, required=True, help='Command to qsub')
	return argParser.parse_args()

def setUpEnv(args):
	if os.path.exists(timestamp): shutil.rmtree(timestamp)
	os.mkdir(timestamp)
	if args.extra_files is not None: 
		for extra_file in args.extra_files:
			shutil.copy(extra_file, timestamp)
	os.chdir(timestamp)

def generateQsubScript(args):
	if args.use_job_array:
		generateQsubArrayScript(args)
	else:
		generateQsubSingleJobScript(args)

def generateQsubArrayScript(args):
	outfile = file('SGEWrapper.sh', 'w')
	inputDirectory = '../' + args.input_directory[0]
	(inputFile, nJobs) = scanInputs(inputDirectory)
	script = '#!/bin/bash\n# Qsub script generated from SGEWrapper.py\n#$ -S /bin/bash\n#$ -cwd\n#$ -N SGEWrapper\n#$ -q bigmem1.q\n'
	script += '#$ -t 1-' + nJobs + '\n'
	script += '#$ -tc ' + args.max_running_tasks[0] + '\n\n'
	script += '#$ -o ' + args.stdout[0] + '\n'
	script += '#$ -e ' + args.stderr[0] + '\n\n'
	script += 'module load ' + ' '.join(args.modules) + '\n\n'
	script += 'INPUT_FILE=`sed -n "${SGE_TASK_ID}p" ' + inputFile + '`\n'
	script += args.command[0] + '\n'
	outfile.write(script)

def scanInputs(inputDirectory):
	inputFile = 'SGEWrapperInputFiles.txt'
	cmd = 'ls ' + inputDirectory + '/* > ' + inputFile
	proc = subprocess.call(cmd, shell=True)
	lineCounter = 0
	inputFileHandle = open(inputFile)
	nJobs = sum(lineCounter + 1 for line in inputFileHandle)
	return(inputFile, str(nJobs))

def generateQsubSingleJobScript(args):
	outfile = file('SGEJobWrapper.sh', 'w')
	script = '''#!/bin/bash
	#$ -S /bin/bash
	#$ -cwd
	#$ -N SGEJobWrapper
	#$ -q bigmem1.q
	#$ -M ''' + args.email + '''
	#$ -o ''' + args.stdout + '''
	#$ -e ''' + args.stderr + '\n'
	for module in args.modules:
		script += 'module load ' + module + '\n'
	script += args.command[0] + '\n'
	outfile.write(script)

if __name__ == '__main__':
	main()

#RepeatMasker [input] -pa 10 -lib pier-1.3.fa (-q) -nolow -gff
