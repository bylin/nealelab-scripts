#!/usr/bin/python
# Author: Brian Lin
# SGE wrapper. Generates a helper file and runs it.

import subprocess, argparse, os, textwrap, datetime, shutil
global timestamp
timestamp = 'run-'+datetime.datetime.today().strftime("%a-%b-%Y-%H%M%S")

def main():
	args = parseArgs()
	setUpEnv(args)
	generateQsubScript(args)
	cmd = 'qsub '+getScriptName(args)
	subprocess.call(cmd, shell=True)

def parseArgs():
	mydesc = '''
	SGE wrapper. Contains options for qsub's built-in job arrays.
	
	Example usage: 
		SGEWrapper.py -m gcc repeatmasker -cmd 'RepeatMasker input_file.fa -pa 1 -lib library_file.fa -q -nolow -gff'
	
	or, if you want to run a command on multiple files in an input directory using SGE job arrays, 5 slots:
		SGEWrapper.py -j -tc 5 -d [input_directory] -m repeatmasker -f library_file.fa -cmd 'RepeatMasker $INPUT_FILE -lib library_file.fa -q -nolow -gff' '''
	argParser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(mydesc))
	argParser.add_argument('-j', '--use_job_array', help='Use a SGE job array', action='store_true')
	argParser.add_argument('-tc', '--max_running_tasks', help='Number of concurrent jobs to run at a time')
	argParser.add_argument('-d', '--input_directory', help='Directory containing all input files to be scheduled with the command in the job array')
	argParser.add_argument('-m', '--modules', nargs='+', help='Modules to load')
	argParser.add_argument('-f', '--extra_files', nargs='+', help='Extra files that need to be used by the program')
	argParser.add_argument('-o', '--stdout',  help='Standard output file', default='SGEWrapper.stdout')
	argParser.add_argument('-e', '--stderr',  help='Standard error file', default='SGEWrapper.stderr')
	argParser.add_argument('-N', '--name',  help='Job name', default='SGEWrapper')
	argParser.add_argument('-M', '--email', help='Email account user wishes to send job start and end notifications to')
	argParser.add_argument('-cmd', '--command',  required=True, help='Command to qsub')
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

def getName(args):
	if args.name is not None: return args.name
	else: return 'SGEWrapper'

def getScriptName(args):
	return getName(args)+'.sh'

def generateQsubArrayScript(args):
	outfile = file(getScriptName(args), 'w')
	(inputFile, nJobs) = scanInputs(args)
	script = addHeadings()
	script += addOptionals(args)
	script += '#$ -t 1-' + nJobs + '\n'
	script += addModules(args)
	#script += 'for f in '+args.input_directory+'*\n'
	#script += 'do\n\t'+args.command+'\ndone'
	script += addDate()
	script += 'INPUT_FILE=`sed -n "${SGE_TASK_ID}p" ' + inputFile + '`\n'
	script += args.command + '\n'
	script += addDate()
	outfile.write(script)

def generateQsubSingleJobScript(args):
	outfile = file(getScriptName(args), 'w')
	script = addHeadings()
	script += addOptionals(args)
	script += addModules(args)
	script += addDate()
	script += args.command
	script += addDate()
	outfile.write(script)

def scanInputs(args):
	inputFile = getName(args)+'InputFiles.txt'
	cmd = 'ls ../' + args.input_directory + '/* > ' + inputFile
	proc = subprocess.call(cmd, shell=True)
	lineCounter = 0
	inputFileHandle = open(inputFile)
	nJobs = sum(lineCounter + 1 for line in inputFileHandle)
	return(inputFile, str(nJobs))

def addHeadings():
	script = '#!/bin/bash\n'
	script += '#$ -S /bin/bash\n'
	script += '#$ -cwd\n'
	script += '#$ -q bigmem1.q\n'
	return script

def addOptionals(args):
	script = ""
	script += addEmail(args)
	script += addName(args)
	script += addStdout(args)
	script += addStderr(args)
	script += addMaxRunningTasks(args)
	return script	

def addEmail(args):
	script = ""
	if args.email is not None: 
		script += '#$ -m abe\n'
		script += '#$ -M '+args.email+'\n'
	return script


def addName(args):
	script = '#$ -N '+getName(args)+'\n'
	return script
		
def addStdout(args):
	script = ""
	if args.stdout is not None: 
		script += '#$ -o '+args.stdout+'\n'
	return script

def addStderr(args):
	script = ""
	if args.stderr is not None: 
		script += '#$ -e '+args.stderr+'\n'
	return script

def addMaxRunningTasks(args):
	script = ""
	if args.max_running_tasks is not None:
		script += '#$ -tc ' + args.max_running_tasks + '\n'
	return script

def addModules(args):
	script = ""
	if args.modules is not  None:
		script += 'module load '+ ' '.join(args.modules)+'\n' 
	return script

def addDate():
	return '\ndate\n'

if __name__ == '__main__':
	main()

