from subprocess import PIPE,Popen

# most of the time we'll be passing in the command as a list.
# list format can be useful for dynamic construction of a command

# all Popen commands (check_output, call, check_call, etc...) are WRAPPERS for Popen.
# so I'm just going to go through Popen(), not sure how the other functions handle concurrency.

print "////////Part 1, basic Popen//////////"
command = ['ls', '-l']
# Popen() returns a Popen object. The stdout flag is needed so that output is redirected to Popen.stdout.
process = Popen(command, stdout=PIPE)
# as we can see: this is of type 'file'
print type(process.stdout)
# we can do all file operations on it, such as:
for line in process.stdout:
	print line.strip()
print "////////Part 2//////////"
command = ['python', 'timer.py']
# once a Popen object is created, the process runs in the background.
# Most of the time, we need to wait for it to finish before continuing
process = Popen(command, stdout=PIPE)
print "This line prints immediately"
# querying process.stdout before the subprocess is finished stops this script until the subprocess finishes.
for line in process.stdout:
	print line.strip()
# the file object is not closed.

print "////////Part 3//////////"
command = ['python', 'timer.py']
# Popen.communicate() waits for the process to terminate, and also adds support for I/O with the process
process = Popen(command, stdout=PIPE)
output = process.communicate()
print "This line waits to print"
# if we use communicate(), the file object is closed, and the output is returned as a (stdout string, stderr string)
print type(output[0])
print output[0]


print "////////Part 4//////////"
# There are two (easy) ways to pipe:
print "VERY EASY, but UNSAFE. Applicable for 99% of situations:"
command = 'ls -l | grep .py'
process = Popen(command, stdout=PIPE, shell=True)
for line in process.stdout:
	print line.strip()
print "EASY, and SAFE. Applicable for 100% of situations:"
firstCommand = ['ls', '-l']
firstProcess = Popen(firstCommand, stdout=PIPE)
secondCommand = ['grep', '.py']
secondProcess = Popen(secondCommand, stdin=firstProcess.stdout, stdout=PIPE)
firstProcess.stdout.close() # important. see http://stackoverflow.com/questions/8369506/why-does-sigpipe-exist
output = secondProcess.communicate()[0]
print output

# Extra: auto login with ssh, passwords, etc. not possible with Popen alone (shell expects TTY).
# check out http://www.noah.org/python/pexpect/ for that kind of functionality
