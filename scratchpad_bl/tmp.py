from subprocess import PIPE,Popen


firstCommand = ['ls', '-l']
firstProcess = Popen(firstCommand, stdout=PIPE)
secondCommand = ['grep', '".py"']
print "hmmm"
secondProcess = Popen(secondCommand, stdin=firstProcess.stdout, stdout=PIPE)
firstProcess.stdout.close() # important. see http://stackoverflow.com/questions/8369506/why-does-sigpipe-exist
for line in secondProcess.stdout:
        print line.strip()
