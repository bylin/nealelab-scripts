#!/usr/bin/python
import cPickle as pickle

def sendToPickleJar(source,outfile):
	pickle.dump(source,open('/home/jjzieve/Pita_Genome-0.9_Repeats/pickle_jar/'+outfile,'w'))
	return

def getFromPickleJar(infile):
	return pickle.load(open('/home/jjzieve/Pita_Genome-0.9_Repeats/pickle_jar/'+infile,'r'))

#def genPickler(outfile):
#	fh = open
#	pickler = pickle.Pickler(fh)
#	return pickler
#
#def genPickle(outfile,source):
#	pickle = genPickler(outFile)
#	for item in source:
#		yield pickle.dumps(item)
#
#def genUnpickle(infile):
#	while True:
#		try:
#			item = pickle.load(infile)
# 			yield item
#		 except EOFError:
#			 return

