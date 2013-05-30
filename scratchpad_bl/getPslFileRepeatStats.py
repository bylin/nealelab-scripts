#!/usr/bin/python

def checkAndReturnArguments():
  if (len(sys.argv) != 1):
    print "Usage: getPslFileRepeatStats.py [input psl file]"
  return sys.argv[1]

def convertCSVsToIntList(csvString):
  stringList = csvString.split(',')
  intList = []
  for value in stringList:
    intList.append(int(value))
  return intList

class repeatElement:
  namedRepeats = 
    def __init__(self, string):
      if string == 
      self.CLASS = 
    def parseString(string):
      return 

class repeatStats:
  def __init__(self):
    self.classes = ['I', 'II', 'Unclassified']
    self.orders = ['']

  def addHitWithNbps(hit, nbps):


if __name__ == '__main__':
  filename = checkAndReturnArguments()

  for line in filename:
    tabs = filename.split('\t')
    annotation = tabs[9]
    blockSizesString = tabs[18]
    currentHit = repeatElement(annotation)
    blockSizes = convertCSVsToIntList(blockSizesString)
    for blockSize in blockSizes:

