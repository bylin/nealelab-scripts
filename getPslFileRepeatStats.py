#!/usr/bin/python
from classes import PierRepeatElement

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

class RepeatStats(object):
  BpsByRepeatClassifs = {}

  def addRepeatCopy(repeat, blockSize): # repeat is of type PierRepeatElement
    if repeat.CLASS in BpsByRepeatClassifs:
      BpsByRepeatClassifs[repeat.CLASS] += blockSize
    if repeat.SUPER in BpsByRepeatClassifs:
      BpsByRepeatClassifs[repeat.SUPER] += blockSize
    if repeat.ORDER in BpsByRepeatClassifs:
      BpsByRepeatClassifs[repeat.ORDER] += blockSize

  def returnRepeatStats(self):
    outputString = ''
    for repeat in BpsByRepeatClassifs:
      outputString += repeat + ': ' + BpsByRepeatClassifs[repeat] + '\n'
    return outputString


if __name__ == '__main__':
  filename = checkAndReturnArguments()
  for line in filename:
    tabs = filename.split('\t')
    seqName = tabs[9]
    currentHit = PierRepeatElement(seqName)
    blockSizes = convertCSVsToIntList(tabs[18]) # tabs[18] is a CSV string
    for blockSize in blockSizes:
      RepeatStats.addRepeatCopy(currentHit)
  print RepeatStats.returnRepeatStats()
