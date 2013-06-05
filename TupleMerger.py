# TupleMerger.py
# Authors: Hans Vasquez-Gross, John Liechty, Brian Lin

def merge(tupleList, newTuple):
	return list(mergeGenerator(tupleList, newTuple))

def mergeGenerator(tupleList, newTuple):
	tupleList.append(newTuple)
	sortedTuples = sorted([sorted(t) for t in tupleList])
	currentTuple = sortedTuples[0]
	for nextStart, nextEnd in sortedTuples[1:]:
		if nextStart <= currentTuple[1]:
			currentTuple[1] = max(currentTuple[1], nextEnd)
		else:
			yield tuple(currentTuple)
			currentTuple = [nextStart, nextEnd]
	yield tuple(currentTuple)

def totalLength(tupleList):
	return sum((end-start) for (start,end) in tupleList)
