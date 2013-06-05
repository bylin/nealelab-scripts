def merge(tuples, newTuple):
	return list(mergeGenerator(tuples, newTuple))

def mergeGenerator(tuples, newTuple):
	tuples.append(newTuple)
	sortedTuples = sorted([sorted(t) for t in tuples])
	currentTuple = sortedTuples[0]
	for nextStart, nextEnd in sortedTuples[1:]:
		if nextStart <= currentTuple[1]:
			currentTuple[1] = max(currentTuple[1], nextEnd)
		else:
			yield tuple(currentTuple)
			currentTuple = [nextStart, nextEnd]
	yield tuple(currentTuple)

def totalLength(tuples):
	return sum((end-start) for (start,end) in tuples)
