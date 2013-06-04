def merge(listOfLocationTuples):
		sortedTuples = sorted([sorted(t) for t in listOfLocationTuples])
		currentTuple = sortedTuples[0]
		for nextStart, nextEnd in sortedTuples[1:]:
				if nextStart <= currentTuple[1]:
						currentTuple[1] = max(currentTuple[1], nextEnd)
				else:
						yield tuple(currentTuple)
						currentTuple = [nextStart, nextEnd]
		yield tuple(currentTuple)