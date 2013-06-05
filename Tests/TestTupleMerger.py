import unittest
import TupleMerger

class TestTupleMerger(unittest.TestCase):

	def setUp(self):
		self.testTupleList1 = [(0, 20), (40, 60), (70, 75)]
		self.testTupleList2 = [(14, 24), (20, 29), (40, 49), (51, 60)]

	def testMergeTupleFromNewRange(self):
		newTuple = (30,35)
		result1 = TupleMerger.merge(self.testTupleList1, newTuple)
		result2 = TupleMerger.merge(self.testTupleList2, newTuple)
		self.assertEqual([(0, 20), (30, 35), (40, 60), (70, 75)], result1)
		self.assertEqual([(14, 29), (30, 35), (40, 49), (51, 60)], result2)

	def testMergeTupleWithinRange(self):
		newTuple = (42, 44)
		result1 = TupleMerger.merge(self.testTupleList1, newTuple)
		result2 = TupleMerger.merge(self.testTupleList2, newTuple)
		self.assertEqual([(0, 20), (40, 60), (70, 75)], result1)
		self.assertEqual([(14, 29), (40, 49), (51, 60)], result2)

	def testMergeTupleFromOverlap(self):
		newTuple = (25, 45)
		result1 = TupleMerger.merge(self.testTupleList1, newTuple)
		result2 = TupleMerger.merge(self.testTupleList2, newTuple)
		self.assertEqual([(0, 20), (25, 60), (70, 75)], result1)
		self.assertEqual([(14, 49), (51, 60)], result2)

	def testTotalLengthReturnsCorrectLength(self):
		result1 = TupleMerger.totalLength(self.testTupleList1)
		result2 = TupleMerger.totalLength(self.testTupleList2)
		self.assertEqual(45, result1)
		self.assertEqual(37, result2)

if __name__ == '__main__':
	unittest.main()
