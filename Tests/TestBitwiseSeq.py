from bitwiseSeq import *
import unittest

class TestBitwiseSeq(unittest.TestCase):
	
	def setUp(self):
		seq = BitwiseSeq(1000)

	def testBasePairIsEmpty(self):
		result = seq.BasePairIsEmpty(2)
		self.assertEqual(result, True)

	def testAddBlockToEmptyRegion(self):
		seq.addBlockToEmptyRegion(10, 20)

if __name__ == '__main__':
	unittest.main()

