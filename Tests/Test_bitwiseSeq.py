from bitwiseSeq import *
from memory_profiler import profile
import unittest, copy

class Test_bitwiseSeq(unittest.TestCase):
	
	def setUp(self):
		self.seq1 = BitwiseSeq(30)
		self.seq2 = BitwiseSeq(3000)
		self.seq3 = BitwiseSeq(300000)

	def test__init__(self):
		seq = BitwiseSeq(10)
		self.assertEqual(seq.bitmax, 1023)

	def test_convertIntToBitint(self):
		result = convertIntToBitint(5)
		self.assertEqual(result, 16)
		result = convertIntToBitint(3)
		self.assertEqual(result, 4)

	def test_convertBigIntToBitint(self):
		result = convertIntToBitint(5015280)
		answer = pow(2, 5015279)
		self.assertEqual(result, answer)
		result = convertIntToBitint(21000)
		answer = pow(2, 20999)
		self.assertEqual(result, answer)

	def test_convertBlockToBitint(self):
		result = convertBlockToBitint((2, 5))
		self.assertEqual(result, 14)
		result = convertBlockToBitint((12, 20))
		self.assertEqual(result, 522240)
		
	def test_convertEmptyBlockToBitint(self):
		result = convertBlockToBitint((4, 4))
		self.assertEqual(result, 0)

	def test_countFlippedBits(self):
		result = countFlippedBits(29)
		self.assertEqual(result, 4)
		result = countFlippedBits(11)
		self.assertEqual(result, 3)
		result = countFlippedBits(0)
		self.assertEqual(result, 0)

	def test_verifyValidBlockRange(self):
		blocks = [(24, 30), (150, 250), (150, 30000)]
		results = []
		results.append(self.seq1.verifyBlockRange(blocks[0]))
		results.append(self.seq2.verifyBlockRange(blocks[1]))
		results.append(self.seq3.verifyBlockRange(blocks[2]))
		self.assertEqual(results, [True, True, True])
		block = (14, 14)
		result = self.seq1.verifyBlockRange(block)
		self.assertEqual(result, True)

	def test_verifyInvalidBlockRange(self):
		blocks = [(-3, 12), (24, 21), (30, 1000), (0, 3)]
		result = self.seq1.verifyBlockRange(blocks[0])
		self.assertEqual(result, False)
		result = self.seq1.verifyBlockRange(blocks[1])
		self.assertEqual(result, False)
		result = self.seq1.verifyBlockRange(blocks[2])
		self.assertEqual(result, False)
		result = self.seq1.verifyBlockRange(blocks[3])
		self.assertEqual(result, False)

	def test_addValidBlock(self):
		self.seq1 += (10, 15)
		result = lambda: self.seq1.bitint
		answer = convertBlockToBitint((10, 15))
		self.assertEqual(result(), answer)
		self.seq1 += (12, 20)
		answer = convertBlockToBitint((10, 20))
		self.assertEqual(result(), answer)
		self.seq1 += (8, 25)
		answer = convertBlockToBitint((8, 25))
		self.assertEqual(result(), answer)
		self.seq1 += (12, 15)
		self.assertEqual(result(), answer)

#	def test_addEmptyBlock(self):

if __name__ == '__main__':
	unittest.main()


