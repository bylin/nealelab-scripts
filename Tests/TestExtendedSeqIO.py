import unittest
import ExtendedSeqIO

class TestExtendedSeqIO(unittest.TestCase):

	def setUp(self):
		self.testFile = ExtendedSeqIO.ExtendedFastaFile('cprd.fa')

	def test_getSeqUsingDescription_returns_correct_seq(self):
		resultSeqRecord = self.testFile.getSeqUsingDescription('PtIFG7\tGypsy\tPinus taeda')
		self.assertEqual('PtIFG7\tGypsy\tPinus taeda', resultSeqRecord.description)

	def test_getSeqsUsingRegex_returns_correct_single_seq(self):
		resultSeqRecords = self.testFile.getSeqsUsingRegex('Gymny')
		firstRecord = resultSeqRecords[0]
		self.assertEqual('RLG_Gymny-1', firstRecord.name)

	def test_getSeqsUsingRegex_returns_correct_seqs(self):
		resultSeqRecords = self.testFile.getSeqsUsingRegex('IFG')
		for record in resultSeqRecords:
			print record.description
			self.assertRegexpMatches(record.description, 'IFG')

if __name__ == '__main__':
	unittest.main()
