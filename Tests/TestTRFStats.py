import unittest,trfStats,classes

class TestTRFStats(unittest.TestCase):
#Need to implement
	def setUp(self):
		pass
'''
	def testFilterTRFDict(self):
		trfDict = {'seq1': [classes.trfHit(['seq1','12','224','55','2.1','55','91','0','300','33','19','15','30','1.93','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATG','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATGCTCTAAAATCCTTGATTCATTTTTTAATGAGCTTGACACAAGGTAAGGAACTATGCTC']),classes.trfHit(['seq1','225','337','55','2.1','55','91','0','181','33','19','15','30','1.93','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATG','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATGCTCTAAAATCCTTGATTCATTTTTTAATGAGCTTGACACAAGGTAAGGAACTATGCTC'])]}
		trfDict = classes.trfHit(['seq2','225','337','55','2.1','55','91','0','181','33','19','15','30','1.93','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATG','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATGCTCTAAAATCCTTGATTCATTTTTTAATGAGCTTGACACAAGGTAAGGAACTATGCTC'])]}
	
		newTandem = classes.trfHit(['seq1','12','224','55','2.1','55','91','0','300','33','19','15','30','1.93','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATG','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATGCTCTAAAATCCTTGATTCATTTTTTAATGAGCTTGACACAAGGTAAGGAACTATGCTC'])
'''		

	def testOverlap(self):
		t0 = classes.trfHit(['seq1','12','224','55','2.1','55','91','0','300','33','19','15','30','1.93','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATG','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATGCTCTAAAATCCTTGATTCATTTTTTAATGAGCTTGACACAAGGTAAGGAACTATGCTC'])
		t1 = classes.trfHit(['seq1','225','337','55','2.1','55','91','0','181','33','19','15','30','1.93','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATG','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATGCTCTAAAATCCTTGATTCATTTTTTAATGAGCTTGACACAAGGTAAGGAACTATGCTC'])
		t2 = classes.trfHit(['seq1','230','330','55','2.1','55','91','0','180','33','19','15','30','1.93','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATG','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATGCTCTAAAATCCTTGATTCATTTTTTAATGAGCTTGACACAAGGTAAGGAACTATGCTC'])
		t3 = classes.trfHit(['seq2','225','337','55','2.1','55','91','0','181','33','19','15','30','1.93','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATG','CTCTAAAAACCTTGATTCATTGTTTAACGAGCTTGACACAAGACAAGGAACTATGCTCTAAAATCCTTGATTCATTTTTTAATGAGCTTGACACAAGGTAAGGAACTATGCTC'])

		self.assertFalse(trfStats.overlap(t0,t1))
		self.assertTrue(trfStats.overlap(t1,t2))
		self.assertFalse(trfStats.overlap(t0,t3))

if __name__ == '__main__':
	unittest.main()
