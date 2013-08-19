#!/usr/bin/python
# Author: Brian Lin
import unittest
from GCclassif import *
from classes import HmmHit

class Test_GCclassif(unittest.TestCase):

	def setUp(self):
		self.gagHmmHit = HmmHit('GAG', 100, 200, 20)
		self.apHmmHit = HmmHit('AP', 150, 250, 20)
		self.rtHmmHit = HmmHit('RT', 200, 250, 20)
		self.rhHmmHit = HmmHit('RH', 300,400, 20)
		self.intHmmHit = HmmHit('INT', 400,500,20)
		self.intHmmHit2 = HmmHit('INT', 0, 100, 40)
		self.rhHmmHit2 = HmmHit('RH', 50, 100, 40)
	
	def test_Gypsy_sequential(self):
		hitList = constructHitList(self.gagHmmHit, self.apHmmHit, self.rtHmmHit, self.rhHmmHit, self.intHmmHit)
		classif, score = classify(hitList)
		self.assertEqual('Gypsy', classif)
		self.assertEqual(100, score)
	
	def test_Copia_sequential(self):
		hitList = constructHitList(self.intHmmHit2, self.rhHmmHit2)
		classif, score = classify(hitList)
		self.assertEqual('Copia', classif)
		self.assertEqual(80, score)

	def test_Copia_nondeterministic(self):
		hitList = constructHitList(self.gagHmmHit, self.apHmmHit, self.intHmmHit2, self.rhHmmHit, self.rhHmmHit2)
		classif, score = classify(hitList)
		self.assertEqual('Copia', classif)
		self.assertEqual(80, score)
	
	def test_nondeterministic_behavior(self):
		hitList = constructHitList(self.intHmmHit2, self.gagHmmHit, self.apHmmHit, self.rtHmmHit, self.rhHmmHit, self.intHmmHit, self.intHmmHit2, self.rhHmmHit2)
		for hit in hitList: print hit
		classifs = []
		currentScore = 0
		classify_GAG(hitList, classifs, currentScore)
		print classifs
		self.assertEqual(len(classifs), 6) 
	
	'''
	def test_Gypsy_nondeterministic_scoring_diffClassifs(self):

	def test_Copia_nondeterministic_scoring_diffClassifs(self):

	def test_Gypsy_sequential_scoring_diffClassifs(self):

	def test_Copia_sequential_scoring_diffClassifs(self):

	def test_Gypsy_nondeterministic_scoring_sameClassif(self):

	def test_Copia_nondeterministic_scoring_sameClassif(self):
	'''

def constructHitList(*hits):
	hitList = []
	for hit in hits: 
		hitList.append(hit)
	return sorted(hitList)

if __name__ == '__main__':
	unittest.main()
