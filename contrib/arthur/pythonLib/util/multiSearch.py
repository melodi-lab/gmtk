"""A fast implementation of searching a list for multiple sublists simultaneuously.
@author: 	Arthur Kantor
@contact: 	akantorREMOVE_THIS@uiuc.edu
@copyright:	Arthur Kantor 2009
@license: 	GPL version 3
@date: 		3/06/2009
@version: 	0.1


examples and unit test:
>>> 
"""

from util.PrefixTree import *

class MultiSearch:
	"""A container of sub-lists which will be searched for via the search method."""
	def __init__(self):
		self._t=PrefixTree()

	def insert(self,key,val):
		"""inserts a new sublist, and an object identifying it.  Inserting the same key multiple times
		   is allowed, as long as val is different.
		   @param key: a tuple sub-list
		   @param val: the object associated with the key
		"""
		try:
			self._t.get(key).append(val)
		except(KeyError):
			self._t.insert(key,[val])

	def search(self,li):
		"""searches li for sublists contained in this object.
			@param li: a sequence to be searched
			@return: a dict of val -> [int, int, ...] where val is the object
			associated with a found sublist,, and [int, int,...] is a list of one or
			more start indexes into li, where the sublist associated with val was found"""
		d={}
		for i in range(len(li)-1):
			(hitNodes, rest) = self._t.getUpToLongestPrefix(li[i:])
			for hits in [node.value for node in hitNodes]:
				for hit in hits:
					if hit in d:
						#print relevantMistakes[curM]
						d[hit].append(i)
					else:
						d[hit]=[i,]
		return d

if __name__ == "__main__":
    import doctest
    doctest.testmod()

