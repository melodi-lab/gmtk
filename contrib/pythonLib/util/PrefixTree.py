"""A prefix tree (or a trie) module.
@author: 	Arthur Kantor
@contact: 	akantorREMOVE_THIS@uiuc.edu
@copyright:	Arthur Kantor 2008
@license: 	GPL version 3
@date: 		1/27/2009
@version: 	0.1


examples and unit test:
>>> trie = PrefixTree()
>>> trie.insert("abcd", "wordABCD")
>>> trie.insert("a", "wordA")
>>> trie.insert([1, 3, 5], "number135")

This returns the longest existing prefix, and the unprocessed remained of the key,
which in this case is 'b'
>>> trie.getLongestPrefix('ab')
(< ['a'] : wordA > , 'b')

>>> trie.getLongestPrefix('ABC')
(None, 'ABC')

>>> trie.getUpToLongestPrefix('ab')
([< ['a'] : wordA > ], 'b')

>>> trie.getUpToLongestPrefix('abcd')
([< ['a', 'b', 'c', 'd'] : wordABCD > , < ['a'] : wordA > ], '')


The following succeeds
>>> trie.get("abcd")
'wordABCD'

The following errors, unlike getLongestPrefix()
>>> trie.get("abc")
Traceback (most recent call last):
  File "<ipython console>", line 1, in <module>
  File "/home/arthur/thesis/fisher/scripts/pythonLib/util/PrefixTree.py", line 53, in get
    raise KeyError(restKey)
KeyError: 'bc'
"""

class PrefixTreeNode:
	""" A Prefix tree node.  
		The key can be anything that is an element in a sequence. 
		Each node can contain a value in addition to to child PrefixTreeNodes"""

	def __init__(self,parent=None,keyInParent=None):
		""" @param parent the pointer to the parent, if this is a non-root node"""
		self.parent=parent
		self.keyInParent=keyInParent
		self.k=dict()
		self.value=None

	def __repr__(self):
		return "< %s : %s > "%(self.getKey(),str(self.value))

	def getKey(self):
		"""The key required to get to this node.  The root node will return an empty list, so depth of the node is 1+len(returned key)"""
		if self.parent==None:
			return []
		else:
			l = self.parent.getKey()
			l.append(self.keyInParent)
			return l


	def insert(self,key,val):
		""" @param key a sequence object, to be used as a path into the prefix tree
			@param val the value to store for the key, overwriting the earlier value if any.
			@return the overwritten value if it was stored, or None otherwise."""

		if key:
			curk=key[0]
			if not self.k.has_key(curk):
				self.k[curk]=PrefixTreeNode(self,curk)
			return self.k[curk].insert(key[1:len(key)],val)
		else:
			ret = self.value
			self.value=val
			return ret

	def get(self, key):
		""" @param key a sequence object, to be used as a path into the prefix tree
			@return the value associated with the key
			@throws KeyError if the key is not found, with the value of the remaining key sequence"""
		(p,restKey)=self.getLongestPrefix(key)
		if len(restKey)>0:
			raise KeyError(restKey)
		else:
			return p.value

	def getLongestPrefix(self, key):
		""" Returns the PrefixTreeNode for the longest existing prefix, and the leftover key sequence
			@param key a sequence object, to be used as a path into the prefix tree
			@return the tuple (prefixTreeNode, sequence), where the prefixTreeNode is the longest existing prefix
			if no prefix exists, PrefixTreeNode can be None"""
		#print "getLongestPrefix(%s)"%(key)

		#by default we return us if we have a value, or None if we dont
		if not self.value == None:
			ret = (self,key)
		else:
			ret = (None, key)
		
		#we check deeper iff we want to (there is more key left) AND we can (we have the first element of the key in our hash)
		if key and self.k.has_key(key[0]):
			branch = self.k[key[0]].getLongestPrefix(key[1:len(key)])
			if branch[0] != None: #a longer prefix was found
				ret = branch 

		return ret
	
	def getUpToLongestPrefix(self, key):
		""" Same as getLongestPrefix, but also returns all existing shorter prefixes, in addition to
		    the longest prefix.
		    @return: the tuple ([prefixTreeNodeList], sequence), where the prefixTreeNodeList is the list
		    of all existing prefixes up to the longest existing prefix
			if no prefix exists, prefixTreeNodeList will be empty"""
		
		(n,rest) = self.getLongestPrefix(key)
		nl=[]
		while n != None:
			if n.value != None:
				nl.append(n) 
			n=n.parent
		return (nl,rest)

class PrefixTree(PrefixTreeNode):
	"""A PrefixTree is PrefixTreeNode with no parent."""
	def __init__(self):
		PrefixTreeNode.__init__(self,None,None)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

