"""A list that contains only unique elements.  If an element that already exists is put 
into the list, the action succeeds but only the already existing element is kept

@author: 	Arthur Kantor
@contact: 	akantorREMOVE_THIS@uiuc.edu
@copyright:	Arthur Kantor 2008
@license: 	GPL version 3
@date: 		11/18/2008
@version: 	0.9
"""
class UniqueList(list):
	"""A complete implementation of the list interface.  All functions succeed iff the same 
	   function would succeed in a typical list, but no duplicates are ever inserted.
	"""
	
	def __init__(self, seq=[]):
		self.s=dict()
		self.extend(seq)
	
	def __delitem__(self, i):
		del self.s[self[i]]
		list.__delitem__(self, i)

	def __delslice__(self,i,j):
		for k in range(i,j):
			del self.s[self[k]]
		list.__delslice__(self,i,j)
	
	def __iadd__(self, seq):
		self.extend(seq)

	def __imul__(self, n):
		raise NotImplementedError()
	
	def __setitem__(self,i,y):
		if y not in self.s:
			del self.s[self[i]]
			self.s[y]=i
			list.__setitem__(self,i,y)
 
	def __setslice__(self,i, j, y):
		for k in range(i,j):
			self[k]=y[k]
	
	def __contains__(self, el):
		return el in self.s
	
	def append(self, el):
		if el not in self:
			self.s[el]=len(self)
			list.append(self, el)
		
	def extend(self, seq):
		for el in seq:
			self.append(el)
			
	def index(self, value, start=0, stop=None):
		if stop == None:
			stop=len(self)
		if value in self and self.s[value]<stop and self.s[value] >= start:
			return self.s[value]
		else:
			raise ValueError("UniqueList.index(x): x not in list")
		

	def insert(self, index, object):
		if not object in self.s:
			list.insert(index, object)
			for i,v in enumerate(self):
				self.s[v]=i
 
	def pop(self, index=None):
		v=list.pop(index)
		del self.s[v]
		return v

	def remove(self,value):
		del self[self.s[value]]

	def reverse(self):
		list.reverse()
		for i,v in enumerate(self):
			self.s[v]=i

	def sort(self, cmp=None, key=None, reverse=False):
		list.sort(cmp,key,reverse)
		for i,v in enumerate(self):
			self.s[v]=i
	
	#set functions
	def issuperset(self, other):
		for i in other:
			if not i in self.s:
				return False
		return True
	
	def issubset(self, other):
		for i in self.s:
			if not i in other:
				return False
		return True
	
	def union(self, other):
		return set(self.s).union(other)
	
	def intersection(self, other):
		return set(self.s).symmetric_difference(other)
	
	def difference(self, other):
		return set(self.s).difference(other)
	
	def symmetric_difference(self, other):
		return set(self.s).symmetric_difference(other)
	