''' units (or contexts), tied-units (or replacements) and their containers  
@author: 	Arthur Kantor
@contact: 	akantorREMOVE_THIS@uiuc.edu
@copyright:	Arthur Kantor 2008
@license: 	GPL version 3
@date: 		3/04/2009
@version: 	0.1


examples and unit test:

'''

from util.Counter import *
from util.PrefixTree import *
from util.multiSearch import *
import logging
import math
import cPickle
from asr.Stream import *

			
class Unit(tuple):
	""" A contiguous string of phonemeIds(integers), plus context, each stored as a tuple."""

	(SUBWORD_TYPE, WHOLEWORD_TYPE, MULTIWORD_TYPE)=range(3)

	def __new__(cls,other,hypPhones=None):
		return super(Unit,cls).__new__(cls, other)

	def __init__(self,other,hypPhones=None):
		'''can be initialized either from another Unit, (a copy is made) or
		from a tuple and hypPhones.
		@param other: a Unit or a (center,left,right,wordIDlist) tuple
		@param hypPhones: a tuple of most frequent hypothesized phoneme IDs instead of the ref phoneme IDs
		''' 
				
		#self.count is the number of mistakes in the initial pass through utterances which 
		#may or may not be replaced 
		#it's not the same as sum(self.tot) because self.tot keeps track of units actually replaced

		#The following gets updated only by countCorrect() for estimating the replacement quality 
		self.tot=[0,0]
		self.completelyIncorrect=0
		#use these to get mean and variance of hypothesized phoneme starts per mistake for correct and incorrect
		self._totHypStarts=[0.0,0.0] 
		self._totHypStartsSq=[0.0,0.0]
		self._countHypPhones=Counter()
		
		#the following count actuall mistakes from utterance
		self.count=0
		self._hypPhones=Counter()

		if isinstance(other,Unit):
			self.addFrom(other)
			if hypPhones:
				raise ValueError('hypPhones is not allowed when initializing with unit')
		else:
			self.count=1
			self._hypPhones[hypPhones]=1
			
		#sets the unit type: 0 for subWord, 1 for wholeWord, 2 for multiword. Uses 
		#context EOW To distinguish between subWord and wholeWord 
		if Phoneme.EOW in self.center:
			self.type = self.MULTIWORD_TYPE
		elif len(self.left)>0 and len(self.right)>0 and self.left[-1] == self.right[0] == Phoneme.EOW:  
			self.type = self.WHOLEWORD_TYPE
		else:
			self.type = self.SUBWORD_TYPE
		

	#we also have self.multiWordName in the case self.type == self.MULTIWORD_TYPE

	def addFrom(self,otherUnit):
		"""Update this unit from otherUnit the (in)correct counts and hypPhone statistics are merged
			@param otherUnit: Where to augment statistics from"""
		for i in range(2):
			self.tot[i] += otherUnit.tot[i]
			self._totHypStarts[i] += otherUnit._totHypStarts[i]
			self._totHypStartsSq[i] += otherUnit._totHypStartsSq[i]
		self.count += otherUnit.count
		self.completelyIncorrect += otherUnit.completelyIncorrect
		self._hypPhones.update(otherUnit._hypPhones)
		self._countHypPhones.update(otherUnit._countHypPhones)
	"""The chances of making the recognition units better/worse after the replacement when the occurance of the replaced
	phonemes was incorrect/correct"""
	betterChance=.50
	worseChance= .50
	
	@property
	def center(self):
		return self[0]

	@property
	def left(self):
		return self[1]

	@property
	def right(self):
		return self[2]
  
	@property
	def wordlist(self):
		return self[3]

	def __repr__(self):
		st=[]
		for s in self[:3]:
			ss=[Phoneme.names[p] for p in s] + []
			st.append(tuple(ss))
		st.append(self[3])
		return str(tuple(st))
	
	@property
	def improvement(self):
		"""returns net improvement on the unit occurrences in the corpus
		@precondition: countCorrect() has been run"""
		return (self.betterChance*self.tot[1]-self.worseChance*self.tot[0])#*len(self.center)

	@property
	def stats(self):
		"""returns a stats string about this unit"""
		try:
			precision=self.tot[1]/float(sum(self.tot))
		except(ZeroDivisionError):
			precision=float('NaN')
		return "%s \t%d:%.2f(%.2f) \t%d:%.2f(%.2f), %d %f %f" \
			%(self,\
			  self.tot[0], self.meanHypStarts()[0],self.varHypStarts()[0], \
			  self.tot[1], self.meanHypStarts()[1],self.varHypStarts()[1], \
			  self.count, precision, self.improvement)#,self.countHypPhoneHist())

	def countHypPhoneHist(self):
		s = sorted([ (v, tuple(Phoneme.names[n] for n in k)) for (k,v) in self._countHypPhones.items()])
		r=[ (k,v) for (v,k) in s]
		r.reverse()
		return r
	
	@classmethod
	def mistakesFromUtterance(cls,utt, contextType):
		"""@return  a list of mistakes from in a given utterance.  A mistake is a Unit where the center 
		is a consecutive sequence of incorrect phones.   
		
		@todo: NOT DONE The context for each mistake is the correct phonemes for the words enclosing the mistake,
		plus one more word to the left and to the right.  So if the utterance is 
		<s> word1 worXX XXrd3 word4 </s> The left context will be the phones for word1 wor, and
		right context is rd3 word4.  The context is shortened if it's too close to the utterance boundary
		"""

		mistakes=[]
		#find longest incorrect subseq
		end=-1
		while True:
			start=utt.findCorrect(False,end+1)
			if start==None:
				break
			end=utt.findCorrect(True,start+1)
			if end == None:
				end=len(utt)
			#clean up the EOW's on the edges - make them always part of context
			#regardless of whether they are part of the mistake or not
			(corrStart,corrEnd)=(start,end)
			#print (corrStart,corrEnd)
			if utt[start].phonemeId==Phoneme.EOW:
				corrStart += 1
			if utt[end-1].phonemeId==Phoneme.EOW:
				corrEnd -= 1
			#by now, both the start and end of sequence are guaranteed not to be EOW
			if corrStart<corrEnd:
				left=[]
				right=[]
				if contextType == 'none':
					pass
				elif contextType == 'eow':
					if corrStart>0 and utt[corrStart-1].phonemeId == Phoneme.EOW:
						left.append(Phoneme.EOW)
					if corrEnd<len(utt) and utt[corrEnd].phonemeId == Phoneme.EOW:
						right.append(Phoneme.EOW)
				elif contextType == 'containingWords':
					for i in range(corrStart-1,-1,-1):
						if utt[i].phonemeId == Phoneme.EOW:
							break
					left =[p.phonemeId for p in utt[i:corrStart]]
					for i in range(corrEnd,len(utt)):
						if utt[i].phonemeId == Phoneme.EOW:
							break
					right= [p.phonemeId for p in utt[corrEnd:i+1]]
					
				else:
					assert False
					
				#get the list of word ids spanned by the center.
				#we may have different words generate identical phone 
				#sequences, and we may want to track that.
				prevWordId=-1
				wordIds=[]
				for i in range(corrStart,corrEnd):
					if prevWordId == -1:
						wordIds.append(utt[i].wordId)
						prevWordId = utt[i].wordId
					if utt[i].phonemeId == Phoneme.EOW:
						prevWordId = -1
				
				#find the hypothesized phones for each mistake 
				hypPhoneIds=[]
				for i in range(corrStart,corrEnd):
					hypPhoneIds.append(utt[i].hypMode)

				mistakenPhonemeIds=tuple(p.phonemeId for p in utt[corrStart:corrEnd])
				#print ((mistakenPhonemeIds,tuple(left),tuple(right),tuple(wordIds)),hypPhoneIds)
				mistakes.append(Unit((mistakenPhonemeIds,tuple(left),tuple(right),tuple(wordIds)),tuple(hypPhoneIds)))
		return mistakes
	
	def projectContext(self, contextType):
		"""returns a copy of the unit with context reduced to a smaller (more general) type.
		The copy does not have countCorrect() run on it,and it has no
		(correct,incorrect,completelyIncorrect) statistics.
		@todo: only contextType none is implemented for now."""
		
		if contextType == 'removeWordId':
			wordlist=()
			left=self.left
			right=self.right
		elif contextType == 'projectToMistake':
			left=()
			right=()
			wordlist=()
		elif contextType == 'none':
			wordlist=self.wordlist
			left=self.left
			right=self.right
		else:
			raise NotImplementedError()
		
		return Unit((self.center,left,right,wordlist))
	
		 
	def countCorrect(self, utt, starts):
		"""Update statistics for the number of correct/incorrect occurrences of self in context on the given utterance
		(tot[0],tot[1],completelyIncorrect) gather the counts of phone substrings
		which are either all correct, or incorrect and the last is the completely incorrect subset of incorrect
		side effect: also sets self._totHypStarts and self._totHypStartsSq for mean and variance calculations"""
		

		#FIXME:  if context is more complicated, do this differently
		(unitStart,unitEnd)=(len(self.left),len(self.left)+len(self.center))
		isMulti = self.type == self.MULTIWORD_TYPE

		for match in starts:
			(start,end)=(match+unitStart,match+unitEnd)
			firstMark=utt.findMarked(True,start,end)
			if firstMark != None: 
				#the match was previously marked, so don't count it again
				#because it will be replaced by some other unit.
				#Skip forward to the end of marks and continue search
				nextStart=utt.findMarked(False,firstMark)
				if nextStart == None:
					nextStart=len(utt)
			else: #the match is not marked. 
				c=utt.isCorrect(start,end)
				hypStarts=sum(p.hypStarts for p in utt[start:end])
				hypStartsSq = hypStarts**2
				if c == True:
					#already is correct - there is a chance replacing unit will harm score
					ci=0
				else:
					ci=1
					if c == False:
						self.completelyIncorrect += 1
					hypPhones=tuple(p.hypMode for p in utt[start:end])
					if not self._countHypPhones[hypPhones]:
						self._countHypPhones[hypPhones] =0 
					self._countHypPhones[hypPhones]+=1
				self.tot[ci] += 1
				self._totHypStarts[ci] += hypStarts
				self._totHypStartsSq[ci] += hypStartsSq
				
				#normally we mark just replaced part, but for multi-words we mark the multi-word and 
				#the entire context.  This removes the units which occur only as a part longer units
				if isMulti:
					utt.mark(True,match,match+len(self.left)+len(self.center)+len(self.right))
				else:
					utt.mark(True,start,end)
				nextStart=end

	def meanHypStarts(self):
		meanHypStarts=[0.0,0.0]
		for i in range(2):
			if self.tot[i]>0:
				meanHypStarts[i]=self._totHypStarts[i]/self.tot[i]
			else:
				meanHypStarts[i]=-1.0
		return meanHypStarts
	
	def varHypStarts(self):
		varHypStarts=[0.0,0.0]
		for i in range(2):
			if self.tot[i]>0:
				m=self.meanHypStarts()
				varHypStarts[i]=self._totHypStartsSq[i]/self.tot[i]-m[i]**2
			else:
				varHypStarts[i]=-1.0
		return varHypStarts

		
	
class TiedUnit(Unit):
	""" An actual physical replacement unit which may be used to replace phonemes in multiple contexts """

	def _subUnitCountSumComponentPhonemes(self):
		c = sum(Phoneme.subPhones[p] for p in self.center)
		return c
	
	def _subUnitCountHeur(self):
		"""heuristically guess an appropriate number of sub-unit states for this unit"""
		c = sum(Phoneme.subPhones[p] for p in self.center if p != Phoneme.EOW)
		#it seems that most of the mistakes are likely reductions, so reduce the number of substates
		if c>=6:
			c=int(c*.75)
		return c
		
	def _subUnitCountMeanObserved(self):
		"""Estimate the number of units by looking at the number of hypothesized incorrect units for this unit"""
		#I am guessing may be 2.8 sub-units per unit - dipthongs are shorter some feel like they can be shorter.
		#print "***"+str(self)+str(self.meanHypStarts())
		if self.meanHypStarts()[1] == 0:
			logging.warning("subUnitCount is estimated heuristics for %s, %s %s",self.meanHypStarts(),self.tot,self)
			#a special unit for something weird happening 
			return self._subUnitCountHeur()
		u = self.meanHypStarts()[1]
		f = 3
		return round(u*f)

	subUnitCountMethods={'asBefore' : _subUnitCountSumComponentPhonemes,'meanObserved' : _subUnitCountMeanObserved}
	'''How to guess the number of substates in the units'''
	@classmethod 
	def setSubUnitCountMethod(cls, meth):
		cls.subUnitCount=cls.subUnitCountMethods[meth]
		cls.methName=meth

	@classmethod 
	def getSubUnitCountMethod(cls):
		return cls.methName

	def __init__(self,t):

		"""All the full-context units which will be replaced by this one."""
		super(TiedUnit, self).__init__(t)
		self.sourceUnits=UnitDict()


	@property
	def name(self):
		"""The name constructed from phonemes replaced"""
		parts = ['_'.join(Phoneme.names[p] for p in s) for s in self[0:3]]
		return '.'+'.'.join(parts)+'.'
	

	def addSource(self,sourceUnit):
		"""Add the unit which will be replaced by this  unit. the in/correct counts are merged
		@param tiedUnit: The source unit to be replaced by this unit.  This unit can replace multiple source units"""
		sourceUnit.replacement=self
		self.sourceUnits[sourceUnit]=sourceUnit
		self.addFrom(sourceUnit)

	def removeSource(self,sourceUnit):
		"""inverse of addSource
			@raise KeyError: if the source does not use this as tied unit """
		if sourceUnit not in self:
			raise KeyError('%s does not use %s as tied unit'%(sourceUnit, self))
		
		del sourceUnit.replacement
		del self.sourceUnits[sourceUnit]
		for i in range(2):
			self.tot[i] -= sourceUnit.tot[i]
			self._totHypStarts[i] -= sourceUnit._totHypStarts[i]
			self._totHypStartsSq[i] -= sourceUnit._totHypStartsSq[i]
		
		
	@property	
	def stats(self):
		r =super(TiedUnit,self).stats + '\t %s %d, %f'%(self.name,self.subUnitCount(), self.risk)
		for s in self.sourceUnits.sortedByImprovement():
			r += '\n\t'+s.stats
		return r
	
	@property
	def risk(self):
		"""Takes into account the improvement and the cost"""
		return self.improvement/self.subUnitCount()
	
		
class UnitContainer:
	""" A bag (multi-set) of Units."""

	#def growWithSubstrings(self):
	#	""" @return a UnitBag containing all the substrings of the Units.  The identical
	#		units/substrings of units are merged and their counts are summed.  The removed
	#		prefix and suffix are moved into the left and right context"""

	def sortedByCoverage(self):
		"""@return a list of Units m sorted in descending order by len(m.center)*count	"""
		def sortKey(item):
			return -len(item[0].center)*item[1],item[0]
			
		return sorted(self.iteritems(),key=sortKey)
	
	def sortedByLength(self):
		"""@return a list of Units m sorted in descending order first by longest units, and if units are the same length
		then by the most specific context."""
		def sortKey(item):
			return -len(item.center), -len(item.left), -len(item.right), -item.count, item.center, item.left, item.right, item.wordlist
			
		return sorted(self,key=sortKey)
	
	def sortedByImprovement(self,worseChance=.05, betterChance=.5):
		"""@precondition: countCorrect() has been run"""
		
		def sortKey(item):
			return -item.improvement

		return sorted(self,key=sortKey)

	def countCorrect(self, utterances):
		"""Gets statistics on how many correct/incorrect phoneme strings would be affected in utterances
		by making the substitution with this set of units.  The substitutions are performed in order,
		first substituting the longest units, and if units are the same length, then substituting the
		most specific context.
		@return: 
		@precondition: phonemes in utterances are marked False
		@postcondition: same as precondition"""
		
		for m in self:							
			m.tot=[0,0]
			m._countHypPhones.clear()
		
		#build a trie of units and their contexts		
		t=MultiSearch()
		for m in self:
			k = m.left+m.center+m.right
			t.insert(k,m)
		logging.info("mistake trie constructed")
		#this should be the same as running all the utterances
		#through each mistake in turn, where the mistakes are
		#in sortedByLength() order

		for u in utterances:
			phones =[p.phonemeId for p in u]
			relevantMistakes=UnitDict(t.search(phones))
			for m in relevantMistakes.sortedByLength():
				m.countCorrect(u,relevantMistakes[m])
			u.mark(False)
			
	def removeSpurious(self):
		"""@precondition: countCorrect must have been run"""
		logging.info("removing superfluous units from among %d units",len(self))
		misKeys=self.keys()
		for m in misKeys:
			if sum(m.tot)==0:
				#all of these mistakes were a substring of a longer mistake already added - so drop this mistake
				del self[m]
		logging.info("got %d units after removing 0-occurrence units",len(self))
	
		misKeys=self.keys()
		ignore=[(i,) for i in range(Phoneme.SIL,Phoneme.EOW+1)]
		for m in misKeys:
			if m.center in ignore:
				del self[m]
			elif Phoneme.SIL in m.center:
				del self[m]
			elif Phoneme.EOW in m.center and (Phoneme.SIL in m.left or Phoneme.SIL in m.right):
				del self[m]
		logging.info("got %d units after removing non-speech single-phoneme units and units containing SIL",len(self))


	@property
	def improvement(self):
		"""returns net improvement on the corpus from replacing all the unit occurrences in this set
		@precondition: countCorrect() has been run"""
		return sum(m.improvement for m in self)
		
		
	def getTiedUnits(self, contextType):
		"""returns a UnitDict of TiedUnits where the sources are from this bag.
		side effect: sets the replacement for each contained unit"""
		tied=TiedUnitDict()
		for u in self: 
			t=u.projectContext(contextType)
			if t not in tied:
				newT=TiedUnit(t)
				tied[newT]=newT
			tied[t].addSource(u)
		return tied				

	def stats(self):
		totalImprovement=0.0
		totalReplacements=0
		unitTypes=[0,0,0]
		for m in self:
			totalImprovement  += m.improvement
			totalReplacements += sum(m.tot) 
			unitTypes[m.type] += 1

		return "Total improvement: %f\tUnit tokens: %d \tUnit Types %d: (%d subWords, %d words, %d multiwords)"%\
			tuple([totalImprovement,totalReplacements,sum(unitTypes)]+unitTypes)

#class UnitBag(Counter,UnitContainer):
#	pass

class UnitSet(dict,UnitContainer):
	'''A set of units, with the property that if a unit is added, and it exists,
	the hypPhonemes and other statistics are updated in the existing unit.
	The key is the base tuple of the object, and the value is the object itself
	'''
	def __init__(self, iterable=None):
		if iterable:
			self.update(iterable)
	
	def update(self,other):
		for o in other:
			self.updateItem(o)
			
	def updateItem(self,o):
		#print o
		if o in self:
			self[o].addFrom(o)
		else:
			newO=Unit(o)
			self[newO]=newO
	
	#'''you probably don't want to use these.'''
	#def __setitem__(self,k,v):
	#def fromkeys(self,S,v=None):
	
class UnitDict(dict,UnitContainer):
	pass

class TiedUnitDict(UnitDict):
	'''The Keys are the base-class tuple of the TiedUnit object, 
	   and the values are the objects themselves.'''

	def totalStates(self):
		return sum(r.subUnitCount() for r in self)

	def stats(self):
		(totalStates, totalStatesHeur,diff)=(0,0,0)
		for r in self.values():
			totalStates += r.subUnitCount()
			#totalStatesHeur += r.subUnitCountHeur()
			#diff += math.fabs(r.subUnitCount() - r.subUnitCountHeur())
		return super(TiedUnitDict,self).stats()+"\tstates: %d"%totalStates



class FinalUnitsPickle(object):
	''' read and write finalUnits pickle which is used by a number of programs
		contexts is of class UnitSet
		repl is of class TiedUnitDict
	'''
	def __init__(self):
		self.contexts=UnitSet()
		self.repl=TiedUnitDict()
	
	@classmethod
	def readFromFile(cls,fname):
		o=FinalUnitsPickle()
		inuf=file(fname)
		subUnitsHeur=cPickle.load(inuf)
		TiedUnit.setSubUnitCountMethod(subUnitsHeur)
		o.contexts=cPickle.load(inuf)
		o.repl=cPickle.load(inuf)
		inuf.close()
		return o
	
	def writeToFile(self,fname):
		f=open(fname,'w')
		self.writeToIO(f)
		f.close()
		
	def writeToIO(self,outf):
		subUnitsHeur=TiedUnit.getSubUnitCountMethod()
		cPickle.dump(subUnitsHeur,outf,-1)
		cPickle.dump(self.contexts,outf,-1)
		cPickle.dump(self.repl,outf,-1)

	def info(self):
		info = 'subUnitsHeur: %s\nreplacements:\n'%(TiedUnit.getSubUnitCountMethod())
		info += self.repl.stats()
		info += 'replacement details:\n'
		for r in self.repl.sortedByImprovement():
			info += '%s\n'%r.stats	
		return info


if __name__ == "__main__":
    import doctest
    doctest.testmod()

