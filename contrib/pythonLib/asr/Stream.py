"""Phoneme and Utterance IO streams
@author: 	Arthur Kantor
@contact: 	akantorREMOVE_THIS@uiuc.edu
@copyright:	Arthur Kantor 2008
@license: 	GPL version 3
@date: 		3/04/2009
@version: 	0.1


examples and unit test:

"""

from operator import *

class Phoneme:
	"""A reference phoneme, and some statistics about it."""

	#these are from fihser/scr/phone.state file
	names=( 'AA', 'AE', 'AH', 'AO', 'AW', 'AY', 'B',  'CH', 'D',  'DH',\
			'EH', 'ER', 'EY', 'F',  'G',  'HH', 'IH', 'IY', 'JH', 'K',\
			'L',  'M',  'N',  'NG', 'OW', 'OY', 'P',  'R',  'S',  'SH',\
			'T',  'TH', 'UH', 'UW', 'V',  'W',  'Y',  'Z',  'ZH', 'SIL',\
			'NOISE', 'LAUGH', 'BREATH', 'COUGH', 'LIPSMACK', 'SIGH', 'SNEEZE', 'EOW')

	subPhones=(3, 3, 3, 3, 2, 2, 3, 3, 3, 3,\
			   3, 3, 2, 3, 3, 3, 3, 3, 3, 3,\
			   3, 3, 3, 3, 2, 2, 3, 3, 3, 3,\
			   3, 3, 3, 3, 3, 3, 3, 3, 3, 3,\
			   3, 3, 3, 3, 3, 3, 3, 1)

	@classmethod
	def phonemeIdListToName(cls, li):
		return '_'.join([cls.names[p] for p in li])
	
	#The special EOW phoneme
	EOW = 47

	#The special Silence phoneme
	SIL = 39

	def __init__(self,phonemeId,wordId,uttId,frames,isCorrect,hypMode, hypStarts):
		"""id of the phoneme in the corpus"""
		self.phonemeId=phonemeId
		self.wordId=wordId
		self.uttId=uttId
		
		"""The number of frames it represents"""
		self.frames=frames
		
		"""At least one frame in the hypothesis matches the self.phonemeId"""
		self.isCorrect=isCorrect
		
		"""The mode of hypothesis phonemeIds"""
		self.hypMode=hypMode
		
		"""The number of hypothesis phoneme starts for the duration of the reference phoneme.
		   If ref and hypothesis phonemes change on the same frame, it is considered as a new 
		   hypothesis phoneme start.  If only the reference phoneme begins and hyp phoneme 
		   remains the same as in the previous frame, it is not a new hypothesis phoneme start."""
		self.hypStarts=hypStarts
        
		"""a convenience variable to mark the phoneme during certain traversals"""
		self.mark=False
        
	def __repr__(self):
		return "Phoneme(%d, %d, %d, %d, %d, %d, %d, %s)"% \
			(self.uttId,self.phonemeId,self.hypMode, self.hypStarts, self.frames,self.wordId,self.isCorrect, self.mark)
		
	@classmethod
	def _phonemesFromFrames(cls, f):
		"""	@param f the file with the corpus
			@return A generator function reading lines in the corpus and spits out phones.
			Assumes that words are delimited, by some special phoneme otherwise not occuring within a word.
			Also assumes that two identical phonemes can appear consecutively within a word."""

		prevPhone=-1
		hypPhonesHist=dict()
		lastHypPhone=-1
		for l in f:
			(uttId,frameNo,hypPhone,refPhone,word)= [int(x) for x in l.strip().split()]
			if(refPhone != prevPhone):
				if prevPhone>-1:
					(modeCount,mode) = max([(c,p) for (p,c) in hypPhonesHist.items()])
					yield cls(prevPhone,prevWord,prevUttId,frames,correct,mode, hypStarts)

				(prevPhone, prevWord, prevUttId) = (refPhone, word, uttId)
				frames=0
				hypStarts=0
				hypPhonesHist.clear()
				correct=False
			
			frames += 1
			if hypPhone != lastHypPhone:
				hypStarts += 1
				lastHypPhone= hypPhone
			hypPhonesHist[hypPhone] = hypPhone in hypPhonesHist and hypPhonesHist[hypPhone]+1 or 1
			correct = correct or refPhone == hypPhone
		
		if prevPhone>-1:
			(modeCount,mode) = max([(c,p) for (p,c) in hypPhonesHist.items()])
			yield cls(prevPhone,prevWord,prevUttId,frames,correct,mode, hypStarts)

	@classmethod
	def _phonemesFromFile(cls, f):
		""" @read phonemes from file, phoneme per line, in format 
			phoneme  word utteranceId frameCount isCorrect mostFrequentHypothesized phone"""
		for l in f:
			(phone,word,uttId,frames,correct,mode,hypStarts)= [int(x) for x in l.strip().split()]
			yield cls(phone,word,uttId,frames,correct,mode,hypStarts)
			
	_inf = None
	@classmethod
	def PhonemeStream(cls, opts):
		if cls._inf != None:
			cls._inf.seek(0)
		else:
			cls._inf = open(opts.inFile)
			 
		if opts.inputFormat == 'frames':
			ing = cls._phonemesFromFrames(cls._inf)
			ing = cls._silEOWSilToSilFilter(ing)
		elif opts.inputFormat == 'phonemes':
			ing = cls._phonemesFromFile(cls._inf)
		else:
			assert False
		return ing

	@classmethod
	def _silEOWSilToSilFilter(cls,phonemeIter):
		""" Merges the phoneme sequence SIL EOW SIL to SIL.  This I suspect is some (harmless?) bug in the GMTK 
		forced alignment script with the </s> word"""
		phoneStack=[]
		for p in phonemeIter:
			if len(phoneStack)==0:
				if p.phonemeId==cls.SIL:
					phoneStack.append(p)
				else:
					yield p
			elif len(phoneStack)==1: #SIL is pushed
				if p.phonemeId==cls.EOW and p.uttId==phoneStack[0].uttId and p.wordId==phoneStack[0].wordId:
					phoneStack.append(p)
				else:
					yield phoneStack.pop(0)
					yield p
			elif len(phoneStack)==2: #SIL and EOW are pushed
				if p.phonemeId==cls.SIL and p.uttId==phoneStack[0].uttId and p.wordId==phoneStack[0].wordId:
					#merge SIL EOW SIL into SIL by summing the frames, summing hypStarts, ORing the correction and 
					#taking the mode from the one with more frames.  (This is not entirely correct, but we don't 
					#really use it)
					phoneStack.append(p)
					frames=0
					hypStarts=0
					isCorrect=False
					for ph in phoneStack:
						frames += ph.frames
						hypStarts += ph.hypStarts 
						isCorrect = isCorrect or ph.isCorrect
					hypModeIdx=max(enumerate([ph.frames for ph in phoneStack]), key=itemgetter(1))[0]
					hypMode=phoneStack[hypModeIdx].hypMode
					newP=Phoneme(cls.SIL,p.wordId,p.uttId,frames,isCorrect,hypMode,hypStarts)
					phoneStack=[newP]
				else:
					yield phoneStack.pop(0) #empty the stack
					yield phoneStack.pop(0)
					if p.phonemeId==cls.SIL: #and push the new SIL if needed
						phoneStack.append(p)
					else:
						yield p
			else:
				assert len(phoneStack)>2

		while phoneStack:
		  yield phoneStack.pop(0)
		
	def writeToFile(self, f):
		l= "%d\t%d\t%d\t%d\t%d\t%d\t%d\n"%(self.phonemeId,self.wordId,self.uttId,self.frames,int(self.isCorrect),self.hypMode, self.hypStarts)
		f.write(l)
		
class Utterance(list):
	"""A list of references phonemes denoting an utterance"""

	@classmethod
	def utteranceStream(cls, phonemes):
		"""@return: a generator of utterances 
		   @param phonemes: a sequence of phonemes"""

		u = cls()
		u.append(phonemes.next())
		for p in phonemes:
			if p.uttId != u[-1].uttId:
				yield u
				u=Utterance([p])
			else:
				u.append(p)
		yield u

	def findCorrect(self, correctness=True, start=0,end=None):
		"""Find the first element with correctness value in [start,end)
		@return the index of the element or None otherwise"""

		if end == None:
			end=len(self)

		for i in range(start,end):
			if self[i].isCorrect == correctness:
				return i
		return None

	def findMarked(self, mark=True, start=0,end=None):
		"""Find the first element with correctness value in [start,end)
		@return the index of the element or None otherwise"""

		if end == None:
			end=len(self)

		for i in range(start,end):
			if self[i].mark == mark:
				return i
		return None


	def find(self, seq, start=0, end=None):
		"""Find the phoneme string seq, as a substring of self, optionally starting at start, and ending at end.
		@return The minimum index of the start of a matching substring if found, None otherwise"""
		
		if end==None:
			end=len(self)
		seq=tuple(seq)
		phones =tuple(p.phonemeId for p in self)
		for i in range(start,end-len(seq)+1):
			if phones[i:i+len(seq)] == seq:
				return i
		
		return None

	def isCorrect(self, start=0, end=None):
		"""@return: True if the whole range is Correct, False if the whole range is Incorrect, None otherwise"""
		correct=True
		incorrect=True
		
		if end == None:
			end=len(self)
		
		for p in self[start:end]:
			correct = correct and p.isCorrect
			incorrect = incorrect and not p.isCorrect
		if correct:
			return True
		elif incorrect:
			return False
		else:
			return None
	
	def mark(self, value, start=0,end=None ):
		"""Mark the range of phones with the specified value"""

		if end == None:
			end=len(self)

		for p in self[start:end]:
			p.mark=value
			
	@property
	def uttId(self):
		return self[0].uttId
	
	def __repr__(self):
		return "%d: %s" %(self.uttId,str(self.phonemeString))
	
	@property
	def phonemeString(self):
		return [Phoneme.names[p.phonemeId] for p in self]


if __name__ == "__main__":
    import doctest
    doctest.testmod()

