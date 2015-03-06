#!/usr/bin/python
"""

%InsertOptionParserUsage%


@author: 	Arthur Kantor
@contact: 	akantorREMOVE_THIS@uiuc.edu
@copyright:	Arthur Kantor 2008
@license: 	GPL version 3
@date: 		11/18/2008
@version: 	0.9

"""

import sys, re, os;
from optparse import OptionParser
from gmtkParam import *
from util import *
from util.UniqueList import *


def usage():
 return """
genSubState.py [--help] [options] 
Define a random variable (RV) that remains constant over a sequence of frames in terms of sub-states of that variable.

This can be used for:
	1) Defining an utterance in terms of it's component words
	2) Defining a word in terms of it's component phones, (or sub word units), with one pronuciation per word
	3) Defining a phoneme in terms of its component sub-phone states

If you are defining a RV STATE, as a sequence of up to K sub states assigned to RV subSTATE,
with the RV subSTATEindex iterating from 0 to k, the max number of sub states for the given state , with k<K.
The following is determined and reported:
	1) K, the cardinality of subSTATEindex RV
	2) the cardinality of subSTATE RV
	3) the decision tree (DT) for subSTATE given subSTATEindex and STATE
	4) (optionally) the DT for binary RV STATEtransition, given subSTATEtransition, subSTATEindex and STATE,
		which is 1 only if subSTATEtransition is 1 and subSTATEindex is 'pointing' to the last sub-state for STATE
		and false otherwise

The dictonary file defines STATE in terms of a sequence of subSTATEs, and contains the one definition per line 
with  typcal format
	LINE := STATE ' '  DEFN
	STATE := word 
	DEFN := word | word ' ' DEFN 

if --autoGenerateSubStateNames is specified, the DEFN in the typical format is changed to 
	DEFN := integer
If the entry in the dictonary is (for example)
	AA	3
then it is equivalent to 
	AA AA_0 AA_1 AA_2


"""

def makeParser():
	parser = OptionParser()
	parser.usage=usage()
	
	# parser.add_option("-i", "--inTrainedParamFile", type="string", 
	#				  help="<REQUIRED>")
	# parser.add_option("-v", "--inOrderedVocabFile", type="string", 
	#				  help="<REQUIRED>")
	parser.add_option("-d", "--dictionary", type="string", metavar="FILE",	help="<REQUIRED> The dictonary file ")

	parser.add_option("--noHead", action="store_false", default=True, dest='hasHead',
					help="""If false, generates subSTATE given only subSTATEindex (subSTATE does not not depend 
							on STATE), and generates a decision tree per every definition in the dict file, instead 
							of one large decision tree. The STATE becomes the name of the DT. Use for defining 
							utterances as strings of words. Both --subStateDT and --stateTransitionDT are
							affected. (default: %default)""")
							
	parser.add_option("-a", "--autoGenerateSubStateNames", action="store_true", default=False,
					help="""If specified, auto-generates names for subSTATE. Use for defining phonemes as 
							strings of sub-phoneme-states.  (default: %default)""")
					
	parser.add_option("-v", "--varMap", type="string", metavar="FILE",
					help="""if specified, creates a varMap file, one subSTATE name per line, corresponding to the index used in the DT,
							so the name on the first line corresponds to the value 0 assigned to subSTATE""")
	
	parser.add_option("-u", "--useVarMap", type="string", metavar="FILE",
					help="""if specified, checks the definitions in --dictionary against the values in 
							the varMap file FILE, (same format as --varMap).  If a definition uses a word not in FILE, exits with error""")
	
	parser.add_option("-s", "--subStateDT", type="string", metavar="FILE",
					help="""<REQUIRED> Where to store the DT for subSTATE given subSTATEindex and STATE.  In case of --noHead, the multiple DTs are stored here """)
					
	parser.add_option("-t", "--stateTransitionDT", type="string", metavar="FILE",
					help="""If specified, generates STATEtransition and stores it in FILE.  In case of --noHead, the multiple 
							DTs are stored here.""")
					
	parser.add_option("-e", "--enforceConsistentLastSubState", action="store_true", default=False,
					help="""If specified, verifies that the last subState in each definition is the same.
							If this is not the case, exits with error.
							This is useful e.g. if utterances should have a special END_OF_UTTERANCE word </s>.
							In this case you probably don't need to generate --stateTransitionDT.""")
	
	parser.add_option("-c", "--cardinalities", type="string", metavar="FILE",
					help="""If specified, stores cardinalities of STATE, subSTATEindex, subSTATE and the 
							index of the last subState if it is consistently used in all definitions to FILE""")
					
	#parser.add_option("-n", "--stateName", type="string", metavar="STR", default='STATE',
					#help="""State name, used in naming decision trees""")
					
	#parser.add_option("-m", "--subStateName", type="string", metavar="STR",
					#help="""derived from stateName, as sub<stateName>, by default, its subSTATE""")
					
	return parser

def parseOpts():						
	opts, args = parser.parse_args()
	if not (opts.dictionary and opts.subStateDT):
		parser.error("need all of the following defined: dictionary subStateDT")
	
	#if not opts.subStateName:
		#opts.subStateName='sub'+opts.stateName
	
	return opts
	
def main(argv):
	cmd=' '.join(argv)
	opts = parseOpts()
	
	comment= '%'+" generated with command: %s\n\n"%cmd
	
	if opts.hasHead:
		name = os.path.basename(opts.subStateDT).rsplit('.',1)[0]
		ssDt = HeadedDTBuilder( name, 2, opts.subStateDT, comment,)
	else:
		ssDt = HeadlessDTBuilder(opts.subStateDT, comment)
	
	if opts.stateTransitionDT:
		if opts.hasHead:
			name = os.path.basename(opts.stateTransitionDT).rsplit('.',1)[0]
			stDt = HeadedDTBuilder( name, 3, opts.stateTransitionDT, comment,)
		else:
			stDt = HeadlessDTBuilder(opts.stateTransitionDT, comment)
	
	if opts.useVarMap:
		f=open(opts.useVarMap)
		varMap = ValidatingVarMap([l.strip() for l in f])
		f.close()
	else:
		varMap = PermissiveVarMap()
	
 	dp = opts.autoGenerateSubStateNames and AutoGenDefinitionParser() or TypicalDefinitionParser()
	
	df=open(opts.dictionary)
	defns = DictReader(df,dp,varMap)

	print("opened files...");

	lastSubState=[None,True] #The last subState and whether it has been consistently used
	for d in defns:
		#print d
		if not lastSubState[0] == d[1][-1] and lastSubState[1]:
			if lastSubState[0] == None:	#d is the first line of dict
				lastSubState[0] = d[1][-1]
			else: #inconsitency
				lastSubState[1] = False
				if opts.enforceConsistentLastSubState:
					raise ValueError("The definition %s does not end in %s as all previous definitions do." %(d,lastSubState[0]))

		ssDt.appendDT(SubStateDT(d,varMap))
		if opts.stateTransitionDT:
			stDt.appendDT(StateTransitionDT(d,varMap))

	df.close()

	print("generated DTs...");
	
	ssDt.finalize()
	
	if opts.stateTransitionDT:
		stDt.finalize()
		
	
	if opts.varMap:
		f=open(opts.varMap,'w')
		varMap.writeToFile(f)
		f.close()

	if opts.cardinalities:
		f=open(opts.cardinalities,'w')
		f.write("STATE_CARD %d\nsubSTATE_CARD %d\nsubSTATEindex_CARD %d\n"%(defns.stateCount, len(defns.subStateSet), defns.maxDefLen))
		if lastSubState[1]:
			f.write("LAST_subSTATE_TOKEN %s\nLAST_subSTATE %d\n" % (lastSubState[0],defns.subStateSet.index(lastSubState[0])))
		f.close()

	print("wrote DTs...");

	#print (len(defns.subStateSet), defns.maxDefLen)
	sys.exit(0)

class SubStateDT(DT):
	"""Parent 0 is subStateIndex"""
	def __init__(self, defn, subStateSet):
		default=TreeBranch(-1,TreeLeaf(subStateSet.index(defn[1][-1])),defn[1][-1])
		tree = TreeBranch(0,default, defn[0])
		for i,d in enumerate(defn[1][0:-1]):
			tree[i]=TreeBranch(-1,TreeLeaf(subStateSet.index(d)),d)
			
		DT.__init__(self, dict(),defn[0],1,tree)
		
class StateTransitionDT(DT):
	"""Parent 0 is subStateIndex, Parent 1 is subStateTransition"""
	def __init__(self, defn, subStateSet):
		lastIndex=len(defn[1])-1
		if lastIndex == 0:
			tree = TreeBranch(-1,TreeLeaf('{p1}'))		
		else:
			default=TreeBranch(-1,TreeLeaf(0))
			tree = TreeBranch(0,default, defn[0])
			tree[lastIndex]=TreeBranch(-1,TreeLeaf('{p1}'))
			
		DT.__init__(self, dict(), defn[0],  2, tree)
		
class DTBuilder:
	"""abstract builder of DT or DTs objects"""
	def __init__(self, fileName, comment=''):
		self.f=file(fileName,'w')
		self.ws=WorkspaceIO(self.f)
		self.ws.writelnWord(comment)
	
	def appendDT(self, dt):
		raise NotImplementedError()
	
	def finalize(self):
		self.f.close()
	
		
	
class HeadlessDTBuilder(DTBuilder):
	def __init__(self,  fileName, comment=''):
		DTBuilder.__init__(self,  fileName, comment)
		self.dtCount=0
		self.f.write('\n')
		self.dtCountOffset = self.f.tell()
		self.f.write(' '*10+'%%num of DTs')

	
	def appendDT(self, dt):
		self.f.write('\n%d\n'%self.dtCount)
		self.dtCount +=1
		dt.writeToFile(self.ws)
		
	def finalize(self):
		dtCountStr=str(self.dtCount)
		if len(dtCountStr)>10:
			raise Exception("Too many DTs.  grow the blank space buffer")
		
		self.f.seek(self.dtCountOffset)
		self.f.write(dtCountStr)
		DTBuilder.finalize(self)


class HeadedDTBuilder(HeadlessDTBuilder):
	"""A DTs with a single DT, which gets new children via appendDT.
	   This DT is in an inconsistant state until finalize or writeToFile is called"""
	def __init__(self, name, numParents, fileName, comment=''):
		HeadlessDTBuilder.__init__(self,  fileName, comment)
		
		"""The last parent is the one on which The head decision is made on the last parent (numParents-1)"""
		tree=TreeBranch(numParents-1)
		self.headDT=DT(dict(), name, numParents, tree)
	
		self.tree=self.headDT.getTree()
	
	def finalize(self):
		lastBranch=self.tree.numQuestions()-1
		self.tree.default=self.tree[lastBranch]
		del self.tree[lastBranch]
	
		HeadlessDTBuilder.appendDT(self, self.headDT)
		HeadlessDTBuilder.finalize(self)
	
	def appendDT(self, dt):
		self.tree.append(dt.getTree())
		
	#def writeToFile(self, stream):
		#self.finalize()
		#HeadlessDTBuilder.writeToFile(self, stream)

	#def saveToFile(self, fileName, comment=''):
		#"""requires writeToFile to be implemented"""
		#f=WorkspaceIO(file(fileName,'w'))
		#f.writelnWord(comment)
		#self.writeToFile(f)
		#f.close()

class DictReader:
	"""An iteratable returning a definition (state, [substate_0,substate_1,..]) whenever a line() is called.
		Also gathers statistics for cardinalities of subSTATEindex and subSTATE """
	def __init__(self,dict, defnParser, subStateSet):
		self.dict = dict
		self.defnParser = defnParser
		self.subStateSet=subStateSet
		self.maxDefLen=0
		self.stateCount=0
	
	def __iter__(self):
		for l in self.dict:
			words = l.split()
			state = words.pop(0)
			defn=self.defnParser.parseDef(state, words)
			self.maxDefLen=max(self.maxDefLen,len(defn))
			self.stateCount +=1
			if not self.subStateSet.check(defn):
				raise ValueError("definition '%s' has words not in listed in the varMap file"%defn)
			yield (state, defn)
		raise StopIteration



		
		
		
		
class VarMap(UniqueList):
	"""Base class for a VarMap."""
	
	def check(self, subStates):
		"""returns true if subStates is a sequence of allowable var names, otherwise raises an exception"""
		raise NotImplementedError()

	def writeToFile(self,file):
		for v in self:
			file.write("%s\n"%v)
	
class PermissiveVarMap(VarMap):
	def check(self, subStates):
		"""always returns true, and extends self with previously unseen subStates"""
		self.extend(subStates)
		return True
	
class ValidatingVarMap(VarMap):
	def __init__(self, allowedNames):
		UniqueList.__init__(self, allowedNames)
		
	def check(self, subStates):
		"""returns true if subStates are in the allowedNames, false otherwise"""
		return self.issuperset(subStates)
		
		
		
class DefinitionParser:
	def parseDef(self, state, defn):
		raise NotImplementedError()
	
class TypicalDefinitionParser(DefinitionParser):
	def parseDef(self, state, defn):
		return defn

class AutoGenDefinitionParser(DefinitionParser):
	def parseDef(self, state, defn):
		return ["%s_%d" % (state,i) for i in range(int(defn[0]))]

#the parser is used for generating documentation, so create it always, and augment __doc__ with usage info  
#This messes up epydoc a little, but allows us to keep a single version of documentation for all purposes
parser = makeParser()
__doc__ = __doc__.replace("%InsertOptionParserUsage%\n", parser.format_help())


if __name__ == "__main__":
	
	main(sys.argv)

