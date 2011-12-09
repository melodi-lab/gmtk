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
from string import Template

def usage():
 return """
multiPronToSinglePron.py [--help] [options] 
Reads the pronunciation densePDF matrix, inOrderedVocabFile, with each 
with the inOrderedVocabFile sorted by word and pronunciation, where word 
and within-word-pronunciation order matches the pronunciation matrix's 
row and column orders respectively.
Picks the pronunciation with highest probability and generates a phonetic 
dictionary, with the one most probable pronunciation per line.

See generate_phoneStateDT.pl for a script that generates a suitable 
inWord_pronunciation_wordStateCounter_2_state_DT_file
"""

	
	 
def makeParser():
 parser = OptionParser()
 parser.usage=usage()
 
 parser.add_option("-i", "--inTrainedParamFile", type="string", 
				  help="<REQUIRED>")
 parser.add_option("-v", "--inOrderedVocabFile", type="string", 
				  help="<REQUIRED>")
 parser.add_option("-d", "--inDictFile", type="string", 
				  help="<REQUIRED>  ")
 parser.add_option("-o", "--outSinglePronDictFile", type="string", 
				  help="<REQUIRED>  ")
 return parser

def main(argv):
 opts, args = parser.parse_args()
 if not (opts.inTrainedParamFile and opts.inOrderedVocabFile and opts.outSinglePronDictFile):
	parser.error("need all of the following defined: inTrainedParamFile inOrderedVocabFile outSinglePronDictFile")
 
 #Precondition: the dictonary is word and pronunciation order matches the pronunciation matrixes row and column orders respectively
 dictFile=open(opts.inDictFile)
 
 lastWord=''
 dic=dict()
 for l in dictFile:
	word,defn = l.split(None,1)
	if word==lastWord:
		dic[word].append(defn)
		#print orderedDict[-1]
	else:
		dic[word]=[defn]
		lastWord=word
 print 'read %d words from %s'% (len(dic),opts.inDictFile)
 
 orderedVocabFile=open(opts.inOrderedVocabFile)
 orderedWords= [l.strip() for l in  orderedVocabFile]
 print 'read %d ordered words from %s'% (len(orderedWords),opts.inOrderedVocabFile)
 
 wk = Workspace()
 wk.readTrainableParamsFile(opts.inTrainedParamFile);
 pronMat=wk.objects[DCPT]['pronunciationMatrix']
 print('read trainable params file: %s' % opts.inTrainedParamFile)
 
 outDict=open(opts.outSinglePronDictFile,'w')
 for i in range(len(dic)):
	maxProb,maxInd=max( (prob,i) for i,prob in enumerate(pronMat[:i]))
	#print (i,orderedWords[i],"%s %s"%(orderedWords[i],dic[orderedWords[i]][maxInd]),pronMat[:i])
	outDict.write("%s %s"%(orderedWords[i],dic[orderedWords[i]][maxInd]))
 outDict.close()
 
 sys.exit()


#the parser is used for generating documentation, so create it always, and augment __doc__ with usage info  
#This messes up epydoc a little, but allows us to keep a single version of documentation for all purposes
parser = makeParser()
__doc__ = __doc__.replace("%InsertOptionParserUsage%\n", parser.format_help())

if __name__ == "__main__":
	main(sys.argv)

