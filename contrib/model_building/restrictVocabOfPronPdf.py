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



     
def makeParser():
 usage="""
 restrictVocabOfPronPdf [--help] [options]  <inGMTKParamFile> <outGMTKParamFile>
 Finds the pronunciation densePDF matrix, in inGMTKParamFile, with each 
 ith row corresponding to the ith row in inSortedVocabFile , trims it so
 that the only the rows remaining are the ones in outSortedVocabFile, 
 sorted as in outSortedVocabFile.
 if UNK is missing from the in-words, it is added to the in-words, and a 
 corresponding pronunciation line [1 0 0 ...] is added to the incoming pronMatrix.
 """

 parser = OptionParser()
 parser.usage=usage
 
 parser.add_option("-i", "--inGMTKParamFile", type="string", 
                  help="<REQUIRED>")
 parser.add_option("-o", "--outGMTKParamFile", type="string", 
                  help="<REQUIRED>")
 parser.add_option("-v", "--inOrderedVocabFile", type="string", 
                  help="<REQUIRED>")
 parser.add_option("-w", "--outOrderedVocabFile", type="string", 
                  help="<REQUIRED>  This file is read to determine the order of the words in the pronunciation matrix of outGMTKParamFile")
 parser.add_option("-m", "--minPronProbThreshold", type="float", default="0.0001",
                  help="If a pronunciation variant's probability is below this threshold, it gets floored to 0, with the probability being equaly distributed to the remaining pronunciation variants [default: %default]")
 return parser

def main(argv):
 opts, args = parser.parse_args()
 if not (opts.inGMTKParamFile and opts.outGMTKParamFile and opts.inOrderedVocabFile and opts.outOrderedVocabFile):
    parser.error("need all of the following defined: inGMTKParamFile outGMTKParamFile inOrderedVocabFile outOrderedVocabFile")
 
 #m=DCPT(Workspace(),'pronunciationMatrixTemp',[2, 3],2)
 
 #for i in range(2):
     #for j in range(3):
        #m[:i,j] = [.1*i+.01*j]*2
 
 #m.writeToFile(WorkspaceIO(sys.stdout))
 #m.changeParentCards([3,2])
 #m.writeToFile(WorkspaceIO(sys.stdout))
 #m.changeParentCards([4,4])
 #m[:3,3]=[50, 500]
 #m.writeToFile(WorkspaceIO(sys.stdout))
 #sys.exit()
 
 inOrderedVocab=file(opts.inOrderedVocabFile)
 l=inOrderedVocab.readlines()
 #l=('a','b','c')
 wordToIndex=dict(zip(l,range(len(l))))
 #print wordToIndex
 inOrderedVocab.close()
 
 outOrderedVocab=file(opts.outOrderedVocabFile)
 l=outOrderedVocab.readlines()
 #l=('a','c')
 orderedOutWords=dict(zip(range(len(l)),l))
 outOrderedVocab.close()
 #print orderedOutWords
 print('Read %d in-ordered-words, %d out-ordered-words' % (len(wordToIndex), len(orderedOutWords)))


 wk = Workspace()
 print('Reading master file: %s' % opts.inGMTKParamFile)
 wk.readMasterFile(opts.inGMTKParamFile);
 pronMat=wk.objects[DCPT]['pronunciationMatrix']
 
 if "unk\n" not in wordToIndex:
     print 'adding unk (and its pronunciation) to the in-words...'
     unkInd=len(wordToIndex)
     wordToIndex['<unk>\n']=unkInd
     pronMat.changeParentCards([len(wordToIndex)+1])
     #pronMat[:unkInd]=[0]*pronMat.getSelfCards()[0]
     pronMat[0:unkInd]=1
     
 if "<s>\n" not in wordToIndex and "<S>\n" in wordToIndex: 
     print 'renaming <S> to <s> in the in-words set...'
     wordToIndex['<s>\n']=wordToIndex['<S>\n']
     del wordToIndex['<S>\n']
     print 'renaming </S> to </s> in the in-words set...'
     wordToIndex['</s>\n']=wordToIndex['</S>\n']
     del wordToIndex['</S>\n']    
 newPronMat=DCPT(wk,'pronunciationMatrixTemp',len(orderedOutWords),pronMat.getSelfCards()[0])
 
 for i,w in orderedOutWords.items():
     if w in wordToIndex:
         newPronMat[:i]=thresholdProb(pronMat[:wordToIndex[w]], opts.minPronProbThreshold)
     else :
         raise LookupError("the desired out-word %s is not in the list of in-words" % w )
 newPronMat.name='pronunciationMatrix'
 wk.objects[DCPT]['pronunciationMatrix'] = newPronMat
 del wk.objects[DCPT]['pronunciationMatrixTemp'] 
 
 #print pMat.__repr__()
 #print pMat[:1]
 #pMat[:1]=[10, 10, 10, 10, 10]
 #print pMat[:1]
 #del pMat[:1]
 #w.objects[DCPT]['newpronunciationMatrix']=DCPT(w,
 print('Writing master file: %s' % opts.outGMTKParamFile)
 wk.writeMasterFile(opts.outGMTKParamFile);
 #print w.__str__();

def thresholdProb(pIn, thresh):
 p=list(pIn)
 leftovermass=0.0
 nonZeroProbs=0
 for i,v in enumerate(p):
     if v<thresh:
         leftovermass += v
         p[i]=0
     else:
         nonZeroProbs += 1
         
 extra=leftovermass/nonZeroProbs
 
 for i,v in enumerate(p):
     if v>=thresh:
         p[i] += extra
 return p



#the parser is used for generating documentation, so create it always, and augment __doc__ with usage info  
#This messes up epydoc a little, but allows us to keep a single version of documentation for all purposes
parser = makeParser()
__doc__ = __doc__.replace("%InsertOptionParserUsage%\n", parser.format_help())

if __name__ == "__main__":
    main(sys.argv)
