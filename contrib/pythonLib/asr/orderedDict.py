"""
@author:     Arthur Kantor
@contact:     akantorREMOVE_THIS@uiuc.edu
@copyright:    Arthur Kantor 2009
@license:     GPL version 3
@date:         3/16/2009
@version:     0.1

"""

from util.markableList import MarkableList
import logging

class OrderedDict(list):
    '''a list of (word, defn) tuples.  defn is usually a list of str's or Units, but could be anything.  
    This is useful for holding a lexicon with optional definitions.  The sequential order of the words
     in the lexicon is tracked and preserved.
     
     If present defn is a MarkableList, otherwise None '''
    @classmethod
    def readDict(cls,fname):
    	'''Read the ordered dictionary from a file of format:
    	word1 phone1 phone2 ... phoneLast
    	word2 phone1 phone2 ... phoneLast
    	
    	phone1 phone2 ... phoneLast is optional, and if missing, the definition is None.
    	'''  
        phoneDict=cls()
        dictF=file(fname)
        i=0
        for l in dictF:
            ents=l.strip().split(None,1)
            word=ents[0]
            if len(ents)==2:
                phoneDict.append((word, MarkableList(ents[1].split())))
            else:
                phoneDict.append((word,None))           
            i += 1
        dictF.close()
        logging.info('Read %d word definitions from dictionary %s',len(phoneDict),fname)
        return phoneDict
    
    def joinOnWord(self,otherDict):
        '''@return a list l such that self[i][0]==otherDict[l[i]][0]'''
        wToI=otherDict.wordToIndex()
        return [ wToI[d[0]] for d in self]
    
    def wordToIndex(self):
        d= dict([ (v[0],i) for (i,v) in enumerate(self)])
        if len(d) != len(self):
            raise Exception("words are not unique: len(d) = %d and len(self) = %d"%(len(d), len(self)))
        return d 
       
    def update(self,moreDefsList):
    	'''If a def in moreDefsList is already in self, replaces the definition, otherwise, appends the new def'''
    	#d=Dict([(w,i) for (i,(w,d)) in enumerate(self)])
    	oldDict=self.wordToIndex()
    	for (w,d) in moreDefsList:
    		if w in oldDict:
    			self[oldDict[w]]=(w,d)
    		else:
    			oldDict[w]=len(self)
    			self.append((w,d))
