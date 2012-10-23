"""
@author:     Arthur Kantor
@contact:     akantorREMOVE_THIS@uiuc.edu
@copyright:    Arthur Kantor 2009
@license:     GPL version 3
@date:         3/16/2009
@version:     0.1

"""

class MarkableList(list):
    ''' A list where elements can be marked with an object'''
    
    def __init__(self,li):
        super(MarkableList,self).__init__(li)
        self._mrk=[None,]*len(li)

    def replace(self,replObj, length, starts):
        """replace a seq of items of given length at given starts, with replObj,
        ensuring that replaced seqs don't overlap.
        @param replObj: The replacement object
        @param length: non-negative int
        @param starts: a list of non-negative ints""" 
        for s in starts:
            if self._findMarked(True, s, s+length)==None:
                self._mark((replObj,), s, s+length)

    def getWithReplacements(self):
    	''' @return: a list of objects, where the object is the original element in the list, if the element is not marked
    	    otherwise, the mark object.  If a sequence of elements is marked with the same object, only a single mark object
    	    is returned, thus the returned list can be shorter than len(self)
    	''' 
        r=[]
        lastRep=None
        for i in range(len(self)):
            if lastRep==None or not (lastRep is  self._mrk[i]):
                if self._mrk[i]==None:
                    r.append(self[i])
                else:
                    r.append(self._mrk[i][0])
                lastRep=self._mrk[i]
        return r
    
    def getMark(self,i):
    	''' @return: the mark for ith element in the list'''
    	return self._mrk[i]           
	   
    def _findMarked(self, mark=True, start=0,end=None):
        """Find the first marked element with boolean value mark in [start,end)
        @return the index of the element or None otherwise"""

        if end == None:
            end=len(self)

        for i in range(start,end):
            if bool(self._mrk[i]) == mark:
                return i
        return None
        
    def _mark(self, value, start=0,end=None ):
        """Mark the range of phones with the specified value"""

        if end == None:
            end=len(self)

        for i in range(start,end):
            self._mrk[i]=value

