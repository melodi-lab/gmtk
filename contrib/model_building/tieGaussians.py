#!/usr/bin/python
"""
%InsertOptionParserUsage%


@author: 	Arthur Kantor
@contact: 	akantorREMOVE_THIS@uiuc.edu
@copyright:	Arthur Kantor 2009
@license: 	GPL version 3
@date: 		4/23/2009
@version: 	0.1

"""

import sys, re, os;
import logging
import time
import argparse
import operator 
from gmtkParam import *
from util import *
from util.UniqueList import *
from asr.Unit import *
from asr.Stream import Phoneme
import copy
import textwrap
import numpy
from numpy import *




def makeParser():
	usage = """\
		tieGaussians.py [--help] [options] 
		Takes a trainable parameters file, with one gaussian per mixture, and ties the mixtures via 
		bottom up clustering. Gaussians i ~ N(mu_i,sigma_i) and j ~ N(mu_j,sigma_j) are tied if 
		d(i,j) is less than some threshold, where d(i,j) is eq 17.1 in the HTK book.
		
		Also reads/writes the named collections corresponding to the in/out gmtk param files.
		Only the non-phoneme units (of form e.g. gm:plp:.V.EOW_Y_UW.EOW.:0 1 1 are affected).
		In the above case, only mixtures having the same value between the first dots (.V.) 
		are considered
		
		The name format for input gaussian mixtures is:
			gm:featureKind:.mistake.leftContext.rightContext.:subUnitState e.g.
			gm:plp:.V.EOW_Y_UW.EOW.:0
		
		The output will be of form 
			gm:plp:.V.EOW_Y_UW.EOW.:0
		"""

	def readableFile(fname):
		if os.access(fname,os.R_OK):
			return fname
		else:
			print "cannot read file %s"%fname
			raise TypeError()

	parser = argparse.ArgumentParser(prog='tieGaussians.py', description=textwrap.dedent(usage),formatter_class=argparse.RawDescriptionHelpFormatter)
	
	parser.add_argument( "-i", "--inTrainedParamFile",type=readableFile, metavar="FILE", required=True, 
					  help="A trainable parameters file, with one gaussian per mixture")
	
	parser.add_argument( "-m", "--maxStates",type=float, metavar="float", required=True, 
					  help="Cluster until at most maxStates are left")

	parser.add_argument( "-r", "--allowUnitRemoval",type=bool, metavar="float", default=True, 
					  help="Also consider the default phonetic pronunciaion as the default context, and remove the unit if all contexts get merged into the default context")

	parser.add_argument("-u", "--inUnits", type=readableFile, metavar="FILE", required=True, 
					  help="Reads the units from a pickle in FILE. Format is the same as errorDrivenUnits.py finalContextsPickle")
	
	parser.add_argument("-U", "--outUnits", type=str, metavar="FILE", required=True, 
					  help="Writes the resulting units pickle to FILE.  Format is the same as errorDrivenUnits.py finalContextsPickle")

	#parser.add_argument( "-o", "--outTrainedParamFile", type=str, metavar="FILE", required=True, 
	#				  help="The output trainable parameters file, into which the tied triphone mixtured will be written.")
	
	#parser.add_argument("-m", "--inNamedCollectionFile", type=readableFile, metavar="FILE", required=True, 
	#				  help="The output namedCollection file, into which the gm names will be written in the correct order.")

	#parser.add_argument("-n", "--outNamedCollectionFile", type=str, metavar="FILE", required=True, 
	#				  help="The output namedCollection file, into which the gm names will be written in the correct order.")

	#parser.add_argument("-u", "--unitStateFile", type=readableFile, metavar="FILE", required=True, 
	#				  help="The units states file used to define the component phones for the units.  The new units start and end with a . (eg .IH_N.) and they must follow the original phonetic units")

	parser.add_argument("-v", "--verbosity", type=int, metavar="INT", default = 51-logging.INFO,
		help="Prints debug info to STDOUT ranges from 1 (critical) to 50 (everything) default: %(default)s")
	return parser

class TiableUnit(dict):
	'''A class for tracking the tying of units.  The items are gmtkIO gaussian components for the sub-unit states.'''
	def __init__(self,context, unit, tokCount):
		self.tokCount=tokCount
		self.active=True #becomes inactive when merged into some other unit
		self.tieTarget=None
		self.unit=unit
		self.context=context	


	def weightedAverageUpdate(self,other):
		'''Set self to be the weighted mean of self and other with mean weights being the token counts.
		   the variance is the element-wise maximum of self and other.
		   The tokCount is updated to be the sum of self and other''' 
		A=self
		B=other
		S=len(A) #number of Sub-unit states
		assert S==len(B)
		sNames=A.keys()
		
		#the gmtkio collections mean and covar collections
		means=A.values()[0].parent[MEAN]
		vars=A.values()[0].parent[COVAR]
		
		for k in sNames:
			mu_A=array(means[A[k].meanName])
			mu_B=array(means[B[k].meanName])
			#print means[A[k].meanName]
			#print [i for i in means[A[k].meanName]]
			#print means[B[k].meanName]
			#print [i for i in means[B[k].meanName]]
			means[A[k].meanName][:]=(mu_A*A.tokCount+mu_B*B.tokCount)/(A.tokCount+B.tokCount)
			#print means[A[k].meanName]
			#print [i for i in means[A[k].meanName]]
			sigma_A=array(vars[A[k].varName])
			sigma_B=array(vars[A[k].varName])
			vars[A[k].varName][:]=maximum(sigma_A,sigma_B)
		
		self.tokCount += other.tokCount

	def htkGaussianDistance(self,other):
		'''computes the distance between two gaussians with diagonal covariance, as described HHed function 
		of the HTK book.  See eq. 17.1 in the HHEd Reference section.  We extend the distance metric to multiple states simply 
		by concatinating them'''

		
		S=len(self) #number of Sub-unit states
		assert S==len(other)
		sNames=self.keys()
		
		#the gmtkio collections mean and covar collections
		means=self.values()[0].parent[MEAN]
		vars=self.values()[0].parent[COVAR]
		
		mu_A=array([])
		mu_B=array([])
		sigma_A=array([])
		sigma_B=array([])
		for k in sNames:
			mu_A=append(mu_A,means[self[k].meanName])
			mu_B=append(mu_B,means[other[k].meanName])
			sigma_A=append(sigma_A,vars[self[k].varName])
			sigma_B=append(sigma_B,vars[other[k].varName])
			
		r = ((mu_A-mu_B)**2/(sigma_A*sigma_B)**.5).mean()
		#print r
		#print sigma_A
		#print vstack((mu_A,mu_B,sigma_A,sigma_B)).T
		#sys.exit()
		return float(r)
		
class TieSet(dict):
	'''A set of units which can be tied with each other.  Each unit can have >1 sub-unit'''
	
	def __init__(self):
		self._li=[]
		self._d=None

			
	def insert(self,context,subUnit,comp,unit,tokCount):
		'''@param comp: a gmtk gaussian component
		   @param unit: a Unit object from the finalUnits
		   @param context: the left/right context combined into a string
		'''  
		if not context in self:
			self[context]=TiableUnit(context,unit,tokCount)
			

		assert not subUnit in self[context]
		assert self[context].tokCount == tokCount
		assert self[context].unit == unit
 
		self[context][subUnit]=comp

		#minOCc is not used  now - occupancy is assumed to be high enough
		
	def calcAllPairwiseDist(self):
		'''calculate the pairwise distance between all units in the tie set. Distance is assumed symmetric.
		   Use it to initialize the distance matrix.  This ignores any merging that may have happened.'''
		li=tuple(v for k,v in sorted(self.items())) #name original, tie target if any
		N=len(li)
		d=ones((N,N))*numpy.inf
		
		for i in range(N):
			for j in range(i+1,N):
				d[i,j]=li[i].htkGaussianDistance(li[j])
		#print d
		self._d=d
		self._li=li
		
	def minDist(self):
		'''@return: the minimum pairwise distance between all units in the tie set.'''
		return numpy.min(self._d)
	
	def stateCount(self):
		'''@return: The number of sub-units (states) in each unit in this tieSet''' 
		return len(self.values()[0])
	
	def tieClosestPair(self):
		'''Tie the two nearest units together.  The merged unit is the weighted average of the merged units,  
			and the weights are token counts of the merged units.'''
		mi=argmin(self._d)
		N=len(self._li)		
		i,j=self.ind2sub(N,N,mi)
		(jw,iw)=(self._li[j].tokCount, self._li[i].tokCount)
		s=float(jw+iw)
		jw=jw/s
		iw=iw/s
		logging.info('Merging %d (%s, weight %f) into %d (%s weight %f)',j,self._li[j].context,jw,i,self._li[i].context,iw)
		logging.debug('\n%s',self._d)
		
		#we change ith unit to be the weighed average of i and j, and make j inactive
		self._li[i].weightedAverageUpdate(self._li[j])
		self._li[j].tieTarget=i
		self._li[j].active=False
		
		#if anything has been tied to j, we tie it to i instead
		for u in self._li:
			if u.tieTarget==j:
				u.tieTarget=i
		
		#set the distance for all units in the merged cluster to inf, so they are 
		#not considered for merging again			
		self._d[j,j:]=numpy.inf
		self._d[:j,j]=numpy.inf
		
		#and recompute distances to i, if they are not merged
		for l in range(i):
			if self._d[l,i] != numpy.inf:
				self._d[l,i]=self._li[l].htkGaussianDistance(self._li[i])
		for l in range(i,N):
			if self._d[i,l] != numpy.inf:
				self._d[i,l]=self._li[i].htkGaussianDistance(self._li[l])

		logging.debug('\n%s',self._d)
		
	@staticmethod
	def ind2sub(row,col,ind):
	    """ Converts row,col indices into one index for .flat """
	    i = ind/col
	    j = ind - i* row
	    return i,j    
		
	def depositActiveUnits(self,finalUnits):
		'''deposit active units and their replacements into finalUnits'''
		#0th element of _li is always the default
		for u in self._li:
			if u.context == '!!!default':
				continue
				
			#we have the property that inactive units are may only be tied
			#to lower-index active units in _li 
			if u.unit in finalUnits.contexts:
				raise ValueError('unit %s is already in the units set, tied to another replacement unit: some naming conflict?'%u.unit)
			if u.tieTarget != None and self._li[u.tieTarget].context=='!!!default':
				logging.debug('context %s was tied to default phonetic trascription and will not be used',u.context)
				continue
			
			finalUnits.contexts[u.unit]=u.unit

			if u.active:				
				if u.unit in finalUnits.repl:
					raise ValueError('unit %s is already in the replacements set: some naming conflict?'%u.unit)
				finalUnits.repl[u.unit]=TiedUnit(u.unit)
				finalUnits.repl[u.unit].addSource(u.unit)
				logging.debug('context %s is a replacement target',u.context)
			else:
				finalUnits.repl[self._li[u.tieTarget].unit].addSource(u.unit)
				logging.debug('context %s is tied to %s',u.context, self._li[u.tieTarget].context)
				
			
	
def main(argv):
	cmd=' '.join(argv)
	opts = parser.parse_args()

	#set up logging
	logging.basicConfig(stream=sys.stdout,format='%(levelname)s %(message)s',level=51-opts.verbosity)
	logging.info('Program started on %s as %s', time.ctime(), cmd)
	
	comment= '%'+" generated with command: %s\n\n"%cmd
	
	ws=Workspace()
	ws.readTrainableParamsFile(opts.inTrainedParamFile)
	logging.info('read readTrainableParamsFile %s', opts.inTrainedParamFile)
	gc=ws[GC]

	finalUnits=FinalUnitsPickle.readFromFile(opts.inUnits)
	logging.info('read finalUnits pickle from %s',	opts.inUnits)
	
	tieSets={}
	for con in finalUnits.contexts:
		phoNames=[Phoneme.phonemeIdListToName(plist) for plist in con[0:3]]
		(center,context)=(phoNames[0],'%s.%s'%(phoNames[1],phoNames[2]))
		unit= '.%s.%s.'%(center,context)
		occCount=sum(con.tot)
		#THE GCs should be named like gc:plp_Col:.AA_B_V.EOW.IY_AH_S_L_IY_EOW.:0
		#here we are assuming that the number of subunits is calculated by adding number
		# of components, (--asBefore flag) 
		for subUnit in range(con.replacement.subUnitCount()):
			c='gc:plp_Col:%s:%d'%(unit,subUnit)
			if not center in tieSets:
				tieSets[center]=TieSet()
			tieSets[center].insert(context,int(subUnit),gc[c],con,occCount)
		#do an error check
		subUnit = con.replacement.subUnitCount()
		c='gc:plp_Col:%s:%d'%(unit,subUnit)
		if c in gc:
			raise ValueError('GC has more substates than expected from finalPickle: %s should not be present.'%c)

	#we also add the ___default context, which always be the 0th element in the list.
	#This represents the phonetic spelling of the unit 
	#if ___default is the only remaining unit, then the unit was essentially tied to the phonetic pronunciation
	#since merging is done always by merging the higher index into the lower one, ___default will never be removed
	if opts.allowUnitRemoval: 
		if TiedUnit.subUnitCount !=TiedUnit._subUnitCountSumComponentPhonemes:
			raise ValueError('finalUnits pickle TiedUnit.subUnitCount is not _subUnitCountSumComponentPhonemes')

		for center,ts in tieSets.items():
			componentPhones=center.split('_')
			ssCount=0
			for p in componentPhones:
				for pss in range(3):
					gcName='gc:plp_Col:%s:%d'%(p, pss)
					if gcName in gc:
						logging.debug('inserting %s as %s:%d',gcName,center,ssCount)
						ts.insert('!!!default',ssCount,gc[gcName],None,999999999)
						ssCount +=1
		
						
	#print tieSets
	
	#calc the pairwise distances for all tie sets 
	units = [[-1,u,name] for name,u in tieSets.items()]
	minStates=0
	totStates=0
	for u in units:
		u[1].calcAllPairwiseDist()
		u[0]=u[1].minDist()
		minStates += u[1].stateCount()
		#if len(u[1])>1:
		#	print (u[1].stateCount(),u[1].stateCount() * len(u[1]))
		totStates += u[1].stateCount() * len(u[1])
	logging.info('Number of tie sets: %d', len(units))
	logging.debug([u[2] for u in units])
	
	if opts.allowUnitRemoval: 
		totStates -= minStates
		minStates -= minStates

	
	logging.info('Starting with total states that are candidate for removal: %d', totStates)
	logging.info('Minimum possible number of states: %d', minStates)
	
	maxDist=9999999999
	m = min(units)
	while totStates>opts.maxStates and m[0]<maxDist:
		logging.info('tie set %s has min distance %f',m[2],m[0])
		m[1].tieClosestPair()
		m[0]=m[1].minDist()
		totStates -= m[1].stateCount()
		logging.debug('totStates is now %d', totStates)
		m = min(units)

	logging.info('final totStates is %d', totStates)
	
	logging.info('creating the tied replacement units')
	newFinalUnits = FinalUnitsPickle()
	for u,v in tieSets.items():
		logging.debug('creating tied units for %s',u) 
		v.depositActiveUnits(newFinalUnits)
		
	logging.info('writing new finalUnits pickle to %s',	opts.outUnits)
	newFinalUnits.writeToFile(opts.outUnits)

		
#the parser is used for generating documentation, so create it always, and augment __doc__ with usage info  
#This messes up epydoc a little, but allows us to keep a single version of documentation for all purposes
parser = makeParser()
__doc__ = __doc__.replace("%InsertOptionParserUsage%\n", parser.format_help())

if __name__ == "__main__":
	main(sys.argv)
	logging.info('Program finished on %s', time.ctime())
