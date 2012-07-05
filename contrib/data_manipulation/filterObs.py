#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
%InsertOptionParserUsage%


@author: 	Arthur Kantor
@contact: 	akantorREMOVE_THIS@uiuc.edu
@copyright:	Arthur Kantor 2009
@license: 	GPL version 3
@date: 		7/8/2009
@version: 	0.1

"""

import sys, re, os;
import logging
import time
import argparse
import operator 

from util import *
from util.UniqueList import *
import copy
import textwrap

from asr.orderedDict import OrderedDict
from asr.UnitSet import UnitSet
from gmtkParam.obsIO import *


def makeParser():
	usage = """\
		filterObs.py [--help] [options] 

		Filter a set of Utterance (both the accoustic observations and their transcriptions, 
		according to specified criteria.
		
		obs file is in feacat ascii format
		transcription file is in one transcription per line, 0th line corresponding to 0th observation utterance

		For now, the only filtering allowed is consistUttDur which allows through only those utterances with enough frames
		for all the words specified in the utterance
		"""

	def readableFile(fname):
		if os.access(fname,os.R_OK):
			return fname
		else:
			print "cannot read file %s"%fname
			raise TypeError()

	parser = argparse.ArgumentParser(prog='filterObs.py', description=textwrap.dedent(usage),formatter_class=argparse.RawDescriptionHelpFormatter)
	
	parser.add_argument("-d", "--dict", type=readableFile, metavar="FILE", required=True,
		help="Reads the phonetic dictionary from FILE.")

	parser.add_argument("-p", "--inSubUnitState", type=readableFile, metavar="FILE", required=True,
		help="read initial sub-unit states from FILE. ")

	parser.add_argument("-r", "--outDurations", type=str, metavar="FILE", required=False,
		help="if given, writes the minimum word durations to FILE for words ordered as in --dict file")
		
	#parser.add_argument("-w", "--wordIdCol", type=int, metavar="INT", required=True,
	#	help="The 0-based column in inObservations specifying the word Id.")



	parser.add_argument("-i", "--inTranscription", type=readableFile, metavar="FILE", required=False,
		help="The input transcription file, one line per utterance.")
	parser.add_argument("-o", "--outTranscription", type=str, metavar="FILE", required=False,
		help="Where to write transcriptions after filtering.")

	parser.add_argument("-I", "--inObservations", type=readableFile, metavar="FILE",
		help="The input observation data in feacat ASCII format. Defaults to STDIN")
	parser.add_argument("-O", "--outObservations", type=str, metavar="FILE",
		help="Where to write transcriptions after filtering in feacat ASCII format. Defaults to STDOUT.")

	parser.add_argument("-f", "--filterType",choices=['consistentUttDur'], metavar="filter",
					default=('consistentUttDur',), action='append', type=str,
					help="Multiple filters can be given, and are applied in order. See the filter types described above. Default: consistUttDur")

	parser.add_argument("-v", "--verbosity", type=int, metavar="INT", default = 51-logging.INFO,
		help="Prints debug info to STDERR ranges from 1 (critical) to 50 (everything) default: %(default)s")
	return parser


def main(argv):
	cmd=' '.join(argv)
	opts = parser.parse_args()

	#set up logging
	logging.basicConfig(stream=sys.stderr,format='%(levelname)s %(message)s',level=51-opts.verbosity)
	logging.info('Program started on %s as %s', time.ctime(), cmd)

	unitSet=UnitSet.readUnitStates(opts.inSubUnitState)
	dict=OrderedDict.readDict(opts.dict)
		
	wordDurations={}
	for (w, defn) in dict:
		dur=sum([unitSet[unitName][2] for unitName in defn])
		wordDurations[w]=dur
		#print (w,defn,dur)
	logging.info('got word durations')

	if opts.outDurations:
		outDur =open(opts.outDurations,'w')
		for (w, defn) in dict:
			outDur.write(str(wordDurations[w])+'\n')
		outDur.close()
		logging.info('wrote word durations')
		return
		
	transF=file(opts.inTranscription)
	if opts.inObservations:	
		obsF=file(opts.inObservations)
	else:
		obsF=sys.stdin
	obsGen=obsSource(obsF)
	
	logging.info('construct the chain of filters')
	filters=[]
	for fstr in opts.filterType:
		if fstr == 'consistentUttDur':
			filters.append(ConsistentUttDurationFilter(wordDurations, -1))
		else:
			raise ValueError ('Unknown filter %s'%fstr)

	transSink = TransSink(open(opts.outTranscription,'w'))
	if opts.outObservations:
		obsSink = ObsSink(open(opts.outObservations,'w'))
	else:
		obsSink = ObsSink(sys.stdout)
	
	trCounter=0
	for tr in transSource(transF):
		try:
			uttObs= obsGen.next()
		except StopIteration:
			raise ValueError('No observations remain after processing the first %d transcriptions.'%trCounter )
		trCounter += 1
		
		#now we have the current utterance observation and transcription
		
		for f in filters:
			passes =f.filterUtt(uttObs,tr)
			if not passes:
				#don't apply later filters if an earlier filter didn't pass an utterance
				break
		if passes:
			transSink.writeTrans(tr)
			obsSink.writeObs(uttObs)
	try:
		obsGen.next()
		raise ValueError('All %d transcriptions have been processed, while Observations remain.'%trCounter )
	except StopIteration:
		pass  #this is normal - no observations left

	logging.info('finished processing utterances')
	for (i,fstr) in enumerate(opts.filterType):
		logging.info('%s: %s'%(fstr,filters[i].stats()))
	
def transSource(transF):
	for t in transF:
		trans=t.split()
		yield trans

class TransSink(object):
	'''writes transcriptions to a file'''
	def __init__(self,ofile):
		self._f = ofile
		
	def writeTrans(self, trans):
		self._f.write(' '.join(trans))
		self._f.write('\n')
		

#various filters go here

class Filter(object):
	def __init__(self):
		'''base filter class - just keeps track of utterances and passes-through everything'''
		self._uttCount=0
		self._emitCount=0

	def filterUtt(self,obs,trans):
		self._uttCount += 1
		self._emitCount +=1
		return True
		
	def stats(self):
		return 'Saw %d utterances, emitted  %d utterances'%(self._uttCount,self._emitCount)

class ConsistentUttDurationFilter(Filter):
	def __init__(self,wordDurations, wordIdCol):
		super(ConsistentUttDurationFilter,self).__init__()
		self._wd=wordDurations
		self._wordIdCol=wordIdCol

	def filterUtt(self,obs,trans):
		self._uttCount += 1
		
		minFramesReq=sum(self._wd[word] for word in trans[1:] )
		if len(obs)>=minFramesReq:
			self._emitCount +=1
			return True
		else:
			return False
		
#the parser is used for generating documentation, so create it always, and augment __doc__ with usage info  
#This messes up epydoc a little, but allows us to keep a single version of documentation for all purposes
parser = makeParser()
__doc__ = __doc__.replace("%InsertOptionParserUsage%\n", parser.format_help())

if __name__ == "__main__":
	main(sys.argv)
	logging.info('Program finished on %s', time.ctime())
