#!/usr/bin/python
# -*- coding: utf-8 -*-
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
import argparse
import textwrap
import logging
import time

from gmtkParam.obsIO import *
def usage():
 return """
"""


def makeParser():
	usage = """\
		collectionToPhysicalState.py [--help] [options] 

		Replace the RV values in the observation file with the new specified values.
		Useful for replacing obeserved triphone state with the actual physical state.
		"""

	def readableFile(fname):
		if os.access(fname,os.R_OK):
			return fname
		else:
			print "cannot read file %s"%fname
			raise TypeError()

	parser = argparse.ArgumentParser(prog='collectionToPhysicalState.py', description=textwrap.dedent(usage),formatter_class=argparse.RawDescriptionHelpFormatter)
	
	parser.add_argument("-s", "--sourceValues", type=readableFile, metavar="FILE", required=True)
	parser.add_argument("-t", "--targetValues", type=readableFile, metavar="FILE", required=True,
		help="Ordered list of source and target RV values (GM names) from FILE.")

	parser.add_argument("-I", "--inObservations", type=readableFile, metavar="FILE", required=True)
	parser.add_argument("-O", "--outObservations", type=str, metavar="FILE", required=True,
		help="The input and output observation data in feacat ASCII format.")

	parser.add_argument("-c", "--column", type=int, required=True, default=0,
		help="on which column to perform value replacement (0-based, not counting the utterance and frame Ids)")
	
	parser.add_argument("-v", "--verbosity", type=int, metavar="INT", default = 51-logging.INFO,
		help="Prints debug info to STDERR ranges from 1 (critical) to 50 (everything) default: %(default)s")
	return parser


def readValueList(fname):
	f=open(fname)
	dic = dict(enumerate([s.strip() for s in f]))
	logging.info('read %d values from %s'% (len(dic),fname))
	f.close()
	return dic

def main(argv):
	cmd=' '.join(argv)
	opts = parser.parse_args()

	#set up logging
	logging.basicConfig(stream=sys.stderr,format='%(levelname)s %(message)s',level=51-opts.verbosity)
	logging.info('Program started on %s as %s', time.ctime(), cmd)
 
	sourceDic=readValueList(opts.sourceValues)
	targDic=readValueList(opts.targetValues)
	targDic=dict([(v,str(k)) for (k,v) in targDic.items()])
	#print targDic.items()[0]
	#print sourceDic.items()[0]
	valMap=dict([(s,targDic[sName]) for (s,sName) in sourceDic.items()])
	logging.info('valMap has %d items'% (len(valMap),))

	ifile = open(opts.inObservations)
	oSource = obsSource(ifile)
	of = open(opts.outObservations,'w')
	obsSink = ObsSink(of)
	for utt in obsSource(ifile):
		for (i,frame) in enumerate(utt):
			els = frame[2].split()
			els[opts.column] = valMap[int(float(els[opts.column]))]
			utt[i]=(frame[0],frame[1],' '.join(els))
		obsSink.writeObs(utt)
	of.close()
	ifile.close()
	logging.info('Processed %d frames', obsSink.uttCount)

	sys.exit()


#the parser is used for generating documentation, so create it always, and augment __doc__ with usage info  
#This messes up epydoc a little, but allows us to keep a single version of documentation for all purposes
parser = makeParser()
__doc__ = __doc__.replace("%InsertOptionParserUsage%\n", parser.format_help())

if __name__ == "__main__":
	main(sys.argv)

