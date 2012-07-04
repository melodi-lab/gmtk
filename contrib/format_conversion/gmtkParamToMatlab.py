#!/usr/bin/python
"""
%InsertOptionParserUsage%


@author: 	Arthur Kantor
@contact: 	akantorREMOVE_THIS@uiuc.edu
@copyright:	Arthur Kantor 2008
@license: 	GPL version 3
@date: 		1/10/2010
@version: 	0.1

"""

import sys, re, os;
import logging
import time
import argparse
import textwrap
from gmtkParam import *
#from util import *
#from util.UniqueList import *
#import copy



def makeParser():
	usage = """\
		gmtkParamToMatlab.py [--help] [options] 
		Takes a trainable parameters file with N-dimensional mixture gaussians with diagonal covariances and writes a M x 2*N+2 matrix of the form
		
		MixtureNumber ComponentNumber means covariances
		...
		
		Any tied parameters are untied.
		The order of the mixtures is given by the named collection in file nameCollectionFile.
		
		"""

	def readableFile(fname):
		if os.access(fname,os.R_OK):
			return fname
		else:
			print "cannot read file %s"%fname
			raise TypeError()

	parser = argparse.ArgumentParser(prog='gmtkParamToMatlab.py', description=textwrap.dedent(usage),formatter_class=argparse.RawDescriptionHelpFormatter)
	
	parser.add_argument( "-i", "--inTrainedParamFile",type=readableFile, metavar="FILE", required=True, 
					  help="An ascii trainable parameters file")
	
	parser.add_argument("-o", "--outMatlabAsciiFile", type=str, metavar="FILE", required=True, 
					  help="Writes the mixtures to FILE.")

	parser.add_argument("-r", "--orderNameCollectionFile",type=readableFile, metavar="FILE", required=True, 
					  help="the named collection of gaussian mixtures specifying the order in which they are written out")

	parser.add_argument("-v", "--verbosity", type=int, metavar="INT", default = 51-logging.INFO,
		help="Prints debug info to STDOUT ranges from 1 (critical) to 50 (everything) default: %(default)s")
	return parser


	
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
	gm=ws[MG]

	ws2=Workspace()
	ws2.readFromFile(NameCollection,opts.orderNameCollectionFile)
	logging.info('read NameCollection %s', opts.orderNameCollectionFile)
	orderColl=ws2[NameCollection].values()[0] #there is only one value

	matf = file(opts.outMatlabAsciiFile,'w')

	for (mixNum,mixName) in enumerate(orderColl):
		#print mixName+"\n"
		mix=gm[mixName]
		for (w,gcName) in zip(ws[DPMF][mix.weightsDpmfName], mix.gcNames):
			gc=ws[GC][gcName]
			(mean,covar) = (ws[MEAN][gc.meanName],ws[COVAR][gc.varName])
			matf.write(' '.join([`mixNum`, `w`, ' '.join(`f` for f in mean), ' '.join(`f` for f in covar), '\n']))
			
	matf.close()

	logging.info('wrote %d mixtures to outMatlabAsciiFile %s', mixNum+1,opts.outMatlabAsciiFile)

		

#the parser is used for generating documentation, so create it always, and augment __doc__ with usage info  
#This messes up epydoc a little, but allows us to keep a single version of documentation for all purposes
parser = makeParser()
__doc__ = __doc__.replace("%InsertOptionParserUsage%\n", parser.format_help())

if __name__ == "__main__":
	main(sys.argv)
	logging.info('Program finished on %s', time.ctime())
