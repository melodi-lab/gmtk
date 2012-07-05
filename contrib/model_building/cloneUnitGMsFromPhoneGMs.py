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
import copy
import textwrap


def makeParser():
	usage = """\
		cloneUnitGMsFrom.py [--help] [options] 
		Takes a trainable parameters file, with one gaussian per mixture, and generates a tri-unit 
		(e.g. tri-phone) trainable parameters file,	tied across left and right contexts.
		
		Also generates a named collection and certain kinds of feature value files
		
		The name format for input gaussian mixtures is:
			gm:featureKind:unit:subUnitState e.g.
			gm:checkwinsizeplp:AA:0
		
		The output will be of form 
			gm:plp:AE:AA:AE:0
		"""

	def readableFile(fname):
		if os.access(fname,os.R_OK):
			return fname
		else:
			print "cannot read file %s"%fname
			raise TypeError()

	parser = argparse.ArgumentParser(prog='cloneUnitGMsFrom.py', description=textwrap.dedent(usage),formatter_class=argparse.RawDescriptionHelpFormatter)
	
	parser.add_argument( "-i", "--inTrainedParamFile",type=readableFile, metavar="FILE", required=True, 
					  help="A trainable parameters file, with one gaussian per mixture")
	
	parser.add_argument( "-o", "--outTrainedParamFile", type=str, metavar="FILE", required=True, 
					  help="The output trainable parameters file, into which the tied triphone mixtured will be written.")
	
	parser.add_argument("-n", "--outNamedCollectionFile", type=str, metavar="FILE", required=True, 
					  help="The output namedCollection file, into which the gm names will be written in the correct order.")

	parser.add_argument("-u", "--unitStateFile", type=readableFile, metavar="FILE", required=True, 
					  help="The units states file used to define the component phones for the units.  The new units start and end with a . (eg .IH_N.) and they must follow the original phonetic units")

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

	phonemes={}
	newUnits=[]
	stateTransitionRows={}
	ncws=Workspace()
	nc=NameCollection(ncws,'plp_Col')
	pref={	DPMF:'mx:plp_Col:',
			MEAN:'mean:plp_Col:',
			COVAR:'covar:plp_Col:',
			GC:'gc:plp_Col:',
			MG:'gm:plp:'}

	s=open(opts.unitStateFile)
	for l in s:
		wl=l.strip().split()
		unit=wl[0]
		count=int(wl[1])
		defn=wl[2:]
		nc.extend(["%s%s:%d"%(pref[MG],unit,i) for i in range(count)])
		if not (unit[0]=='.' and unit[-1]=='.'):
			phonemes[unit]=count
			for i in range(count):
				stateTransitionRows["%s%s:%d"%(pref[MG],unit,i)]=len(stateTransitionRows)
		else:
			newUnits.append((unit,defn))		
	s.close()
	
	#print phonemes
	#print newUnits
	#print sorted(stateTransitionRows.iteritems(), key=operator.itemgetter(1))

	transMat=ws.objects[DCPT]['subPhoneTransitionMatrix']
	for u in newUnits:
		subUCount=0
		for d in u[1]:
			for i in range(phonemes[d]):
				sourceSuff='%s:%d'%(d,i)
				destSuff='%s:%d'%(u[0],subUCount)
				sd={}
				for k in pref:
					sd[k]=(pref[k]+sourceSuff,pref[k]+destSuff)
				subUCount +=1
					
				#now make the actual copies - assuming one gaussian per mixture
				for o in [DPMF,MEAN,COVAR]:
					ws[o][sd[o][0]].copy(ws,sd[o][1])
				GC(ws,sd[GC][1],ws[GC][sd[GC][0]].dimensionality,sd[MEAN][1],sd[COVAR][1])					
				MG(ws, sd[MG][1], ws[MG][sd[MG][0]].dimensionality, sd[DPMF][1], [sd[GC][1]])

				#also add the subPhoneTransitionMatrix row
				cards=transMat.parentCards[0]
				#print cards
				transMat.changeParentCards(cards+1)
				mgName=sd[MG][0]
				#print stateTransitionRows[mgName]
				transMat[:cards]=transMat[:stateTransitionRows[mgName]]
	o=openIOForWrite(opts.outTrainedParamFile, comment)
	ws.writeTrainableParamsIO(o)
	o.close()
	
	o=openIOForWrite(opts.outNamedCollectionFile, comment)
	ncws.writeToIO(NameCollection, o)
	o.close()

	logging.info('wrote  trainable params file %s ', opts.outTrainedParamFile)
	

	

def openIOForWrite(fileName,comment):
	f=file(fileName,'w')
	wsIO=WorkspaceIO(f)
	wsIO.writelnWord(comment)	
	return wsIO

		
#the parser is used for generating documentation, so create it always, and augment __doc__ with usage info  
#This messes up epydoc a little, but allows us to keep a single version of documentation for all purposes
parser = makeParser()
__doc__ = __doc__.replace("%InsertOptionParserUsage%\n", parser.format_help())

if __name__ == "__main__":
	main(sys.argv)
	logging.info('Program finished on %s', time.ctime())
