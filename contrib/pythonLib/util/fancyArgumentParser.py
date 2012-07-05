"""

@author: 	Arthur Kantor
@contact: 	akantorREMOVE_THIS@uiuc.edu
@copyright:	Arthur Kantor 2010
@license: 	GPL version 3
@date: 		1/6/2010
@version: 	0.1

"""

import argparse
from argparse import *
import textwrap
import os


class FancyArgumentParser(argparse.ArgumentParser):
	'''formats the usage nicely for command line and includes some extra parameter validation choices'''
	def __init__(self,**kw):
		kw['description']=textwrap.dedent(kw['description'])
		kw['formatter_class']=argparse.RawDescriptionHelpFormatter
		parser = argparse.ArgumentParser.__init__(self,**kw)

	'''this function can be used as type to make sure that strings are readable filenames'''
	@classmethod
	def readableFile(cls, fname):
		#print "testing %s\n"%fname
		if os.access(fname,os.R_OK):
			return fname
		else:
			raise ValueError("cannot read file %s"%fname)


