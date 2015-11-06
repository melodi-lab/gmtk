/*-
 * GMTK_MaxClique.cc
 *
 *     maxClique support and JT probabilistic inference. Includes the
 *     implementation of a message in the message passing algorithm.
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */




/*
 * TODO: turn this into multiple files
 *   mc, mctable CE, mctable DE, mctable prune, csctable
 *
 */



/*
 * TODO:
 *   1) Maybe: In some of the structures the choice of when
 *   to use a hash table vs. an inlined value is done
 *   by either using an unsigned* or a unsigned. Create
 *   options where we do 'unsigned val[LEN]' and where
 *   only if packed value is > LEN words do we resort to
 *   the hash table.
 *
 *   2) Figure out way to remove the (unsigned**) type casts
 *   in front of call to packer object (see below).   
 *
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <ctype.h>

#include <iterator>
#include <map>
#include <set>
#include <algorithm>

#include "general.h"
#include "error.h"
#include "debug.h"
#include "rand.h"

#include "GMTK_FileParser.h"
#include "GMTK_RV.h"
#include "GMTK_DiscRV.h"
#include "GMTK_HidDiscRV.h"
#include "GMTK_ObsDiscRV.h"
#include "GMTK_ContRV.h"
#include "GMTK_GMTemplate.h"
#include "GMTK_GMParms.h"
#include "GMTK_MaxClique.h"
#include "GMTK_MDCPT.h"
#include "GMTK_MSCPT.h"
#include "GMTK_MTCPT.h"
#include "GMTK_Mixture.h"
#include "GMTK_ObservationSource.h"
#include "GMTK_JunctionTree.h"

VCID(HGID)



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        FactorClique support
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//        Constructors
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


FactorClique::FactorClique(FactorInfo& _factorInfo,
			   vector <RV*>& unrolled_rvs,
			   map < RVInfo::rvParent, unsigned > ppf,
			   const unsigned offset)
{
  nodes = getRVVec(unrolled_rvs,ppf,
		   _factorInfo.variables,
		   _factorInfo.frame + offset);

  orderedNodes = getRVOVec(unrolled_rvs,ppf,
			   _factorInfo.variables,
			   _factorInfo.frame + offset);

  factorInfo = &_factorInfo;

  assert( orderedNodes.size() == nodes.size() );

}



/////////////////////////////////////////////	
/// END OF FILE
/////////////////////////////////////////////
