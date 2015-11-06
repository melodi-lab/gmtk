/*-
 * GMTK_VHashMapUnsignedUnsignedKeyUpdatable.cc
 *
 * Written by Jeff Bilmes <bilmes@ee.washington.edu>
 *
 * Copyright (C) 2003 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 *
 */




#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#if HAVE_HG_H
#include "hgstamp.h"
#endif

#include "general.h"
#include "GMTK_VHashMapUnsignedUnsignedKeyUpdatable.h"

VCID(HGID)



VHashMapUnsignedUnsignedKeyUpdatable::
VHashMapUnsignedUnsignedKeyUpdatable(const unsigned arg_vsize,
				     unsigned approximateStartingSize)
  : VHashMapUnsignedUnsigned(arg_vsize,approximateStartingSize)
{
  // do nothing
}

