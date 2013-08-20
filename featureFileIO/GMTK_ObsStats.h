#ifndef GMTK_OBSSTATS_H
#define GMTK_OBSSTATS_H

/*
 * Copyright (C) 2004 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 */

#if 0
#  include "GMTK_ObservationMatrix.h"
#else
#  include "GMTK_FileSource.h"
#endif

#include "range.h"

void obsStats(FILE *out_fp, FileSource* obs_mat,Range& srrng, Range& frrng, const char*pr_str, const size_t hist_bins, const bool quiet_mode);

#endif
