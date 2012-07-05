#ifndef GMTK_OBSSTATS_H
#define GMTK_OBSSTATS_H

#include "GMTK_ObservationMatrix.h"
#include "range.h"

void obsStats(FILE *out_fp, ObservationMatrix* obs_mat,Range& srrng, Range& frrng, const char*pr_str, const size_t hist_bins, const bool quiet_mode);

#endif
