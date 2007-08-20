/*  Generated header
 *  File Name : GMTK_ObsStats.cc
 *
 *  Created   : 2003-12-05 14:26:06 karim
 *  Author    : Karim Filali (karim@cs.washington.edu)
 *  Time-stamp: <>
*/


#include <limits.h>
#include <float.h>
#include <math.h>
#include "GMTK_ObsStats.h"

typedef struct { 
  size_t sent_no;
  size_t frame_no;
} PfileLocation;



void obsStats(FILE *out_fp, ObservationMatrix* obs_mat,Range& srrng, Range& frrng, const char*pr_str, const size_t hist_bins, const bool quiet_mode) {

    // Feature and label buffers are dynamically grown as needed.

    size_t buf_size = 300;      // Start with storage for 300 frames.
    const size_t n_labs = obs_mat->numDiscrete();
    const size_t n_ftrs = obs_mat->numContinuous();

    float *ftr_buf = new float[buf_size * n_ftrs];
    UInt32* lab_buf = new UInt32[buf_size * n_labs];

    size_t total_frames = 0;
    double *const ftr_sum = new double [frrng.length()];
    double *const ftr_sumsq = new double [frrng.length()];
    double *const ftr_means = new double [frrng.length()];
    double *const ftr_stds = new double [frrng.length()];
    float *const ftr_maxs = new float [frrng.length()];
    PfileLocation *const ftr_maxs_locs = new PfileLocation [frrng.length()];
    float *const ftr_mins = new float [frrng.length()];
    PfileLocation *const ftr_mins_locs = new PfileLocation [frrng.length()];
    size_t * histogram = NULL;

    double *ftr_sum_p;
    double *ftr_sumsq_p;
    double *ftr_means_p;
    double *ftr_stds_p;
    float *ftr_maxs_p;
    float *ftr_mins_p;
    PfileLocation *ftr_maxs_locs_p;
    PfileLocation *ftr_mins_locs_p;


    // Initialize the above declared arrays
    size_t i,j;
    for (i=0;i<frrng.length();i++) {
      ftr_sum[i] = ftr_sumsq[i] = 0.0;
      ftr_means[i] = ftr_stds[i] = 0.0;
      ftr_maxs[i] = -FLT_MAX;
      ftr_mins[i] = FLT_MAX;
    }

    for (Range::iterator srit=srrng.begin();!srit.at_end();srit++) {
      obs_mat->loadSegment(*srit);
      const size_t n_frames = obs_mat->numFrames();

	if (!quiet_mode) {
	  if ((*srit) % 100 == 0)
	    printf("Processing sentence %d\n",(*srit));
	}

        // Increase size of buffers if needed.
        if (n_frames > buf_size)
        {
            // Free old buffers.
            delete ftr_buf;
            delete lab_buf;

            // Make twice as big to cut down on future reallocs.
            buf_size = n_frames * 2;

            // Allocate new larger buffers.
            ftr_buf = new float[buf_size * n_ftrs];
            lab_buf = new UInt32[buf_size * n_labs];
        }

	for(unsigned frame_no = 0;  frame_no < n_frames; ++frame_no) {
	  float* start_of_frame = obs_mat->floatVecAtFrame(frame_no);
	  for(unsigned feat_no = 0;  feat_no < n_ftrs; ++feat_no) {
	    ftr_buf[frame_no * n_ftrs + feat_no] = *(start_of_frame  + feat_no);
	  }
	}

	Range prrng(pr_str,0,n_frames);
	for (Range::iterator prit=prrng.begin();
	     !prit.at_end() ; ++prit) {

	  const float *const ftr_buf_p = ftr_buf + (*prit)*n_ftrs;
	  ftr_sum_p = ftr_sum;
	  ftr_sumsq_p = ftr_sumsq;
	  ftr_maxs_p = ftr_maxs;
	  ftr_mins_p = ftr_mins;
	  ftr_maxs_locs_p = ftr_maxs_locs;
	  ftr_mins_locs_p = ftr_mins_locs;
	  

	  for (Range::iterator frit=frrng.begin();
	       !frit.at_end(); ++frit) {
	    const double val = ftr_buf_p[*frit];
	    *ftr_sum_p++ += val;
	    *ftr_sumsq_p++ += (val*val);
	    if (val > *ftr_maxs_p) {
	      *ftr_maxs_p = val;
	      ftr_maxs_locs_p->sent_no = (*srit);
	      ftr_maxs_locs_p->frame_no = (*prit);
	    }
	    else if (val < *ftr_mins_p) {
	      *ftr_mins_p = val;
	      ftr_mins_locs_p->sent_no = (*srit);
	      ftr_mins_locs_p->frame_no = (*prit);
	    }
	    ftr_maxs_p++;
	    ftr_mins_p++;
	    ftr_maxs_locs_p++;
	    ftr_mins_locs_p++;
	  }
	}
	total_frames += prrng.length();
    }

    if (total_frames == 1) {
      if (!quiet_mode) {
	  printf("WARNING:: Ranges specify using only one frame for statistics.\n");
      }
    }


    ftr_means_p = ftr_means;
    ftr_stds_p = ftr_stds;
    ftr_sum_p = ftr_sum;
    ftr_sumsq_p = ftr_sumsq;
    for (i=0;i<frrng.length();i++) {
      (*ftr_means_p) = (*ftr_sum_p)/total_frames;
      (*ftr_stds_p) = 
	sqrt(
	     ((*ftr_sumsq_p) - (*ftr_sum_p)*(*ftr_sum_p)/total_frames)/
	     total_frames);
      ftr_means_p++;
      ftr_stds_p++;
      ftr_sum_p++;
      ftr_sumsq_p++;
    }

    if (hist_bins > 0) {
      //  go through and do a second pass on the file.
      size_t *hist_p;
      float* ftr_ranges = new float[frrng.length()];
      histogram = new size_t [frrng.length()*hist_bins];
      ::memset(histogram,0,sizeof(size_t)*frrng.length()*hist_bins);
      for (i=0;i<frrng.length();i++)
	ftr_ranges[i] = ftr_maxs[i]-ftr_mins[i];

      if (!quiet_mode) {
	printf("Computing histograms..\n");
      }
      for (Range::iterator srit=srrng.begin();!srit.at_end();srit++) {
	obs_mat->loadSegment(*srit);
	const size_t n_frames = obs_mat->numFrames();
	
	if (!quiet_mode) {
	    if ((*srit) % 100 == 0)
	      printf("Processing sentence %d\n",(*srit));
	  }

	  for(unsigned frame_no = 0;  frame_no < n_frames; ++frame_no) {
	    float* start_of_frame = obs_mat->floatVecAtFrame(frame_no);
	    for(unsigned feat_no = 0;  feat_no < n_ftrs; ++feat_no) {
	      ftr_buf[frame_no * n_ftrs + feat_no] = *(start_of_frame  + feat_no);
	    }
	  }

	  Range prrng(pr_str,0,n_frames);
	  for (Range::iterator prit=prrng.begin();
	       !prit.at_end() ; ++prit) {

	    const float *const ftr_buf_p = ftr_buf + (*prit)*n_ftrs;
	    hist_p = histogram;
	    ftr_maxs_p = ftr_maxs;
	    ftr_mins_p = ftr_mins;
	    float *ftr_ranges_p = ftr_ranges;
	    
	    for (Range::iterator frit=frrng.begin();
		 !frit.at_end(); ++frit) {
	      const double val = ftr_buf_p[*frit];
	      size_t ind = size_t(
				  hist_bins*0.9999*
				  (val-*ftr_mins_p)/(*ftr_ranges_p));
	      hist_p[ind]++;
	      hist_p += hist_bins;
	      ftr_maxs_p++;
	      ftr_mins_p++;
	      ftr_ranges_p++;
	    }
	  }
      }
      delete ftr_ranges;
    }

    ftr_means_p = ftr_means;
    ftr_stds_p = ftr_stds;    
    ftr_sum_p = ftr_sum;
    ftr_sumsq_p = ftr_sumsq;
    ftr_maxs_p = ftr_maxs;
    ftr_mins_p = ftr_mins;
    ftr_maxs_locs_p = ftr_maxs_locs;
    ftr_mins_locs_p = ftr_mins_locs;

    double max_maxs_stds=-FLT_MAX;
    double min_mins_stds=+FLT_MAX;
    size_t *hist_p = histogram;
    for (i=0;i<frrng.length();i++) {
      const double maxs_stds = (*ftr_maxs_p)/(*ftr_stds_p);
      const double mins_stds = (*ftr_mins_p)/(*ftr_stds_p);
      fprintf(out_fp,"%d %f %f %f %d %d %f %d %d %f %f ",i,
	     *ftr_means_p,*ftr_stds_p,
	     *ftr_maxs_p,ftr_maxs_locs_p->sent_no,ftr_maxs_locs_p->frame_no,
	     *ftr_mins_p,ftr_mins_locs_p->sent_no,ftr_mins_locs_p->frame_no,
	     maxs_stds,mins_stds);
      if (hist_bins > 0) {
	for (j=0;j<hist_bins;j++) {
	  fprintf(out_fp,"%d ",*hist_p++);
	}
      }
      fprintf(out_fp,"\n");

      ftr_means_p++;
      ftr_stds_p++;
      ftr_maxs_p++;
      ftr_mins_p++;
      ftr_maxs_locs_p++;
      ftr_mins_locs_p++;
      if (maxs_stds > max_maxs_stds)
	max_maxs_stds = maxs_stds;
      if (mins_stds < min_mins_stds)
	min_mins_stds = mins_stds;
    }
    if (!quiet_mode) {
      printf("total sents used = %d, total frames used = %d\n",
	      srrng.length(),total_frames);
      printf("max_maxs_stds = %f, min_mins_stds = %f\n",
	     max_maxs_stds,min_mins_stds);
    }

    delete ftr_buf;
    delete lab_buf;

    delete ftr_sum;
    delete ftr_sumsq;    
    delete ftr_means;
    delete ftr_stds;
    delete ftr_maxs;
    delete ftr_maxs_locs;
    delete ftr_mins;
    delete ftr_mins_locs;

}

