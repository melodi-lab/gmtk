/*  Generated header
 *  File Name : GMTK_ObsInitmg.cc
 *
 *  Created   : 2003-12-15 20:16:42 karim
 *  Author    : Karim Filali (karim@cs.washington.edu)
 *  Time-stamp: <>

 Create an initializtion .mg file for pfile_mi.cc using a simple
 k-means algorithm.  Adapted from the pfile_initmg by Jeff Bilmes
 <bilmes@ee.washington.edu>

*/

#include <cstdio>
#include <cerrno>
#include <cstring>
#include <limits.h>
#include <float.h>
#include <cmath>
#include <cassert>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>

#include "range.h"
#include "rand.h"
//#include "spi.h"

#include "GMTK_ObservationMatrix.h"

#include "pfile.h"
//#include "parse_subset.h"
#include "error.h"
#include "arguments.h"
#include "GMTK_WordOrganization.h"

#include "GMTK_Kmeans.h"


void initmg(ObservationMatrix* obs_mat,
	    FILE *out_fp,
	    Range& srrng,
	    Range& cfrrng,
	    const char *pr_str,
	    const int num_clusters,
	    const int maxReInits,
	    const int minSamplesPerCluster,
	    const int numRandomReStarts,
	    const bool quiet_mode,
	    const float conv_thres)
{

    // Feature and label buffers are dynamically grown as needed.

    const int n_ftrs = (int)obs_mat->numContinuous();

    float *ftr_buf = NULL;

    kmeans::kmeans_k = num_clusters;
    kmeans::kmeans_vl = 2;

    const int num_kmeans = n_ftrs*(n_ftrs-1)/2;
    kmeans *kms = new kmeans[num_kmeans];
    size_t sent_no;

    double bestVarianceSum=DBL_MAX;
    
    for (int epoch=0;epoch<numRandomReStarts;epoch++) {

      printf("Epoch %d\n",epoch);

      int i,j,iter=0;
      int reInits = 0;
    
      float max_dist = 0;

      for (i=0;i<num_kmeans;i++) {
	kms[i].randomAssignment = true;
	kms[i].done = false;
      }

      do {
	iter++;

	for (i=0;i<num_kmeans;i++) {
	  if (!kms[i].done)
	    kms[i].initNew();
	}

	for (Range::iterator srit=srrng.begin();!srit.at_end();srit++) {
	  sent_no = (*srit);

	  obs_mat->loadSegment(sent_no);
	  const size_t n_frames = obs_mat->numFrames();

	  if (!quiet_mode) {
	    if ((sent_no) % 1000 == 0)
	      printf("Processing sentence %ld\n",(sent_no));
	  }

	  //const int n_read =  in_streamp->read_ftrslabs(sent_no, ftr_buf, lab_buf);
	  for(unsigned frame_no = 0;  frame_no < n_frames; ++frame_no) {
	    const float* start_of_frame = obs_mat->floatVecAtFrame(frame_no);
	    for(int feat_no = 0;  feat_no < n_ftrs; ++feat_no) {
	      ftr_buf[frame_no * n_ftrs + feat_no] = *(start_of_frame  + feat_no);
	    }
	  }
	  

	  Range prrng(pr_str,0,n_frames);

	  for (Range::iterator prit=prrng.begin();
	       !prit.at_end() ; ++prit) {
	    const float *const ftr_buf_p = ftr_buf + (*prit)*n_ftrs;
	    
	    int kmno = 0;
	    float buf2[2];
	    for (i=0;i<n_ftrs;i++) {
	      buf2[0] = ftr_buf_p[i];
	      for (j=i+1;j<n_ftrs;j++) {
		if (!kms[kmno].done) {
		  buf2[1] = ftr_buf_p[j];
		  if (kms[kmno].randomAssignment)
		    kms[kmno].add2newRand(buf2);
		  else 
		    kms[kmno].add2new(buf2);
		}
		kmno++;
	      }
	    }
	  }

	}

	for (i=0;i<num_kmeans;i++) {
	  if (!kms[i].done)
	    kms[i].randomAssignment = false;
	}

	// make sure each kmeans had some data
	max_dist = 0;
	int num_active=0;
	int num_reinits = 0;
	for (i=0;i<num_kmeans;i++) {
	  if (!kms[i].done) {
	    if (kms[i].zeroCounts()) {
	      // This shouldn't happen with uniform skmeans. All kmeans
	      // objects should be getting some data.
	      error("kms[%d] was given no input. Probably an command line argument error",i);
	    } else if (kms[i].someClusterHasLessThanNEntries(minSamplesPerCluster)) {
	      // fprintf(stderr,"Warning: kms %d, some clusters have < %d entries. Re-initializing.\n",i,minSamplesPerCluster);
	    
	      // only re-init this particular kmeans object, not all of them.
	      kms[i].randomAssignment = true;
	      num_reinits++;
	      max_dist = 1.0; // keep this condition from stopping the loop.
	      // break;
	    } else {
	      kms[i].finishNew();
	      float tmp = kms[i].newCurDist();
	      if (tmp > max_dist)
		max_dist = tmp;
	      kms[i].swapCurNew();
	      if (tmp <= conv_thres) {
		kms[i].done = true;
		// printf("Iter %d: kms %d converged\n",iter,i);
	      }
	    }
	    num_active++;
	  }
	}
	if (num_reinits > 0)
	  reInits++;
	printf("Iter %d: max_dist = %e, cur num_reinits = %d, tot itr re_init = %d, num_active = %d\n",iter,max_dist,num_reinits,reInits,num_active);
	fflush(stdout);

      } while (max_dist > conv_thres && reInits <= maxReInits);

      if (reInits > maxReInits) {
	error("Error. %d re-inits and convergence didn't occur.", reInits);
      }

      printf("Epoch %d: Computing Variances\n",epoch);

      // Do one more pass over file to compute the variances.
      for (Range::iterator srit=srrng.begin();!srit.at_end();srit++) {
	sent_no = (*srit);

	obs_mat->loadSegment(sent_no);
	const size_t n_frames = obs_mat->numFrames();

	if (!quiet_mode) {
	  if ((sent_no) % 1000 == 0)
	    printf("Processing sentence %ld\n",(sent_no));
	}

	 for(unsigned frame_no = 0;  frame_no < n_frames; ++frame_no) {
	    const float* start_of_frame = obs_mat->floatVecAtFrame(frame_no);
	    for(int feat_no = 0;  feat_no < n_ftrs; ++feat_no) {
	      ftr_buf[frame_no * n_ftrs + feat_no] = *(start_of_frame  + feat_no);
	    }
	  }


	Range prrng(pr_str,0,n_frames);

	for (Range::iterator prit=prrng.begin();
	     !prit.at_end() ; ++prit) {
	  const float *const ftr_buf_p = ftr_buf + (*prit)*n_ftrs;
	  
	  int kmno = 0;
	  float buf2[2];
	  for (i=0;i<n_ftrs;i++) {
	    buf2[0] = ftr_buf_p[i];
	    for (j=i+1;j<n_ftrs;j++) {
	      buf2[1] = ftr_buf_p[j];
	      kms[kmno].computeVariances(buf2);
	      kmno++;
	    }
	  }
	}
      }

    
      double sumVariances=0;
      for (i=0;i<num_kmeans;i++) {
	sumVariances+=kms[i].finishVariances();
      }

      if (sumVariances < bestVarianceSum) {
	bestVarianceSum = sumVariances;
	for (i=0;i<num_kmeans;i++) {
	    kms[i].save();
	}
      }
      
      printf("End of Epoch %d of %d, variance sum = %e, best = %e\n",epoch,
	     numRandomReStarts,
	     sumVariances,bestVarianceSum);

    }
    
//     kmeans *kmsp = kms;
//     for (int i=0;i<n_ftrs;i++) {
//       for (int j=i+1;j<n_ftrs;j++) {
// 	fprintf(out_fp,"feat %d:%d: \n",i,j);
// 	kmsp->printSaved(out_fp);
// 	kmsp++;
//       }
//     }


    // We must have that to > from
#define KMSPOS(from,to,n) ((from)*(n) - ((from)+1)*((from)+2)/2 + (to))

  for (Range::iterator cnit=cfrrng.begin();
       !cnit.at_end();cnit++) {
    const int l = (*cnit);
    // j-loop is over position in the lagged feature vector, x

    for (int j=0;j<n_ftrs;j++) {
      // k-loop is over position in the current feature vector, y
      if (l != 0) {
	for (int k=0;k<n_ftrs;k++) {
	  if (k == j) {
	    if (k < (n_ftrs-1)) 
	      kms[KMSPOS(j,k+1,n_ftrs)].writeMgDoubleRecord2D(out_fp);
	    else 
	      kms[KMSPOS(j-1,k,n_ftrs)].writeMgDoubleRecord2D(out_fp);
	  } else if (k > j) {
	    kms[KMSPOS(j,k,n_ftrs)].writeMgDoubleRecord2D(out_fp);
	  } else { // k < j
	    kms[KMSPOS(k,j,n_ftrs)].writeMgDoubleRecord2D(out_fp);
	  }
	}
      } else {
	for (int k=(j+1);k<n_ftrs;k++) {
	  kms[KMSPOS(j,k,n_ftrs)].writeMgDoubleRecord2D(out_fp);
	}
      }
    }
  }

  delete [] ftr_buf;
  delete [] kms;
}
