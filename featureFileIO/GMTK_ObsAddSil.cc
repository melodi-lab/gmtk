/*
    $Header$
  
    This program selects an arbitrary subset set of features from each
    frame of a pfile and creates a new pfile with that
    subset in each frame but with "silence" added to the
    beginning and ending of each utterance. Silence is defined
    by computing means and variances from select ranges within each
    utterance (one range each for producing silence at the beginning
    and the end), and then randomly sampling from a Normal distribution
    with the correspondingly computed means and variances.

    Written By:
         Jeff Bilmes <bilmes@cs.berkeley.edu>
    Modified by:
         Katrin Kirchhoff <katrin@ee.washington.edu>
   Modified by:
         Karim Filali <karim@cs.washington.edu>  to support ascii, binary and htk input and output formats.
*/


#include "GMTK_ObsAddSil.h"
#include "rand.h"
// #include <values.h>
#include <math.h>


RAND rnd(true);


void addSil(FILE* out_fp, 
	     ObservationMatrix* obs_mat,
	     //InFtrLabStream_PFile in_stream,
	     //OutFtrLabStream_PFile out_stream,
	     Range& srrng,
	     const int nb,
	     const char *prb_str,
	     const int ne,
	     const char *pre_str,
	     const double mmf,
	     const double maf,
	     const double smf,
	     const double saf,
	     const bool dontPrintFrameID,
	     const bool quiet,
	     unsigned ofmt,
	     int debug_level,
	     bool oswap)
{
    // Feature and label buffers are dynamically grown as needed.
    size_t buf_size = 300;      // Start with storage for 300 frames.
    const size_t n_labs = obs_mat->numDiscrete();
    const size_t n_ftrs = obs_mat->numContinuous();
    
    OutFtrLabStream_PFile* out_stream=NULL;
    if(ofmt==PFILE) {
      out_stream = new OutFtrLabStream_PFile(debug_level,"",out_fp,n_ftrs,n_labs,1,oswap);
    }
    

    float *ftr_buf = new float[buf_size * n_ftrs];
    float *ftr_buf_p;
    float *oftr_buf = new float[(ne+nb+buf_size) * n_ftrs];
    float *oftr_buf_p;

    UInt32* lab_buf = new UInt32[buf_size * n_labs];
    UInt32* lab_buf_p;    
    UInt32* olab_buf = new UInt32[(ne+nb+buf_size) * n_labs];
    UInt32* olab_buf_p;    

    double *sums = new double[n_ftrs];
    double *sumsqs = new double[n_ftrs];

    //
    // Go through input pfile to get the initial statistics,
    // i.e., max, min, mean, std, etc.
    for (Range::iterator srit=srrng.begin();!srit.at_end();srit++) {
      obs_mat->loadSegment((const unsigned)(*srit));
      const size_t n_frames = obs_mat->numFrames();
      
      //      const size_t n_frames = in_stream.num_frames((*srit));

	Range prbrng(prb_str,0,n_frames);
	Range prerng(pre_str,0,n_frames);

	if (prbrng.length() == 0 && nb != 0)
	  error("No frames to compute beginning silence, sentence %d\n",
		(*srit));
	if (prerng.length() == 0 && ne != 0)
	  error("No frames to compute beginning silence, sentence %d\n",
		(*srit));

	if ((*srit) % 100 == 0)
	  printf("Processing sentence %d\n",(*srit));

        // Increase size of buffers if needed.
        if (n_frames > buf_size)
        {
            // Free old buffers.
            delete [] ftr_buf;
            delete [] oftr_buf;
            delete [] lab_buf;
            delete [] olab_buf;

            // Make twice as big to cut down on future reallocs.
            buf_size = n_frames * 2;

            // Allocate new larger buffers.
            ftr_buf = new float[buf_size * n_ftrs];
            oftr_buf = new float[(nb+ne+buf_size) * n_ftrs];
            lab_buf = new UInt32[buf_size * n_labs];
            olab_buf = new UInt32[(nb+ne+buf_size) * n_labs];
        }

	for(unsigned frame_no = 0;  frame_no < n_frames; ++frame_no) {
	  const float* start_of_frame = obs_mat->floatVecAtFrame(frame_no);
	  const UInt32* start_of_unsigned_frame = obs_mat->unsignedAtFrame(frame_no);
	  for(unsigned feat_no = 0;  feat_no < n_ftrs; ++feat_no) {
	    ftr_buf[frame_no * n_ftrs + feat_no] = *(start_of_frame  + feat_no);
	  }
	  for(unsigned unsigned_feat_no = 0;  unsigned_feat_no < n_labs; ++unsigned_feat_no) {
	    lab_buf[frame_no*n_labs + unsigned_feat_no] = *(start_of_unsigned_frame+unsigned_feat_no);
	  }
	}
	

	// compute the beginning means and vars
	oftr_buf_p = oftr_buf;
	olab_buf_p = olab_buf;
	if (nb > 0) {
	  ::memset(sums,0,sizeof(double)*n_ftrs);
	  ::memset(sumsqs,0,sizeof(double)*n_ftrs);	
	  for (Range::iterator prit=prbrng.begin();
	       !prit.at_end(); ++prit) {
	    ftr_buf_p = ftr_buf + (*prit)*n_ftrs;
	    int i=0;
	    //	    for (Range::iterator frit=frrng.begin();!frit.at_end(); ++frit) {
	    for (unsigned frit=0;frit<n_ftrs; ++frit) {
	      double tmp = ftr_buf_p[frit];
	      sums[i] += tmp;
	      sumsqs[i] += tmp*tmp;
	      i++;
	    }
	  }
	  for (unsigned i=0;i<n_ftrs;i++) {
	    sums[i] /= prbrng.length();
	    sumsqs[i] = sumsqs[i] / prbrng.length() - sums[i]*sums[i];
	    if (sumsqs[i] <= DBL_MIN)
	      sumsqs[i] = 0.0;
	    else 
	      sumsqs[i] = sqrt(sumsqs[i]);
	    sums[i] = mmf*sums[i] + maf;
	    sumsqs[i] = smf*sumsqs[i] + saf;
	  }
	  for (int i=0;i<nb;i++) {
	    //	    for (Range::iterator lrit=lrrng.begin(); !lrit.at_end(); ++lrit) {
	    for (unsigned lrit=0; lrit<n_labs; ++lrit) {
	      // copy labels from first frame only.
	      *olab_buf_p++ = lab_buf[lrit];
	    }
	    for (unsigned j=0;j<n_ftrs;j++) {
	      *oftr_buf_p++ = sums[j] + 
		sumsqs[j]*rnd.inverse_normal_func(rnd.drand48pe());
	    }
	  }
	}

	// copy normal frames.
	for (unsigned frame=0;frame<n_frames;frame++) {
	  ftr_buf_p = ftr_buf + (frame)*n_ftrs;
	  lab_buf_p = lab_buf + (frame)*n_labs;
	  for (unsigned frit=0;frit<n_ftrs; ++frit) {
	    *oftr_buf_p++ = ftr_buf_p[frit];
	  }
	  for (unsigned lrit=0; lrit<n_labs; ++lrit) {
	    *olab_buf_p++ = lab_buf_p[lrit];
	  }
	}

	if (ne > 0) {
	  ::memset(sums,0,sizeof(double)*n_ftrs);
	  ::memset(sumsqs,0,sizeof(double)*n_ftrs);	
	  for (Range::iterator prit=prerng.begin();
	       !prit.at_end(); ++prit) {
	    ftr_buf_p = ftr_buf + (*prit)*n_ftrs;
	    int i=0;
	    for (unsigned frit=0;frit<n_ftrs; ++frit) { 
	      double tmp = ftr_buf_p[frit];
	      sums[i] += tmp;
	      sumsqs[i] += tmp*tmp;
	      i++;
	    }
	  }
	  for (unsigned i=0;i<n_ftrs;i++) {
	    sums[i] /= prerng.length();
	    sumsqs[i] = sumsqs[i] / prerng.length() - sums[i]*sums[i];
	    if (sumsqs[i] <= DBL_MIN)
	      sumsqs[i] = 0.0;
	    else 
	      sumsqs[i] = sqrt(sumsqs[i]);
	    sums[i] = mmf*sums[i] + maf;
	    sumsqs[i] = smf*sumsqs[i] + saf;
	  }
	  for (int i=0;i<nb;i++) {
	    for (unsigned lrit=0; lrit<n_labs; ++lrit) {
	      // copy label from last frame only.
	      *olab_buf_p++ = (lab_buf + (n_frames-1)*n_labs)[lrit];
	    }
	    for (unsigned j=0;j<n_ftrs;j++) {
	      *oftr_buf_p++ = sums[j] + 
		sumsqs[j]*rnd.inverse_normal_func(rnd.drand48pe());
	    }
	  }
	}

	// Write output.
	 printSegment(*srit, out_fp, oftr_buf,n_ftrs,lab_buf,n_labs,n_frames+nb+ne, dontPrintFrameID,quiet, ofmt, debug_level, oswap, out_stream);

    }


    delete [] ftr_buf;
    delete [] lab_buf;
    delete [] sums;
    delete [] sumsqs;
}

