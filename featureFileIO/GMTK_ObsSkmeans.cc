/*
    $Header$
  
    This program computes utterance specific segmental k-means.
     
    That is, for each segment of the same utterance (e.g., word),
    N segments of equal length are defined and the k-means 
    algorithm is computed on each of them.
    The resulting means and variances are printed as output.

    This program was hacked together very quickly so there is currently
    little C-code optimization performed.
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

#include "vbyteswapping.h"


#include "GMTK_Kmeans.h"
#include "GMTK_ObsInitmg.h"


ObservationMatrix globalObservationMatrix;

void (*copy_swap_func_ptr)(size_t, const intv_int32_t*, intv_int32_t*)=NULL;

#define SENTS_PER_PRINT 1


/////////////////////////////////////////////////////////////////////////////
// SentIdStream_File returns a stream of sentence id strings from a file.
// TOOD: Make separate abstract base class.
/////////////////////////////////////////////////////////////////////////////

class SentIdStream_File
{
public:
  SentIdStream_File(FILE *fp,const int nw);
  ~SentIdStream_File();

  void rewind() { if (fp != NULL) ::rewind(fp); }
  const char *next();         // Returns pointer to next sentence id, or
                              // 0 if at end of file or error.

  operator char*() { return buf; }
  operator size_t() { return idx; }

  size_t number_read() { return nread; }
private:
  FILE *fp;
  size_t buf_size;
  char *buf;
  size_t nread;
  size_t idx;
  int num_words;
};

SentIdStream_File::SentIdStream_File(FILE *fp_arg,const int nw)
    : fp(fp_arg),num_words(nw)
{
    buf_size = 128;
    buf = new char[buf_size];
    idx = 0;
    ::sprintf(buf,"0");
}

SentIdStream_File::~SentIdStream_File()
{
    delete buf;
}

const char*
SentIdStream_File::next()
{

  // TODO: Doesn't distinguish file I/O errors from end of file condition.
  if (fp == NULL) {
    // if no file, return infinite stream of 0's
    return buf;
  }

  if (fgets(buf, buf_size, fp)==NULL)
    {
      error("Error: EOF in labels file, or I/O error. nread=%ld",
	    nread);
      return 0;
    }
  else
    {
      char *p = strchr(buf,'\n'); // Get pos. of newline character.
      if (p==NULL)
        {
	  error("Sentence id too long.");
        }
      else
        {
	  *p = '\0';          // Trim off newline character.
        }   
      char *ptr;
      int tmp = ::strtol(buf,&ptr,10);
      if (ptr == buf) 
	error("No integer at position %d in file.",nread);
      if (tmp >= num_words) 
	error("Word id (%d) > num_words (%d).",idx,num_words);
      idx = (size_t)tmp;
      nread++;
      return buf;
    }
}


/////////////////////////////////////////////////////////////////////////////
// pfile_skmeans: compute segmental kmeans
/////////////////////////////////////////////////////////////////////////////


static void uniformSkmeans(ObservationMatrix* obs_mat,
			   FILE *out_fp,
			   SentIdStream_File& sid_stream,
			   const int num_words, const int num_segments, 
			   const int num_clusters,const int maxReInits,
			   const int minSamplesPerCluster,
			   const int numRandomReStarts,
			   const char*pr_str,
			   const bool binary,const bool quiet_mode,
			   const float conv_thres,
			   const bool prefetch)
{

    // Feature and label buffers are dynamically grown as needed.

  //    size_t buf_size = 300;      // Start with storage for 300 frames.
    //    const size_t n_labs = in_streamp->n_labs();
  //const size_t n_labs = obs_mat->numDiscrete();
    //const size_t n_ftrs = in_streamp->n_ftrs();
    const size_t n_ftrs = obs_mat->numContinuous();
    const size_t num_stream_sentences = obs_mat->numSegments();

    float *ftr_buf = NULL;
    UInt32* lab_buf = NULL;

    kmeans::kmeans_k = num_clusters;
    kmeans::kmeans_vl = n_ftrs;

    kmeans *kms = new kmeans[num_words*num_segments];
    size_t sent_no;

    double bestVarianceSum=DBL_MAX;
    
    for (int epoch=0;epoch<numRandomReStarts;epoch++) {

      printf("Epoch %d\n",epoch);

      int i,j,iter=0;
      int reInits = 0;
    
      float max_dist = 0;

      for (i=0;i<num_words*num_segments;i++) {
	kms[i].randomAssignment = true;
	kms[i].done = false;
      }
      do {
	iter++;

	sid_stream.rewind();
	for (i=0;i<num_words*num_segments;i++) {
	  if (!kms[i].done)
	    kms[i].initNew();
	}

	pid_t pid = (pid_t)1;
	//	for (sent_no=0;sent_no<in_streamp->n_segs();sent_no++) {
	for (sent_no=0;sent_no<num_stream_sentences;sent_no++) {

	  if (prefetch && pid != 0 && sent_no > 0) {
	    // pid != 0 means this is the parent
	    // sent_no > 0 means there must have been a child since
	    // this is not the first iteration.
	    wait(0); // wait for the child to complete
	  }
	  
	  obs_mat->loadSegment(sent_no);
	  //	  const size_t n_frames = in_streamp->n_frms((sent_no));
	  const size_t n_frames = obs_mat->numFrames();

	  sid_stream.next();
	  size_t word_id = (size_t) sid_stream;
	
	  if (word_id >= (size_t)num_words) {
	    error("Invalid word_id (%d) at position %d\n",
		  word_id,sent_no+1);
	  }

	  if (!quiet_mode) {
	    if ((sent_no) % SENTS_PER_PRINT == 0)
	      printf("Processing sentence %d\n",(sent_no));
	  }

	  //	  const size_t n_read = in_streamp->read_ftrslabs(sent_no,ftr_buf, lab_buf);
	  for(unsigned frame_no = 0;  frame_no < n_frames; ++frame_no) {
	    const float* start_of_frame = obs_mat->floatVecAtFrame(frame_no);
	    for(unsigned feat_no = 0;  feat_no < n_ftrs; ++feat_no) {
	      ftr_buf[frame_no * n_ftrs + feat_no] = *(start_of_frame  + feat_no);
	    }
	  }

	  if (prefetch) {
	    // this is a hack to pre-fetch the next sentence
	    // and optimize Solaris's disk cache strategy.
	    if (pid != (pid_t)0) {
	      // then this is the parent process
	      //	      if ((sent_no+1)<in_streamp->n_segs()) {
	      if ((sent_no+1)<num_stream_sentences) {
		// then there is a next iteration.
		pid = fork();
	      }
	    } else {
	      // then this is the child process
	      // who just read in the next iteration's data into disk cache.
	      exit(0);
	    }
	  
	    if (pid == (pid_t)0) {
	      // then this is the new child,
	      // continue with the next iteration.
	      continue;
	    }
	  }
	  
	  Range prrng(pr_str,0,n_frames);
	  const int frames_per_segment = prrng.length()/num_segments;
	  int cur_seg_no = 0;
	  int segs_so_far = 0;
	  for (Range::iterator prit=prrng.begin();
	       !prit.at_end() ; ++prit) {

	    if (!kms[word_id*num_segments + cur_seg_no].done) {
	      // compute a pointer to the current buffer.
	      const float *const ftr_buf_p = ftr_buf + (*prit)*n_ftrs;
	      if (kms[word_id*num_segments + cur_seg_no].randomAssignment)
		kms[word_id*num_segments + cur_seg_no].add2newRand(ftr_buf_p);
	      else
		kms[word_id*num_segments + cur_seg_no].add2new(ftr_buf_p);
	    }

	    segs_so_far++;
	    if (segs_so_far >= frames_per_segment 
		&& (cur_seg_no+1 < num_segments)) {
	      segs_so_far = 0;
	      cur_seg_no++;
	    }
	  }
	}

	for (i=0;i<num_words*num_segments;i++) {
	  if (!kms[i].done)
	    kms[i].randomAssignment = false;
	}

	// make sure each kmeans had some data
	max_dist = 0;
	int num_active=0;
	int num_reinits = 0;
	for (i=0;i<num_words*num_segments;i++) {
	  if (!kms[i].done) {
	    if (kms[i].zeroCounts()) {
	      // This shouldn't happen with uniform skmeans. All kmeans
	      // objects should be getting some data.
	      error("kms[%d] was given no input. Probably an command line argument error");
	    } else if (kms[i].someClusterHasLessThanNEntries(minSamplesPerCluster)) {
	      // fprintf(stderr,
	      //       "Warning: kms word %d seg %d, some clusters have < %d entries. Re-initializing.\n",i/num_segments,i % num_segments,minSamplesPerCluster);
	    
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
		// printf("Iter %d: kms word %d seg %d converged\n",iter,
		//       i/num_segments, i % num_segments);
	      }
	    }
	    num_active++;
	  }
	}
	if (num_reinits)
	  reInits++;
	printf("Iter %d: max_dist = %e, cur num_reinits = %d, tot itr re_init = %d, num_active = %d\n",iter,max_dist,num_reinits,reInits,num_active);
	fflush(stdout);

      } while (max_dist > conv_thres && reInits <= maxReInits);

      if (reInits > maxReInits) {
	error("Error. %d re-inits and convergence didn't occur.", reInits);
      }

      // Do one more pass over file to compute the variances.
      sid_stream.rewind();
      //      for (sent_no=0;sent_no<in_streamp->n_segs();sent_no++) {
      for (sent_no=0;sent_no<num_stream_sentences;sent_no++) {
	obs_mat->loadSegment(sent_no);
	//	  const size_t n_frames = in_streamp->n_frms((sent_no));
	const size_t n_frames = obs_mat->numFrames();
	//const size_t n_frames = in_streamp->n_frms((sent_no));
	sid_stream.next();
	size_t word_id = (size_t) sid_stream;
	
	if (word_id >= (size_t)num_words) {
	  error("Invalid word_id (%d) at position %d\n",
		word_id,sent_no+1);
	}

	if (!quiet_mode) {
	  if ((sent_no) % SENTS_PER_PRINT == 0)
	    printf("Processing sentence %d\n",(sent_no));
	}

	//	const size_t n_read = in_streamp->read_ftrslabs(sent_no, ftr_buf, lab_buf);
	for(unsigned frame_no = 0;  frame_no < n_frames; ++frame_no) {
	  const float* start_of_frame = obs_mat->floatVecAtFrame(frame_no);
	  for(unsigned feat_no = 0;  feat_no < n_ftrs; ++feat_no) {
	    ftr_buf[frame_no * n_ftrs + feat_no] = *(start_of_frame  + feat_no);
	  }
	}


	Range prrng(pr_str,0,n_frames);
	const int frames_per_segment = prrng.length()/num_segments;
	int cur_seg_no = 0;
	int segs_so_far = 0;
	for (Range::iterator prit=prrng.begin();
	     !prit.at_end() ; ++prit) {

	  // compute a pointer to the current buffer.
	  const float *const ftr_buf_p = ftr_buf + (*prit)*n_ftrs;
	  kms[word_id*num_segments + cur_seg_no].computeVariances(ftr_buf_p);

	  segs_so_far++;
	  if (segs_so_far >= frames_per_segment 
	      && (cur_seg_no+1 < num_segments)) {
	    segs_so_far = 0;
	    cur_seg_no++;
	  }
	}
      }

    
      kmeans *kmsp = kms;
      double sumVariances=0;
      for (i=0;i<num_words;i++) {
	for (j=0;j<num_segments;j++) {
	  sumVariances+=kmsp->finishVariances();
	  kmsp++;
	}
      }
      
      printf("End of Epoch %d of %d, variance sum = %e, best = %e\n",epoch,
	     numRandomReStarts,
	     sumVariances,bestVarianceSum);

      if (sumVariances < bestVarianceSum) {
	bestVarianceSum = sumVariances;
	kmsp = kms;
	for (i=0;i<num_words;i++) {
	  for (j=0;j<num_segments;j++) {
	    kmsp->save();
	    kmsp++;
	  }
	}
      }
    }
    
    kmeans *kmsp = kms;
    for (int i=0;i<num_words;i++) {
      for (int j=0;j<num_segments;j++) {
	fprintf(out_fp,"word %d, seg %d:\n",i,j);
	kmsp->printSaved(out_fp);
	kmsp++;
      }
    }

    delete [] ftr_buf;
    delete [] lab_buf;
    delete [] kms;
}




static void
//pfile_viterbi_skmeans(SPI_base *in_streamp,
viterbiSkmeans(ObservationMatrix * obs_mat,
	       FILE *out_fp,
	       //		      SPI_base* in_lstreamp,
	       const int num_labels,
	       const int num_clusters,const int maxReInits,
	       const int minSamplesPerCluster,
	       const int numRandomReStarts,
	       const char*pr_str,
	       const bool binary,const bool quiet_mode,
	       const float conv_thres,
	       const bool prefetch)
{

    // Feature and label buffers are dynamically grown as needed.

  //size_t buf_size = 300;      // Start with storage for 300 frames.
    //    const size_t n_ftrs = in_streamp->n_ftrs();
    //const size_t n_labs = (in_lstreamp==NULL? in_streamp->n_labs():in_lstreamp->n_labs());
     const size_t n_labs = obs_mat->numDiscrete();
     const size_t n_ftrs = obs_mat->numContinuous();
    
     const size_t num_stream_sentences = obs_mat->numSegments();

    if (n_labs != 1)
      error("File containing labels must have exactly 1 label per frame.  Use -lr option.");

    if (n_ftrs == 0)
      error("Observation file must have more than 0 features.");

    float *ftr_buf = NULL;
    UInt32* lab_buf = NULL;

    kmeans::kmeans_k = num_clusters;
    kmeans::kmeans_vl = n_ftrs;

    kmeans *kms = new kmeans[num_labels];
    size_t sent_no;

    double bestVarianceSum=DBL_MAX;

    //    if (in_lstreamp != NULL) {
    //      if (in_lstreamp->n_segs() != in_streamp->n_segs())
    //	error("Feature and label pfile have differing number of sentences.");
    //}
    
    for (int epoch=0;epoch<numRandomReStarts;epoch++) {

      printf("Epoch %d\n",epoch);

      int i,iter=0;
      int reInits = 0;
    
      float max_dist = 0;

      for (i=0;i<num_labels;i++) {
	kms[i].randomAssignment = true;
	kms[i].done = false;
      }

      do {
	iter++;

	for (i=0;i<num_labels;i++) {
	  if (!kms[i].done)
	    kms[i].initNew();
	}

	pid_t pid = (pid_t)1;
	//	for (sent_no=0;sent_no<in_streamp->n_segs();sent_no++) {
	for (sent_no=0;sent_no<num_stream_sentences;sent_no++) {

	  if (prefetch && pid != 0 && sent_no > 0) {
	    // pid != 0 means this is the parent
	    // sent_no > 0 means there must have been a child since
	    // this is not the first iteration.
	    wait(0); // wait for the child to complete
	  }

	  obs_mat->loadSegment(sent_no);
	  //	  const size_t n_frames = in_streamp->n_frms((sent_no));
	  const size_t n_frames = obs_mat->numFrames();
	  //if (in_lstreamp != NULL) {
	  //if (in_lstreamp->n_frms(sent_no) != n_frames)
	  //  error("Feature and label pfile have differing number of frames at sentence %d.",sent_no);
	  //}

	  if (!quiet_mode) {
	    if ((sent_no) % 1 == 0)
	      printf("Processing sentence %d\n",(sent_no));
	  }

	  //	  const size_t n_read = in_streamp->read_ftrslabs(sent_no,ftr_buf, lab_buf);

	  for(unsigned frame_no = 0;  frame_no < n_frames; ++frame_no) {
	    const float* start_of_frame = obs_mat->floatVecAtFrame(frame_no);
	    const UInt32* start_of_unsigned_frame = obs_mat->unsignedAtFrame(frame_no);
	    for(unsigned feat_no = 0;  feat_no < n_ftrs; ++feat_no) {
	      ftr_buf[frame_no * n_ftrs + feat_no] = *(start_of_frame  + feat_no);
	    }
	    //	    for(unsigned unsigned_feat_no = 0;  unsigned_feat_no < n_labs; ++unsigned_feat_no) {
	    //*lab_buf++ = *(start_of_unsigned_frame+unsigned_feat_no);
	    *lab_buf++ = *(start_of_unsigned_frame);  // will only take first label
	      //}
	  }
	  
	  //if (in_lstreamp != NULL) {
	  //const size_t n_read = in_lstreamp->read_labs(sent_no, lab_buf);
	  //}

	  if (prefetch) {
	    // this is a hack to pre-fetch the next sentence from disk
	    // and optimize any disk cache strategy.
	    if (pid != (pid_t)0) {
	      // then this is the parent process
	      if ((sent_no+1)<num_stream_sentences) {
		// then there is a next iteration.
		pid = fork();
	      }
	    } else {
	      // then this is the child process
	      // who just read in the next iteration's data into disk cache.
	      exit(0);
	    }
	  
	    if (pid == (pid_t)0) {
	      // then this is the new child,
	      // continue with the next iteration.
	      continue;
	    }
	  }

	  Range prrng(pr_str,0,n_frames);
	  for (Range::iterator prit=prrng.begin();
	       !prit.at_end() ; ++prit) {

	    const size_t curLab = lab_buf[(*prit)];
	    if ((int)curLab >= num_labels)
	      error("Label at sentence %d, frame %d is %d and is >= %d.",
		    sent_no,(*prit),curLab,num_labels);

	    if (!kms[curLab].done) {
	      // compute a pointer to the current buffer.
	      const float *const ftr_buf_p = ftr_buf + (*prit)*n_ftrs;
	      if (kms[curLab].randomAssignment)
		kms[curLab].add2newRand(ftr_buf_p);
	      else
		kms[curLab].add2new(ftr_buf_p);
	    }

	  }
	}

	for (i=0;i<num_labels;i++) {
	  if (!kms[i].done)
	    kms[i].randomAssignment = false;
	}


	// make sure each kmeans had some data
	max_dist = 0;
	int num_active=0;
	int num_reinits = 0;
	for (i=0;i<num_labels;i++) {
	  if (!kms[i].done) {
	    if (kms[i].zeroCounts()) {
	      // this is ok, this label perhaps doesn't exist in the file.
	    } else if (kms[i].someClusterHasLessThanNEntries(minSamplesPerCluster)) {
	      // fprintf(stderr,
	      //       "Warning: kms label %d, some clusters have < %d entries. Re-initializing.\n",i,minSamplesPerCluster);
	      
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
		printf("Iter %d: label %d converged\n",iter,i);
	      }
	    }
	  }
	  num_active++;
	}
	if (num_reinits)
	  reInits++;
	printf("Iter %d: max_dist = %e, cur num_reinits = %d, tot itr re_init = %d, num_active = %d\n",iter,max_dist,num_reinits,reInits,num_active);
	fflush(stdout);

      } while (max_dist > conv_thres && reInits <= maxReInits);

      if (reInits > maxReInits) {
	error("Error. %d re-inits and convergence didn't occur.", reInits);
      }

      // Do one more pass over file to compute the variances.
      //      for (sent_no=0;sent_no<in_streamp->n_segs();sent_no++) {
      for (sent_no=0;sent_no<num_stream_sentences;sent_no++) {
	obs_mat->loadSegment(sent_no);
	//	  const size_t n_frames = in_streamp->n_frms((sent_no));
	const size_t n_frames = obs_mat->numFrames();

	//	if (in_lstreamp != NULL) {
	//if (in_lstreamp->n_frms(sent_no) != n_frames)
	//  error("Feature and label pfile have differing number of frames at sentence %d.",sent_no);
	//}

	if (!quiet_mode) {
	  if ((sent_no) % SENTS_PER_PRINT == 0)
	    printf("Processing sentence %d\n",(sent_no));
	}

	//	const size_t n_read = in_streamp->read_ftrslabs(sent_no, ftr_buf, lab_buf);

	//	if (in_lstreamp != NULL) {
	//	  const size_t n_read = in_lstreamp->read_labs(sent_no, lab_buf);
	//}

	for(unsigned frame_no = 0;  frame_no < n_frames; ++frame_no) {
	  const float* start_of_frame = obs_mat->floatVecAtFrame(frame_no);
	  for(unsigned feat_no = 0;  feat_no < n_ftrs; ++feat_no) {
	    ftr_buf[frame_no * n_ftrs + feat_no] = *(start_of_frame  + feat_no);
	  }
	}

	Range prrng(pr_str,0,n_frames);
	for (Range::iterator prit=prrng.begin();
	     !prit.at_end() ; ++prit) {

	  const size_t curLab = lab_buf[(*prit)];
	  // compute a pointer to the current buffer.
	  const float *const ftr_buf_p = ftr_buf + (*prit)*n_ftrs;
	  kms[curLab].computeVariances(ftr_buf_p);
	}
      }

    
      kmeans *kmsp = kms;
      double sumVariances=0;
      for (i=0;i<num_labels;i++) {
	if (!kmsp->zeroCounts())
	  sumVariances+=kmsp->finishVariances();
	kmsp++;
      }
      

      if (sumVariances < bestVarianceSum) {
	bestVarianceSum = sumVariances;
	kmsp = kms;
	for (i=0;i<num_labels;i++) {
	  if (!kmsp->zeroCounts())
	    kmsp->save();
	  kmsp++;
	}
      }
      printf("End of Epoch %d of %d, variance sum = %e, best = %e\n",epoch,
	     numRandomReStarts,
	     sumVariances,bestVarianceSum);


    }

    kmeans *kmsp = kms;
    for (int i=0;i<num_labels;i++) {
      if (!kmsp->zeroCounts()) {
	fprintf(out_fp,"label %d:\n",i);
	kmsp->printSaved(out_fp);
      }
      kmsp++;
    }

    delete [] ftr_buf;
    delete [] lab_buf;
    delete [] kms;
}



#define MAX_OBJECTS 5

char *   input_fname[MAX_OBJECTS] = {NULL,NULL,NULL,NULL,NULL};  // Input file name.
const char *   ifmtStr[MAX_OBJECTS]     = {"pfile","pfile","pfile","pfile","pfile"};
unsigned ifmt[MAX_OBJECTS];

char *   output_fname      = NULL;

char *   input_uname       = NULL;
bool     Print_Stream_Info = false;
bool     Print_Sent_Frames = false;

unsigned nis[MAX_OBJECTS];
unsigned nfs[MAX_OBJECTS];


char  *sr_str               = 0;   // sentence range string
Range *sr_rng;
char  *fr_str[MAX_OBJECTS]  = {NULL,NULL,NULL,NULL,NULL};   // feature range string    
char  *lr_str[MAX_OBJECTS]  = {NULL,NULL,NULL,NULL,NULL};   // label range string  
char  *spr_str[MAX_OBJECTS] = {NULL,NULL,NULL,NULL,NULL};   // per stream per sentence range string 
char  *pr_str               = 0;   // per-sentence range string


const char*    actionIfDiffNumFramesStr[MAX_OBJECTS]={"er","er","er","er","er"};   // 
unsigned actionIfDiffNumFrames[MAX_OBJECTS]={FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR,FRAMEMATCH_ERROR};   // 
const char*    actionIfDiffNumSentsStr[MAX_OBJECTS]={"te","te","te","te","te"}; 
unsigned actionIfDiffNumSents[MAX_OBJECTS]={SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END,SEGMATCH_TRUNCATE_FROM_END};   // 

bool     quiet = false;

#ifdef INTV_WORDS_BIGENDIAN
bool iswap[MAX_OBJECTS]={true,true,true,true,true};
bool oswap = true;
#else
bool iswap[MAX_OBJECTS]= {false,false,false,false,false};
bool oswap             = false;
#endif

char* perStreamTransforms[MAX_OBJECTS]={NULL,NULL,NULL,NULL,NULL};   // 
char* postTransforms=NULL;

unsigned startSkip = 0;
unsigned endSkip   = 0;


const char*    ftrcomboStr = "none";
unsigned ftrcombo    = FTROP_NONE;

bool     cppIfAscii        = true;
char*    cppCommandOptions = NULL;

bool     help              = false;


int num_words = 0;
int num_segments = 0;
int num_clusters = 0;
int maxReInits = 50;
int numRandomReStarts = 1;
int minSamplesPerCluster = 1;
float conv_thres = 0.0;

bool binary = false;
int  debug_level = 0;

// if not true, we do viterbi k-means
bool uniform_skm = true;
bool viterbi     = false;  
bool prefetch    = false;

bool  Init_MG     = false;
char* Init_MG_CFR = NULL;

Arg Arg::Args[] = {
  Arg("i",        Arg::Req, input_fname,  "Input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("ifmt",     Arg::Opt, ifmtStr ,     "Format of input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("o",        Arg::Opt, output_fname, "Output file"),
  Arg("b",        Arg::Opt, binary,       "Binary rather than ascii output"),
  Arg("nf",       Arg::Opt, nfs,          "Number of floats in input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("ni",       Arg::Opt, nis,          "Number of ints (labels) in input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("sr",       Arg::Opt, sr_str,       "Sentence range"),
  Arg("fr",       Arg::Opt, fr_str,       "Feature range",Arg::ARRAY,MAX_OBJECTS),
  Arg("spr",      Arg::Opt, spr_str,      "Per stream per-sentence range",Arg::ARRAY,MAX_OBJECTS),
  Arg("lr",       Arg::Opt, lr_str,       "Label range",Arg::ARRAY,MAX_OBJECTS),
  Arg("startskip",Arg::Opt, startSkip,    "Start skip"),
  Arg("endskip",  Arg::Opt, endSkip,      "End skip"),
  Arg("q",        Arg::Tog, quiet,        "Quiet."),
  Arg("fdiffact", Arg::Opt, actionIfDiffNumFramesStr , "Action if different number of frames in streams: error (er), repeat last frame (rl), first frame (rf), segmentally expand (se), truncate from start (ts), truncate from end (te)",Arg::ARRAY,MAX_OBJECTS),
  Arg("sdiffact", Arg::Opt, actionIfDiffNumSentsStr ,  "Action if different number of sentences in streams: error (er), truncate from end (te), repeat last sent (rl), and wrap around (wa).",Arg::ARRAY,MAX_OBJECTS),
  Arg("trans",    Arg::Opt, perStreamTransforms , "Per stream transformations string",Arg::ARRAY,MAX_OBJECTS),
  Arg("posttrans",Arg::Opt, postTransforms ,      "Final global transformations string"),
  Arg("comb",     Arg::Opt, ftrcomboStr,  "Combine float features (none: no combination, add, sub, mul,div"),
  Arg("k",        Arg::Opt, num_clusters, "Number of clusters (i.e., K)"),
  Arg("r",        Arg::Opt, numRandomReStarts, "Number of epoch random restarts to take best of"),
  Arg("m",        Arg::Opt, maxReInits,        "Maximum number of re-inits before giving up"),
  Arg("x",        Arg::Opt, minSamplesPerCluster, "Min number of samples/cluster for a re-init to occur"),
  Arg("a",        Arg::Opt, conv_thres,           "Convergence threshold"),
  Arg("u",        Arg::Opt, uniform_skm,          "Uniform segmental k-means"),
  Arg("v",        Arg::Opt, viterbi,              "Viterbi segmental k-means"),
  Arg("n",        Arg::Opt, num_words,    "Total number of different words if uniform skmeans (-u option) is used. Otherwise (viterbi skmeans, -v option, denotes the maximum number of labels (range [0:n-1])"),  
  Arg("f",        Arg::Opt, input_uname,  "Input utterance-ID file. (with -u option only)"),
  Arg("s",        Arg::Opt, num_segments, "Number of segments/words. (with -u option only)"),
  Arg("prefetch", Arg::Opt, prefetch,     "Prefetch next sentence at each iteration"),
  Arg("initmg",   Arg::Opt, Init_MG,      "Create an initialization .mg file for bivariate-mi.cc"), 
  Arg("initmgCfr",Arg::Opt, Init_MG_CFR,  "Range of past/future context frames to use when initializing an .mg file"), 
  Arg("iswap",    Arg::Opt, iswap,        "Do byte swapping on the input file",Arg::ARRAY,MAX_OBJECTS),
  Arg("cppifascii",        Arg::Opt, cppIfAscii,        "Pre-process ASCII files using CPP"),
  Arg("cppCommandOptions", Arg::Opt, cppCommandOptions, "Additional CPP command line"),
  Arg("help",     Arg::Tog, help,                       "Print this message"),
  // The argumentless argument marks the end of the above list.
  Arg()
};




int main(int argc, const char *argv[])
{

  int numFiles=0;

  // Figure out the Endian of the machine this is running on and set the swap defaults accordingly
  bool doWeSwap;

  ByteEndian byteEndian = getWordOrganization();
  switch(byteEndian) {
  case BYTE_BIG_ENDIAN:
    doWeSwap=false;
    break;
  case BYTE_LITTLE_ENDIAN:
    doWeSwap=true;
    break;
  default:
    // We weren't able to figure the Endian out.  Leave the swap defaults as they are.
#ifdef INTV_WORDS_BIGENDIAN
    doWeSwap=true;
#else
    doWeSwap=false;
#endif
  }

  oswap=doWeSwap;
  for(int i=0; i<MAX_OBJECTS; ++i) {
    iswap[i]=doWeSwap;
  }
  ///////////////////////////////////////////

  bool parse_was_ok = Arg::parse(argc,(char**)argv);

  if(help) {
    Arg::usage();
    exit(0);
  }

  if(!parse_was_ok) {
    Arg::usage(); exit(-1);
  }

  // oswap might been assigned a new value on the command line
  if(oswap) {
    copy_swap_func_ptr=&swapb_vi32_vi32;
  }
  else {
    copy_swap_func_ptr=&copy_vi32_vi32;
  }

    //////////////////////////////////////////////////////////////////////
    // Check all necessary arguments provided before creating objects.
    //////////////////////////////////////////////////////////////////////


    if (minSamplesPerCluster < 1) {
      error("Minimum samples per cluster (-x) must be >= 1.");
    }

    if(viterbi) uniform_skm = false;


    for (int i=0;i<MAX_OBJECTS;i++) {
      numFiles += (input_fname[i] != NULL);
      if (strcmp(ifmtStr[i],"htk") == 0)
	ifmt[i] = HTK;
      else if (strcmp(ifmtStr[i],"binary") == 0)
	ifmt[i] = RAWBIN;
      else if (strcmp(ifmtStr[i],"ascii") == 0)
      ifmt[i] = RAWASC;
      else if (strcmp(ifmtStr[i],"pfile") == 0)
	ifmt[i] = PFILE;
      //    else if (strcmp(ifmtStr[i],"flatbin") == 0)
      //  ifmt[i]=FLATBIN;
      //else if (strcmp(ifmtStr[i],"flatasc") == 0)
      //  ifmt[i]=FLATASC;
      else
	error("ERROR: Unknown observation file format type: '%s'\n",ifmtStr[i]);
    }
    
    for(int i=0; i < numFiles; ++i)
      if(output_fname!=NULL && strcmp(input_fname[i],output_fname)==0) {
     error("Input and output filenames cannot be the same.");
      }
    
    if (strcmp(ftrcomboStr,"none") == 0)
      ftrcombo = FTROP_NONE;
    else if (strcmp(ftrcomboStr,"add") == 0)
      ftrcombo = FTROP_ADD;
    else if (strcmp(ftrcomboStr,"sub") == 0)
      ftrcombo = FTROP_SUB;
    else if (strcmp(ftrcomboStr,"mul") == 0)
      ftrcombo = FTROP_MUL;
    else if (strcmp(ftrcomboStr,"div") == 0)
      ftrcombo = FTROP_DIV;
    else
      error("ERROR: Unknown feature combination type: '%s'\n",ftrcomboStr);
    
    for(int i=0; i < MAX_OBJECTS; ++i) {
      if(input_fname[i]!=NULL) {
	if (strcmp(actionIfDiffNumFramesStr[i],"er") == 0)
	  actionIfDiffNumFrames[i] = FRAMEMATCH_ERROR;
	else if (strcmp(actionIfDiffNumFramesStr[i],"rl") == 0)
	  actionIfDiffNumFrames[i] = FRAMEMATCH_REPEAT_LAST;
	else if (strcmp(actionIfDiffNumFramesStr[i],"rf") == 0)
	  actionIfDiffNumFrames[i] = FRAMEMATCH_REPEAT_FIRST;
	else if (strcmp(actionIfDiffNumFramesStr[i],"se") == 0)
	  actionIfDiffNumFrames[i] = FRAMEMATCH_EXPAND_SEGMENTALLY;
	else if (strcmp(actionIfDiffNumFramesStr[i],"ts") == 0)
	  actionIfDiffNumFrames[i] = FRAMEMATCH_TRUNCATE_FROM_START;
	else if (strcmp(actionIfDiffNumFramesStr[i],"te") == 0)
	  actionIfDiffNumFrames[i] = FRAMEMATCH_TRUNCATE_FROM_END;
	else
	  error("ERROR: Unknown action when diff num of frames: '%s'\n",actionIfDiffNumFramesStr[i]);
      }
    }
    
    for(int i=0; i < MAX_OBJECTS; ++i) {
      if(input_fname[i]!=NULL) {
	if (strcmp(actionIfDiffNumSentsStr[i],"er") == 0)
	  actionIfDiffNumSents[i] = SEGMATCH_ERROR;
	else if (strcmp(actionIfDiffNumSentsStr[i],"rl") == 0)
	  actionIfDiffNumSents[i] = SEGMATCH_REPEAT_LAST;
	else if (strcmp(actionIfDiffNumSentsStr[i],"wa") == 0)
	  actionIfDiffNumSents[i] = SEGMATCH_WRAP_AROUND;
	else if (strcmp(actionIfDiffNumSentsStr[i],"te") == 0)
	  actionIfDiffNumSents[i] = SEGMATCH_TRUNCATE_FROM_END;
	else
	  error("ERROR: Unknown action when diff num of sentences: '%s'\n",actionIfDiffNumSentsStr[i]);
      }
    }
    
    
    FILE *out_fp=NULL;
    if (output_fname==NULL || !strcmp(output_fname,"-")) {
      out_fp = stdout;
    }
    else {
      if ((out_fp = fopen(output_fname, "w")) == NULL) {
	error("Couldn't open output file for writing.\n");
      }
    }
    

    //////////////////////////////////////////////////////////////////////
    // Create objects.
    //////////////////////////////////////////////////////////////////////

    // If we have a pfile, we can extract the number if features from the file directly
    for(int i=0; i < MAX_OBJECTS; ++i) {
      if(input_fname[i]!=NULL) {
	if(ifmt[i]==PFILE) {
	  FILE *in_fp = fopen(input_fname[i], "r");
	  if (in_fp==NULL) error("Couldn't open input pfile for reading.");
	  InFtrLabStream_PFile* in_streamp = new InFtrLabStream_PFile(debug_level,"",in_fp,1,iswap[i]);
	  nis[i]=in_streamp->num_labs();
	  nfs[i]=in_streamp->num_ftrs();
	  if (fclose(in_fp)) error("Couldn't close input pfile.");
	  delete in_streamp;
	}
	
	if(nis[i]==0 && nfs[i]==0) {
	  error("The number of floats and the number of ints cannot be both zero.");
	}
      }
    }
    
    globalObservationMatrix.openFiles(numFiles,  // number of files.   For now we use only one
				      (const char**)&input_fname,
				      (const char**)&fr_str,
				      (const char**)&lr_str,
				      (unsigned*)&nfs,
				      (unsigned*)&nis,
				      (unsigned*)&ifmt,
				      (bool*)&iswap,
				      startSkip,  // startSkip
				      endSkip,  // endSkip  pr_rng takes care of these two  
				      cppIfAscii,
				      cppCommandOptions,
				      (const char**)&spr_str,
				      actionIfDiffNumFrames,
				      actionIfDiffNumSents,
				      perStreamTransforms,
				      postTransforms,
				      ftrcombo);   
    
    
    sr_rng = new Range(sr_str,0,globalObservationMatrix.numSegments());

    //////////////////////////////////////////////////////////////////////
    // Do the work.
    //////////////////////////////////////////////////////////////////////
    if(Init_MG) {
      Range* Init_MG_CFR_Range = new Range(Init_MG_CFR,-(int)globalObservationMatrix.numFrames(),(int)globalObservationMatrix.numFrames());
      initmg(&globalObservationMatrix,out_fp,
	     *sr_rng, *Init_MG_CFR_Range, pr_str,
	     num_clusters,maxReInits,
	     minSamplesPerCluster,
	     numRandomReStarts,
	     quiet,conv_thres
	     );
      delete Init_MG_CFR;
    }
    else if (uniform_skm) {
      FILE *ut_fp=NULL;;
      if (input_uname != NULL) {
	ut_fp = fopen(input_uname, "r");
	if (ut_fp==NULL)
	  error("Couldn't open input utterance file for reading.");
      } else {
	// all utterances are assumed the same.
      }
      SentIdStream_File* sid_streamp
	= new SentIdStream_File(ut_fp,num_words);
      uniformSkmeans(&globalObservationMatrix,out_fp,*sid_streamp,
		     num_words,num_segments,num_clusters,maxReInits,
		     minSamplesPerCluster,
		     numRandomReStarts,
		     pr_str,
		     binary,quiet,conv_thres,
		     prefetch);
      if (ut_fp) {
	if (fclose(ut_fp))
	  error("Couldn't close utterance file.");
      }
    } 
    else {
      if(globalObservationMatrix.numDiscrete() < 1 ) {
	error("No label found.");
      }      
      viterbiSkmeans(&globalObservationMatrix,out_fp,
		     num_words,num_clusters,maxReInits,
		     minSamplesPerCluster,
		     numRandomReStarts,
		     pr_str,
		     binary,quiet,conv_thres,
		     prefetch);
    }
    
    //////////////////////////////////////////////////////////////////////
    // Clean up and exit.
    //////////////////////////////////////////////////////////////////////

    //delete in_streamp;
    if (fclose(out_fp))
        error("Couldn't close output file.");

    delete sr_rng;

    return EXIT_SUCCESS;
}

