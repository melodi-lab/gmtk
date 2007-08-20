/*  Generated header
 *  File Name : GMTK_ObsNorm.cc
 *
 *  Created   : 2003-12-05 15:23:13 karim
 *  Author    : Karim Filali (karim@cs.washington.edu)
 *  Time-stamp: <>
*/

#ifndef __CYGWIN__
#include <values.h>
#endif

#include <math.h>
#include "GMTK_ObsNorm.h"
#include "general.h"

void obsNorm(FILE*  out_fp,
	     ObservationMatrix* obs_mat,
	     Range& srrng,
	     const double result_mean,
	     const double result_std,
	     char* segment_length_fname,
	     unsigned   segment_length,
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

  const size_t num_segments = obs_mat->numSegments();

  int* seg_markers = new int[num_segments];

  getSegMarkers(seg_markers,num_segments,segment_length_fname,segment_length);


  float *ftr_buf =  new float[buf_size * n_ftrs];
  float *oftr_buf = new float[buf_size * n_ftrs];
  
  UInt32* lab_buf = new UInt32[buf_size * n_labs];
  UInt32* olab_buf = new UInt32[buf_size * n_labs];
  
  double *const ftr_means = new double [n_ftrs];
  double *const ftr_stds  = new double [n_ftrs];
  
  double *ftr_means_p;
  double *ftr_stds_p;
  
  unsigned total_frames = 0;

  // new outer loop: for every group of segments, do this loop once.
  // thus, to imitate the behavior of pfile_normutts, this is done for
  // every sentence
  
  // Set 2 iterators to work in parallel: 1 to do mean/variance calcs,
  // other to do application of means and variances.
  int cur_seg_group;
  Range::iterator srit=srrng.begin();
  Range::iterator srit2=srrng.begin();
  for (cur_seg_group=seg_markers[(*srit)]; !srit.at_end();) {
    
    // Initialize the above declared arrays
    for (size_t i=0;i<n_ftrs;i++) {
      ftr_means[i] = ftr_stds[i] = 0.0;
    }
    total_frames=0;

    if (!quiet)
	printf("Computing means and variances from utt %d onward.\n",(*srit));
    
    // instead of going through entire range, only calculate ranges in 
    // same seg_group
    for (;!srit.at_end() && seg_markers[(*srit)]==cur_seg_group;srit++) {
      obs_mat->loadSegment(*srit);
      const size_t n_frames = obs_mat->numFrames();
   
      if (debug_level > 0) 
	printf("Processing sentence %d\n",(*srit));
      else if (!quiet && (*srit) % 100 == 0)
	printf("Processing sentence %d\n",(*srit));
      
      // Increase size of buffers if needed.
      if (n_frames > buf_size) {
	// Free old buffers.
	delete [] ftr_buf;
	delete [] oftr_buf;
	delete [] lab_buf;
	delete [] olab_buf;
	
	// Make twice as big to cut down on future reallocs.
	buf_size = n_frames * 2;
	
	// Allocate new larger buffers.
	ftr_buf = new float[buf_size * n_ftrs];
	oftr_buf = new float[buf_size * n_ftrs];
	lab_buf = new UInt32[buf_size * n_labs];
	olab_buf = new UInt32[buf_size * n_labs];
      }
      

      for(unsigned frame_no = 0;  frame_no < n_frames; ++frame_no) {
	const float* start_of_frame = obs_mat->floatVecAtFrame(frame_no);
	ftr_means_p = ftr_means;
	ftr_stds_p = ftr_stds;
	for(unsigned feat_no = 0;  feat_no < n_ftrs; ++feat_no) {
	  const double val = *(start_of_frame  + feat_no);
	  *ftr_means_p++ += val;
	  *ftr_stds_p++ += (val*val);
	}
      }

      total_frames += n_frames;
    }
    // update the seg group
    if (!srit.at_end())
      cur_seg_group=seg_markers[(*srit)];
   
    if (total_frames == 1) {
        printf("WARNING:: Only one frame used to compute the statistics.\n");
    }
    if (total_frames == 0) {
      error("ERROR, must have more than 0 frames in file");
    }
    
    ftr_means_p = ftr_means;
    ftr_stds_p = ftr_stds;
    for (unsigned i=0;i<n_ftrs;i++) {
      (*ftr_stds_p) = ( ( (*ftr_stds_p) - (*ftr_means_p)*(*ftr_means_p)/total_frames ) / total_frames );
      if (*ftr_stds_p < DBL_MIN) {
	error("ERROR, computed variance is too small %e\n",*ftr_stds_p);
      } else 
	*ftr_stds_p  = 1.0/sqrt(*ftr_stds_p);
      (*ftr_means_p) = (*ftr_means_p)/total_frames;
      ftr_means_p++;
      ftr_stds_p++;
    }
    
    if (!quiet) {
      if (srit.at_end())
	printf("Normalizing from utt %d to utt %d.\n",(*srit2),srrng.last());
      else
	  printf("Normalizing from utt %d to utt %d.\n",(*srit2),(*srit)-1);
    }
    for (;(*srit2)<(*srit) || (srit.at_end() && !srit2.at_end());srit2++) {
      obs_mat->loadSegment(*srit2);
      const size_t n_frames = obs_mat->numFrames();

      if (debug_level > 0) 
	printf("Processing sentence %d\n",(*srit2));
      else if (!quiet && (*srit2) % 100 == 0)
	printf("Processing sentence %d\n",(*srit2));
	
      // Increase size of buffers if needed.
      if (n_frames > buf_size) {
	// Free old buffers.
	delete ftr_buf;
	delete oftr_buf;
	delete lab_buf;
	delete olab_buf;
	
	// Make twice as big to cut down on future reallocs.
	buf_size = n_frames * 2;
	
	// Allocate new larger buffers.
	ftr_buf = new float[buf_size * n_ftrs];
	oftr_buf = new float[buf_size * n_ftrs];
	lab_buf = new UInt32[buf_size * n_labs];
	olab_buf = new UInt32[buf_size * n_labs];
      }

       
      float *oftr_buf_p = oftr_buf;
      UInt32* olab_buf_p = olab_buf;
      
      for(unsigned frame_no = 0;  frame_no < n_frames; ++frame_no) {
	float* start_of_frame = obs_mat->floatVecAtFrame(frame_no);
	UInt32* lab_buf_p =  obs_mat->unsignedAtFrame(frame_no);
	ftr_means_p = ftr_means;
	ftr_stds_p = ftr_stds;
	for(unsigned feat_no = 0;  feat_no < n_ftrs; ++feat_no) {
	  *oftr_buf_p++ = (*(start_of_frame  + feat_no) - *ftr_means_p)*(*ftr_stds_p)*result_std + result_mean;
	  ftr_means_p++;
	  ftr_stds_p++;
	}
	for(unsigned unsigned_feat_no = 0;  unsigned_feat_no < n_labs; ++unsigned_feat_no) {
	  *olab_buf_p++ = *(lab_buf_p+unsigned_feat_no);
	}
      }
    

      // Write output.
      printSegment(*srit2, out_fp, oftr_buf,n_ftrs,olab_buf,n_labs,n_frames, dontPrintFrameID,quiet, ofmt, debug_level, oswap, out_stream);
      
      //out_stream.write_ftrslabs(prrng.length(), oftr_buf, olab_buf);
      //out_stream.doneseg((SegID) seg_id);
    } // end inner loop... continue with next normalizing group
  } // end outer loop
  
  delete [] ftr_buf;
  delete [] lab_buf;
  
  delete [] ftr_means;
  delete [] ftr_stds;
  
  if(ofmt==PFILE) {
    delete out_stream;
  }
}



class tokFile {
    // A class that accesses a text file containing a list of 
    // white-space separated tokens, which may be interpreted as 
    // plain strings or as integers
public:
    FILE *file;
    char *filename;
    int buflen;
    char *buf;
    int bufpos;
    int linect;
    static const char *WS;

    tokFile(FILE* a_file, const char *a_name = "(unknown)") {
        file = a_file;
        filename = copyToNewStr(a_name);
        buflen = 1024;
        buf = new char[buflen];
        bufpos = 0;
        buf[bufpos] = '\0';
        linect = 0;
    }
    ~tokFile(void) {
        delete [] filename;
        delete [] buf;
    }
    int getNextTok(char **ret);
    int getNextInt(int *ret);
    // Shorthand to close our file (even though we didn't open it)
    void close(void) {
        fclose(file);
    }
};

const char *tokFile::WS = " \t\n";

int tokFile::getNextTok(char **ret) {
    // Pull another space-delimited token from an open file handle.
    // File needs just to be white-space
    // delimited, but can have comment lines starting with "#"
    // Pointer to token string written to *ret. 
    // Return 1 on success, 0 at EOF

    while (buf[bufpos] == '\0') {
        // need a new line
        char *ok = fgets(buf, buflen, file);
        if (!ok) {
            // EOF hit, nothing more got
            return 0;
        }
        int got = strlen(buf);
        while (buf[got-1] != '\n' && !feof(file)) {
            // because line didn't end with EOL, we assume we ran out 
            // of buffer space - realloc & keep going
	    // assert(got == buflen-1);
            int newbuflen = 2*buflen;
            char *newbuf = new char[newbuflen];
            memcpy(newbuf, buf, got);
            delete [] buf;
            buf = newbuf;
            buflen = newbuflen;
            fgets(buf+got, buflen, file);
            got = strlen(buf);
        }
        ++linect;
        // strip the trailing CR
        if (buf[got-1] == '\n') {
            buf[got-1] = '\0';
        }
        // OK, now we've got a new line, continue
        bufpos = strspn(buf, WS);
        // but if it's a comment line, skip it by pretending it's a blank line
        if (buf[bufpos] == '#') {
            buf[bufpos] = '\0';
        }
    }

    // Find the next space after the token
    // (dlbufpos already points to non-WS)
    int toklen = strcspn(buf+bufpos, WS);
    // I think this has to work
    //    assert(toklen > 0);

    // Save the result
    *ret = buf+bufpos;

    if (buf[bufpos+toklen] != '\0') {
        // There's more after the tok
        // terminate the string at terminal WS
        buf[bufpos+toklen] = '\0';      
        // Skip over the terminated, returned token
        bufpos += toklen+1;
        // Advance pointer to look at following non-WS next time
        bufpos += strspn(buf+bufpos, WS);
    } else {
        // This token is last in string - point to end of buffer for next time
        bufpos += toklen;
    }
    return 1;
}

int tokFile::getNextInt(int *ret) {
    // Read another number from the deslen file as the desired length
    // in frames of the next frame.  File needs just to be white-space
    // delimited, but can have comment lines starting with "#"
    // Next read value is put into *pdeslen.  Return 1 on success, 0 at EOF
    int rc = 0;
    // Get the next number in the line buffer & advance pointer
    char *str, *end;

    if ( (rc = getNextTok(&str)) ) {
        int val = strtol(str, &end, 10);

        // If unparseable, end points to start
        if (end == str) {
            fprintf(stderr, "unable to parse token '%s' as an integer "
                    "at line %d in file '%s'\n", str, linect, filename);
            return 0;
        }

        // Got the number
        *ret = val;
    }
    return rc;
}



void getSegMarkers(int* seg_markers, const unsigned size, char* segment_length_fname, unsigned segment_length) {

  int cur_seg=0;
  for (unsigned i=0; i<size; ++i) seg_markers[i] = 0;
  
  // count up total lengths: should match length of pfile
  if (segment_length_fname!=0) {
    FILE *seg_fp;
    tokFile *segTok;
    UInt32 segCount=0;
    int segLen;
    
    if ((seg_fp=fopen(segment_length_fname, "r"))==NULL)
      error("Couldn't open segment length file for reading");
    segTok=new tokFile(seg_fp,segment_length_fname);
    while(segTok->getNextInt(&segLen)) {
      if ((segCount+(UInt32)segLen)<=size) {
	for(UInt32 i=segCount;i<segCount+(UInt32)segLen;i++)
            seg_markers[i]=cur_seg;
	cur_seg++;
      }
      segCount+=(UInt32)segLen;
    }
      segTok->close();
      delete segTok;
      
      if (segCount!=size) {
        char errstr[1024];
        sprintf(errstr,"Mismatch between segment length file (%d utts) and file (%d utts)",segCount,size);
        error(errstr);
      }
  } else if (segment_length<INT_MAX) {
    for(UInt32 i=0;i<size;i+=segment_length) {
      for(UInt32 j=0;j<segment_length && j+i<size; j++)
	seg_markers[i+j]=cur_seg;
      cur_seg++;
      }
  }
  
}
