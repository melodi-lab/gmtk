/*  Generated header
 *  File Name : GMTK_ObsKLT.cc
 *
 *  Created   : 2003-12-03 11:59:48 karim
 *  Author    : Karim Filali (karim@cs.washington.edu)
 *  Time-stamp: <>

 Computes a Karhunen-Loeve transformation on the feature vectors. 

 Originally written by Jeff Bilmes <bilmes@ee.washington.edu>.
 Modified to integrate it with the more general observation tools.

*/

#include "GMTK_ObsKLT.h"
extern "C" {
#include "eig.h"
}

extern "C" void mul_mdmd_md(const int M, const int K, const int N, 
		       const double *const A, const double *const B, double *const C, 
		       const int Astride, const int Bstride, const int Cstride);

void readStats(FILE*f, size_t N, bool ascii, double *cor, double *means, double *vecs, double *vals) {
  size_t i,j;

  for (i=0;i<N;i++) {
    if (!ascii) {
      for (j=0;j<N;j++) 
	if (!fread(&cor[i*N+j],sizeof(cor[0]),1,f))
	  error("input stat file eof error, cor(%d,%d)\n",i,j);
      if (!fread(&means[i],sizeof(means[0]),1,f))
	  error("input stat file eof error, mean(%d)\n",i);
      for (j=0;j<N;j++) 
	if (!fread(&vecs[i*N+j],sizeof(vecs[0]),1,f))
	  error("input stat file eof error, vecs(%d,%d)\n",i,j);
      if (!fread(&vals[i],sizeof(vals[0]),1,f))
	  error("input stat file eof error, vals(%d)\n",i);
    } else {
      for (j=0;j<N;j++) 
	if (!fscanf(f,"%lf",&cor[i*N+j]))
	  error("input stat file eof error, cor(%d,%d)\n",i,j);
      if (!fscanf(f,"%lf",&means[i]))
	  error("input stat file eof error, mean(%d)\n",i);
      for (j=0;j<N;j++) 
	if (!fscanf(f,"%lf",&vecs[i*N+j]))
	  error("input stat file eof error, vecs(%d,%d)\n",i,j);
      if (!fscanf(f,"%lf",&vals[i]))
	  error("input stat file eof error, vals(%d)\n",i);
    }
  }
  // Could check for eof condition to be
  // true here.
}

void writeStats(FILE*f, size_t N, bool ascii, double *cor, double *means, double *vecs, double *vals) {
  size_t i,j;

  for (i=0;i<N;i++) {
    if (ascii) {
      for (j=0;j<N;j++)
	fprintf(f,"%.15f ",cor[i*N+j]);
      fprintf(f,"%.15f ",means[i]);
      for (j=0;j<N;j++)
	fprintf(f,"%.15f ",vecs[i*N+j]);
      fprintf(f,"%.15f\n",vals[i]);
    } else {
      fwrite(&cor[i*N],sizeof(double),N,f);
      fwrite(&means[i],sizeof(double),1,f);
      fwrite(&vecs[i*N],sizeof(double),N,f);
      fwrite(&vals[i],sizeof(double),1,f);
    }
  }
}

void obsKLT(FILE* out_fp, ObservationMatrix* obs_mat, FILE *in_st_fp,FILE *out_st_fp, Range& srrng,Range& ofrrng,const bool unity_variance,const bool ascii,const bool dontPrintFrameID,const bool quiet,unsigned ofmt,int debug_level,bool oswap) {


  // Feature and label buffers are dynamically grown as needed.
  size_t       buf_size     = 300;      // Start with storage for 300 frames.
  const size_t n_labs = obs_mat->numDiscrete();
  const size_t n_ftrs = obs_mat->numContinuous();
  size_t       max_n_frames = 0;

  OutFtrLabStream_PFile* out_stream=NULL;
  if(ofmt==PFILE) {
    out_stream = new OutFtrLabStream_PFile(debug_level,"",out_fp,n_ftrs,n_labs,1,oswap);
  }

  
  float *      ftr_buf  = new float[buf_size * n_ftrs];
  float *      ftr_buf_p;
  UInt32*      lab_buf  = new UInt32[buf_size * n_labs];
  
  size_t       total_frames = 0;
  
  // mean vector E[X}
  double *     const ftr_means = new double [n_ftrs];
  double *     ftr_means_p;
  double *     const ftr_means_endp = ftr_means+n_ftrs;
  double *     ftr_cov;
  
  double * ftr_eigenvecs;
  double * ftr_eigenvals;
  
  // Initialize the above declared arrays
  size_t i,j;
  memset(ftr_means,0,n_ftrs*sizeof(double));
  
  
  if (in_st_fp == NULL) {
    // correlation vector (i.e., E[X X^T]
    double *const ftr_cor = new double [n_ftrs*(n_ftrs+1)];
    double *ftr_cor_p;
    memset(ftr_cor,0,n_ftrs*(n_ftrs+1)*sizeof(double)/2);
    
    printf("Computing feature means and covariance..\n");

    // Go through input pfile to get the initial statistics,
    for (Range::iterator srit=srrng.begin();!srit.at_end();srit++) {

      obs_mat->loadSegment((const unsigned)(*srit));
      const size_t n_frames = obs_mat->numFrames();
      
      if ((*srit) % 100 == 0)
	printf("Processing sentence %d\n",(*srit));
      
      if (n_frames > max_n_frames)
	max_n_frames = n_frames;
      // Increase size of buffers if needed.
      if (n_frames > buf_size) {
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

      float *ftr_buf_base = ftr_buf;
      for (i=0;i<n_frames;i++) {
	size_t j;
	ftr_means_p = ftr_means;
	ftr_cor_p = ftr_cor;
	ftr_buf_p = ftr_buf_base;
	for (j=0;j<n_ftrs;j++) {
	  *ftr_means_p += *ftr_buf_p;
	  
	  register double tmp = *ftr_buf_p;
	  
	  // Only compute the upper triangular part
	  // since matrix is symmetric.
	  float *ftr_buf_pp = ftr_buf_base+j;
	  for (size_t k=j;k<n_ftrs;k++) {
	    *ftr_cor_p++ += (tmp)*(*ftr_buf_pp++);
	  }
	  
	  ftr_means_p++;
	  ftr_buf_p++;
	}
	ftr_buf_base += n_ftrs;
      }
      
      
      
      total_frames += n_frames;
    }
    
    ftr_cov = new double [n_ftrs*n_ftrs];
    double total_frames_inv = 1.0/total_frames;
    // actually compute the means and covariances.
    
    ftr_means_p = ftr_means;
    while (ftr_means_p != ftr_means_endp)
      (*ftr_means_p++) *= total_frames_inv;
    
    ftr_means_p = ftr_means;
    ftr_cor_p = ftr_cor;
    for (i=0;i<n_ftrs;i++) {
      double *ftr_means_pp = ftr_means_p;
      double *ftr_cov_rp = ftr_cov + i*n_ftrs+i; // row pointer
      double *ftr_cov_cp = ftr_cov_rp; // column pointer
      for (j=i;j<n_ftrs;j++) {
	double tmp;
	tmp = (*ftr_cor_p++)*total_frames_inv - 
	  (*ftr_means_p)*(*ftr_means_pp);
	*ftr_cov_rp = *ftr_cov_cp = tmp;
	
	ftr_cov_rp++;
	ftr_cov_cp += n_ftrs;
	ftr_means_pp++;
      }
      ftr_means_p++;
    }
    
    // don't need any longer
    delete ftr_cor;
    
    // now compute the eigen vectors and values
    ftr_eigenvecs = new double[n_ftrs*n_ftrs];
    ftr_eigenvals = new double[n_ftrs];
    eigenanalyasis(n_ftrs,ftr_cov,ftr_eigenvals,ftr_eigenvecs);
    
    // save it
    if (out_st_fp != NULL) 
      writeStats(out_st_fp,n_ftrs,ascii,
		 ftr_cov,ftr_means,
		 ftr_eigenvecs,ftr_eigenvals);
    
  } else {
    // read in the matrix containing the data.
    ftr_cov       = new double [n_ftrs*n_ftrs];
    ftr_eigenvecs = new double [n_ftrs*n_ftrs];
    ftr_eigenvals = new double [n_ftrs];
    readStats(in_st_fp, n_ftrs, ascii, ftr_cov,ftr_means, ftr_eigenvecs,ftr_eigenvals);
    
    
    // God knows why, but if the user
    // ask for this, save it
    if (out_st_fp != NULL) 
      writeStats(out_st_fp,n_ftrs,ascii, ftr_cov,ftr_means, ftr_eigenvecs,ftr_eigenvals);
    }
  
  
  // at this point we no longer need the following
  delete ftr_cov;

  
    
    if (unity_variance) {
      // multily in the eigenvalues into the eigenvectors
      double *valp = ftr_eigenvals;
      double *vecp = ftr_eigenvecs;
      for (i=0;i<n_ftrs;i++) {
	double *vecpp = vecp;
	double valp_inv = sqrt(1.0/(*valp));;
	for (j=0;j<n_ftrs;j++) {
	  (*vecpp) *= valp_inv;
	  vecpp+=n_ftrs;
	}
	valp++;
	vecp++;
      }
    }
    
    
    // allocate space performing the orthogonalization.
    if (max_n_frames > 0)
      buf_size = max_n_frames;
    
    float *oftr_buf = new float[buf_size * ofrrng.length()];
    float *oftr_buf_p;
    double *ftr_dbuf_src = new double[buf_size*n_ftrs];
    double *ftr_dbuf_src_p;
    double *ftr_dbuf_src_endp;
    double *ftr_dbuf_dst = new double[buf_size*n_ftrs];
    double *ftr_dbuf_dst_p;
    
    for (Range::iterator srit=srrng.begin();!srit.at_end();srit++) {
      obs_mat->loadSegment(*srit);
      const size_t n_frames = obs_mat->numFrames();
      
      if ((*srit) % 100 == 0)
	printf("Processing sentence %d\n",(*srit));
      
      // Increase size of buffers if needed.
      if (n_frames > buf_size) {
	// Free old buffers.
	delete ftr_buf;
	delete lab_buf;
	delete oftr_buf;
	delete ftr_dbuf_src;
	delete ftr_dbuf_dst;
	
	// Make twice as big to cut down on future reallocs.
	buf_size = n_frames * 2;
	
	// Allocate new larger buffers.
	ftr_buf = new float[buf_size * n_ftrs];
	lab_buf = new UInt32[buf_size * n_labs];
	oftr_buf = new float[buf_size * ofrrng.length()];
	ftr_dbuf_src = new double[buf_size*n_ftrs];
	ftr_dbuf_dst = new double[buf_size*n_ftrs];
	
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


      // While subtracting off the means, convert the features to 
      // doubles and do the normalization there.
      ftr_dbuf_src_p = ftr_dbuf_src;
      ftr_dbuf_src_endp = ftr_dbuf_src + n_frames*n_ftrs;
      ftr_buf_p  = ftr_buf;
      while (ftr_dbuf_src_p != ftr_dbuf_src_endp) {
	ftr_means_p = ftr_means;
	while (ftr_means_p != ftr_means_endp)
	  *ftr_dbuf_src_p++ = (double) *ftr_buf_p++ - *ftr_means_p++;
      }


      // Normalize the features with matrix multiply.
      // mul_mdmd_md(const int M, const int K, const int N, 
      //             const double *const A, const double *const B, double *const C, 
      //             const int Astride, const int Bstride, const int Cstride)
      
      mul_mdmd_md(n_frames,n_ftrs,n_ftrs,
		  ftr_dbuf_src,ftr_eigenvecs,ftr_dbuf_dst,
		  n_ftrs,n_ftrs,n_ftrs);
      
      
      ftr_dbuf_dst_p = ftr_dbuf_dst;
      oftr_buf_p = oftr_buf;
      for (i=0;i<n_frames;i++) {
	for (Range::iterator frit=ofrrng.begin();!frit.at_end(); ++frit) {
	  *oftr_buf_p++ = (float)ftr_dbuf_dst_p[*frit];
	  }
	ftr_dbuf_dst_p += n_ftrs;
      }
      
      // Write output.
      printSegment(*srit, out_fp, oftr_buf,n_ftrs,lab_buf,n_labs,n_frames, dontPrintFrameID,quiet, ofmt, debug_level, oswap, out_stream);
      
      //out_stream->write_ftrslabs(n_frames, oftr_buf, lab_buf);
      //out_stream->doneseg((SegID) *srit);
      
    }
    
    delete oftr_buf;
    delete ftr_dbuf_src;
    delete ftr_dbuf_dst;

    
  delete ftr_buf;
  delete lab_buf;
  delete ftr_means;
  delete ftr_eigenvecs;
  delete ftr_eigenvals;

  if(ofmt==PFILE) {
    delete out_stream;
  }

}

