/*  Generated header
 *  File Name : GMTK_DataTransformations.cc
 *
 *  Created   : 2003-09-05 18:42:32 karim
 *  Author    : Karim Filali (karim@cs.washington.edu)
 *  Time-stamp: <>
 *
 * $Header$
 *
*/


#include "GMTK_ObservationMatrix.h"

//#define DEBUG

#ifdef DEBUG
#define DBGFPRINTF(_x_) fprintf _x_
#else
#define DBGFPRINTF(_x_)
#endif

/**
 * perform in place mean and variance normalization sentence wise:
 * E[X] = (\sigma X) / num_frames
 * VAR[X]=|E[X^2]-E[X]^2|*(num_frames/(num_frames-1)
 * X = (X - E[X]) / sqrt(VAR[X])
 *---------------------------------------------------------------*/
void ObservationMatrix::inPlaceMeanSubVarNorm(float* x, int vec_size, int stride, int num_frames) {

#define TOO_SMALL (1e-10)

  assert(num_frames > 0);
  float mean_squared,tmp;
  float* mean = new float[vec_size];
  float* var  = new float[vec_size];
  float* x_ptr=x;

  for (int i = 0; i < vec_size; i++){
    mean[i] = 0.0;
    var[i]  = 0.0;
    mean_squared = 0.0;
    for (int j = 0; j < num_frames; j++){
      tmp =  *(x_ptr + i + j*stride);
      mean[i] += tmp;
      mean_squared += tmp*tmp;
    }
    mean[i] /= num_frames;
    mean_squared /= num_frames;
    var[i] = sqrt(num_frames/(num_frames - 1.0)*fabs(mean_squared - mean[i]*mean[i]));
  }

  for (int i = 0; i < vec_size; i++)
    for (int j = 0; j < num_frames; j++){
       *(x_ptr + i + j*stride) -= mean[i];
       if (var[i] > TOO_SMALL)
         *(x_ptr + i + j*stride) /= var[i];
       else
	 warning("WARNING: Value of variance in mean/variance normalization is too small (<=%f).  Will only perform mean substraction.\n", TOO_SMALL);
    }
}

/**
 * perform in place mean substraction sentence wise:
 * E[X] = (\sigma X) / num_frames
 * X = (X - E[X])
 *---------------------------------------------------------------*/
void ObservationMatrix::inPlaceMeanSub(float* x, int vec_size, int stride, int num_frames) {

#define TOO_SMALL (1e-10)

  assert(num_frames > 0);
  float tmp;
  float* mean = new float[vec_size];
  float* x_ptr=x;

  for (int i = 0; i < vec_size; i++){
    mean[i] = 0.0;
    for (int j = 0; j < num_frames; j++){
      tmp =  *(x_ptr + i + j*stride);
      mean[i] += tmp;
    }
    mean[i] /= num_frames;
  }

  for (int i = 0; i < vec_size; i++)
    for (int j = 0; j < num_frames; j++){
      *(x_ptr + i + j*stride) -= mean[i];
    }
}

/**
 * Read a filter from a file
 *
 * Pre-conditions:  filter_file_name needs to have storage alloacted to it already
 * Side effects:    changes filter_len
 */
void ObservationMatrix::readFilterFromFile(char* filter_file_name,float* filter_coeffs, unsigned& filter_len) {

  float tmp;

  FILE* ifp;
  if((ifp=fopen(filter_file_name,"r")) == NULL) {
    error("ERROR: Could not open file \'%s\' for reading.\n",filter_file_name);
  }

  unsigned i=0;
  while(fscanf(ifp,"%f",&tmp) !=EOF) {
    if( i > MAX_FILTER_LEN ) {
      error(" ObservationMatrix::readFilterFromFile: Filter lenght exceeds maximum allowable length (%d)",MAX_FILTER_LEN);
    }
    filter_coeffs[i++]=tmp;
  }

  filter_len=i;
  
}

/**
 * apply a filter
 *
 *
 */
void ObservationMatrix::filter(float* x, unsigned vec_size, unsigned stride, unsigned num_frames,float* filter_coeffs, unsigned filter_len) {

 if(vec_size==0 || stride ==0)
    return;

 assert(num_frames>0);

 if(filter_len > num_frames) {
   error("ERROR: ObservationMatrix::filter: The filter length (%d) is greater than the total number of frames (%d) in sentence.\n",filter_len,num_frames);
 }

 if(filter_len < MIN_FILTER_LEN) {
   warning("WARNING: ObservationMatrix::filter: filter length is less than %d\n", MIN_FILTER_LEN);
 }

  float* tmp_buf=new float[num_frames*vec_size];

  double sum=0.0;
  unsigned half_filter_len=filter_len/2;
  unsigned window_sample;

  // copy buffer into a temporary one
  for(unsigned i=0;i<num_frames;++i)  
    for(unsigned j=0;j<vec_size;++j) {
      tmp_buf[i*vec_size+j]=x[i*stride+j];
    }
  
  // Use leading zeros at the beginning of the data buffer
  for (unsigned sample = 0; sample < half_filter_len ; ++sample)     {
    sum = 0.0;
    window_sample= sample - half_filter_len;
    for (unsigned i = 0; i < filter_len ;++i) 
         sum += filter_coeffs[i] * ((window_sample<0)? 0.0 : tmp_buf[window_sample++]) ;
    // Do something about saturation?
    x[sample] = (float)sum;
  }
  
  for (unsigned sample =  half_filter_len; sample < num_frames -  filter_len +  half_filter_len ; ++sample ) {
    sum = 0.0;
    window_sample= sample -  half_filter_len;
    for (unsigned i = 0; i < filter_len ;i++)
      sum += filter_coeffs[i]*tmp_buf[window_sample++];
    x[sample] = (float)sum;
  }
  
  // end
  for (unsigned sample = num_frames - filter_len +  half_filter_len; sample < num_frames ; ++sample) {
    sum = 0.0;
    window_sample= sample - half_filter_len;
    for (unsigned i = 0; i < filter_len ;i++)
      sum += filter_coeffs[i] * tmp_buf[window_sample++] ;
    x[sample] = (float) sum;
  }
}


/**
 * apply an ARMA filter in place  (from Chia-Ping)
 *
 *
 */
void ObservationMatrix::arma(float* x, unsigned vec_size, unsigned stride, unsigned num_frames,unsigned order) {
  if(vec_size==0 || stride ==0)
    warning("WARNING: Not applying ARMA filter because number of floats is zero.");
  //    return;
  assert(num_frames>0);
  
  for(unsigned i=0;i<vec_size;i++) {
    for(unsigned frame_no=order;frame_no<num_frames-order;frame_no++) {
      float tf = x[stride*frame_no+i];
      for(unsigned tau=1;tau<=order;tau++) {
        tf += (x[stride*(frame_no-tau)+i] + x[stride*(frame_no+tau)+i]);
      }
      x[stride*frame_no+i] = tf/(2*order+1);
   }
  }
  
}

/**
 * multiply each component of the feature matrix by the supplied multiplier
 *
 */
void ObservationMatrix::multiply(float* x, unsigned vec_size, unsigned stride, unsigned num_frames, double multiplier) {
  for (unsigned i = 0; i < num_frames; i++) {
    for (unsigned j = 0; j < vec_size; j++) {
      *(x + i*stride + j) *= multiplier;
    }
  }
}


/**
 * add an offset to  each component of the feature matrix by the supplied multiplier
 *
 */
void ObservationMatrix::addOffset(float* x, unsigned vec_size, unsigned stride, unsigned num_frames, double offset) {
  for (unsigned i = 0; i < num_frames; i++) {
    for (unsigned j = 0; j < vec_size; j++) {
      *(x + i*stride + j) += offset;
    }
  }
}
