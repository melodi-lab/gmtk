/*  Generated header
 *  File Name : GMTK_Kmeans.cc

 Copied from the pfile directory.

 *
 *  Created   : 2003-12-15 21:19:59 karim
 *  Author    : Karim Filali (karim@cs.washington.edu)
 *  Time-stamp: <>
*/

#include "GMTK_Kmeans.h"

#include <string.h>
#include "error.h"
#include "rand.h"


RAND myrand(true);

int kmeans::kmeans_k = 5;
int kmeans::kmeans_vl = 5;

kmeans::kmeans(int _k,int vl)
  : k(_k), vector_length(vl)
{

  if (k <=0) {
    error("ERROR, can't have K=%d clusters\n",k);
  }

  cur_means = new float[k*vector_length];
  new_means = new float[k*vector_length];
  variances = new float[k*vector_length];
  
  saved_means = new float[k*vector_length];
  saved_variances = new float[k*vector_length];
  saved_counts = new int[k];

  done = false;
  randomAssignment = true;

  new_counts = new int[k];
  float *curp = cur_means;
  float *newp = new_means;
  float *varp = variances;
  for (int i=0;i<k;i++) {
    new_counts[i] = 0;
    for (int j=0;j<vector_length;j++) {
      *curp++ = drand48();
      *newp++ = 0.0;
      *varp++ = 0.0;
    }
  }

  //myrand = new RAND(true);

}

kmeans::~kmeans()
{
  delete [] cur_means;
  delete [] new_means;
  delete [] new_counts;
  delete [] variances;
  delete [] saved_means;
  delete [] saved_variances;
  delete [] saved_counts;

  //delete myrand;

}

void kmeans::initNew() {
  float *newp = new_means;
  float *varp = variances;
  for (int i=0;i<k;i++) {
    new_counts[i] = 0;
    for (int j=0;j<vector_length;j++) {
      *newp++ = 0.0;
      *varp++ = 0.0;
    }
  }
}

void kmeans::save() {
  ::memcpy((void*)saved_means,(void*)cur_means,
	   sizeof(float)*k*vector_length);
  ::memcpy((void*)saved_variances,(void*)variances,
	   sizeof(float)*k*vector_length);
  ::memcpy((void*)saved_counts,(void*)new_counts,
	   sizeof(int)*k);
}


inline float kmeans::distance(const float *const v1,const float *const v2) {
  float rc = 0;

  // assumes vector_length > 0
  const float *v1p = v1;
  const float *const v1_endp = v1+vector_length;
  const float *v2p = v2;
  do {
    const float tmp = (*v1p++ - *v2p++);
    rc += tmp*tmp;
  } while (v1p != v1_endp);


  // for (int i=0;i<vector_length;i++) {
  // const float tmp = (v1[i] - v2[i]);
  // rc += tmp*tmp;
  // }
  return rc;
}

inline void kmeans::add2new(const int lk,const float *const v) {
  new_counts[lk]++;
  float *k_mean = &new_means[lk*vector_length];

  // assumes vector_length > 0
  const float *vp = v;
  const float *const v_endp = v+vector_length;
  float *k_meanp = k_mean;
  do {
    *k_meanp++ += *vp++;
  } while (vp != v_endp);

  // for (int i=0;i<vector_length;i++) {
  // k_mean[i] += v[i];
  // }
}

void kmeans::add2new(const float *const v) {
  float *cur_meansp = cur_means;
  float md = distance(cur_meansp,v);
  int inx = 0;

  cur_meansp += vector_length;
  for (int i=1;i<k;i++) {
    const float tmp = distance(cur_meansp,v);
    if (tmp < md) {
      md = tmp;
      inx = i;
    }
    cur_meansp += vector_length;
  }
  add2new(inx,v);
}

void kmeans::add2newRand(const float *const v) {
  //  int inx = myrand->uniform(k-1);
  int inx = myrand.uniform(k-1);
  add2new(inx,v);
}

void kmeans::computeVariances(const float *const v) {

  // first compute the mean this vector is closest to:
  float *cur_meansp = cur_means;
  float md = distance(cur_meansp,v);
  int inx = 0;
  int i;
  cur_meansp += vector_length;
  for (i=1;i<k;i++) {
    float tmp = distance(cur_meansp,v);
    if (tmp < md) {
      md = tmp;
      inx = i;
    }
    cur_meansp += vector_length;
  }

  // mean this vector is closest to is inx.
  cur_meansp = &cur_means[inx*vector_length];
  float *variancesp = &variances[inx*vector_length];
  for (i=0;i<vector_length;i++) {
    const float tmp = v[i] - cur_meansp[i];
    variancesp[i] += tmp*tmp;
  }
}


double kmeans::finishVariances() {
  
  double sum=0.0;
  // mean this vector is closest to is inx.
  float *variancesp = variances;
  for (int i=0;i<k;i++) {
    const float norm = 1.0/new_counts[i];
    for (int j=0;j<vector_length;j++) {
      variancesp[j] *= norm;
      sum += variancesp[j];
    }
    variancesp += vector_length;
  }
  // return sum of variances.
  return sum;
}


bool kmeans::someClusterHasLessThanNEntries(int n) {
  for (int i=0;i<k;i++)
    if (new_counts[i] <n)
      return true;
  return false;
}


bool kmeans::someClusterHasZeroEntries() {
  for (int i=0;i<k;i++)
    if (new_counts[i] == 0)
      return true;
  return false;
}




// return true if no samples were 
// given to this kmeans object.
bool kmeans::zeroCounts() {
  for (int i=0;i<k;i++)
    if (new_counts[i] != 0)
      return false;
  return true;
}



void kmeans::finishNew() {
  float *newp = new_means;
  for (int i=0;i<k;i++) {
    double inv_count = 1.0/new_counts[i];
    for (int j=0;j<vector_length;j++) {
      *newp *= inv_count;
      newp++;
    }
  }
}


float kmeans::newCurDist() {

  float *curp = cur_means;
  float *newp = new_means;
  float dist = 0.0;
  for (int i=0;i<k;i++) {
    dist += distance(curp,newp);
    curp += vector_length;
    newp += vector_length;
  }
  return dist;
}


void kmeans::printSaved(FILE *fp) {
  float *meansp = saved_means;
  float *variancesp = saved_variances;
  for (int i=0;i<k;i++) {
    int j;
    fprintf(fp,"%d means(%d): ",i,saved_counts[i]);
    for (j=0;j<vector_length;j++) {
      fprintf(fp,"%0.5e ",meansp[j]);
    }
    fprintf(fp,"\n");
    fprintf(fp,"%d varns: ",i);
    for (j=0;j<vector_length;j++) {
      fprintf(fp,"%0.5e ",variancesp[j]);
    }
    fprintf(fp,"\n");

    meansp += vector_length;
    variancesp += vector_length;
  }
}

void kmeans::writeMgDoubleRecord2D(FILE *stream) {

  // only writes the first 2

  if (vector_length != 2) 
    error("kmeans::writeMgDoubleRecord2D, vl = %d, this function must have length 2 vectors\n",vector_length);

  char char_k = kmeans_k;
  fwrite(&char_k,sizeof(char_k),1,stream);

  int i;
  int total_count = 0;
  for (i=0;i<kmeans_k;i++) {
    total_count += saved_counts[i];
  }
  float *meansp = saved_means;
  float *varsp = saved_variances;

  for (i=0;i<kmeans_k;i++) {
    double tmp;

    // mean x
    tmp = meansp[0];
    fwrite(&tmp,sizeof(tmp),1,stream);

    // mean y
    tmp = meansp[1];
    fwrite(&tmp,sizeof(tmp),1,stream);

    // var x
    tmp = varsp[0];
    fwrite(&tmp,sizeof(tmp),1,stream);

    // cov xy
    tmp = 0.0;
    fwrite(&tmp,sizeof(tmp),1,stream);

    // var y
    tmp = varsp[1];
    fwrite(&tmp,sizeof(tmp),1,stream);

    // coef
    tmp = saved_counts[i]/(double)total_count;
    fwrite(&tmp,sizeof(tmp),1,stream);

    meansp += vector_length;
    varsp += vector_length;
  }
}
