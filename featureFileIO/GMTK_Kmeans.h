#ifndef GMTK_KMEANS_H
#define GMTK_KMEANS_H

#include <cstdio>
#include <cstdlib>
#include "rand.h"

/////////////////////////////////////////////////////////////////////////////
// kmeans support class
/////////////////////////////////////////////////////////////////////////////


class kmeans {

  const int k;
  const int vector_length;

  // k*vl matrices
  float *saved_means;
  float *saved_variances;
  int   *saved_counts;

  float *cur_means;
  float *new_means;  
  int *new_counts;
  
  float *variances;

  //RAND* myrand;

  float distance(const float *const v1,
		 const float *const v2);
  // add vector v to count against new
  void add2new(const int lk,const float *const v);

public:
  
  static int kmeans_k;
  static int kmeans_vl;

  kmeans(int _k=kmeans_k, int vl=kmeans_vl);
  ~kmeans();

  void initNew();
  void add2new(const float *const v);
  void add2newRand(const float *const v);
  bool someClusterHasZeroEntries();
  bool someClusterHasLessThanNEntries(int n);
  bool zeroCounts();
  void finishNew();

  void save();

  void computeVariances(const float *const v);
  double finishVariances();

  bool done;
  bool randomAssignment;

  void printSaved(FILE *fp);

  // avg distance between new and cur
  float newCurDist();
  // swap the new and current parameters.
  void swapCurNew() { float *tmp=cur_means;cur_means=new_means;new_means=tmp; }

  void writeMgDoubleRecord2D(FILE *stream);

};

#endif
