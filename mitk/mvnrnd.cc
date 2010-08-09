/*  Generated header
 *  File Name : mvnrnd.cc
 *
 *  Created   : 2004-01-09 15:31:20 karim
 *  Author    : Karim Filali (karim@cs.washington.edu)
 *  Time-stamp: <>
*/

#include <cstdio>
#include <cstring>
#include "matrix-ops.h"
#include "rand.h"

#define MAX_DIM 20
#define MAX_NUM_LEN 50

RAND rnd(true); // true to seed the rand num generator; false, otherwise.

unsigned parseMean(char*, double*);
unsigned parseCov(char*, double*);
int parseNumber(char*& pStr, double& num);

bool Print_Chol=false;
bool Print_Rand_Vec=false;

int main(int argc, char* argv[]) {

  if(argc<4) {
    fprintf(stderr,"Need three arguments: mean vector, sigma matrix, and number of samples.\n");
    exit(-1);
  }

  //  srand(time(NULL));

  char * meanStr = argv[1];
  char * covStr = argv[2];
  unsigned cases = atol(argv[3]);
  
  double mean[MAX_DIM];
  double cov[MAX_DIM*MAX_DIM];
  
  unsigned dim = parseMean(meanStr,mean);
  unsigned covDim = parseCov(covStr,cov);


  if(dim != covDim) {
    fprintf(stderr,"Dimensions of covariance matrix and mean vectr do not much.\n"); exit(-1);
  }
  
#ifdef DEBUG
  for(unsigned i=0;i<dim;++i) {
    fprintf(stdout,"%f ",mean[i]);
  }
  fprintf(stdout,"\n");
#endif

  double* chol = new double[covDim*covDim];

  if(choleskyDecomp(cov,chol,covDim)) {
    fprintf(stderr,"The covariance matrix is not positive definite.\n");
    exit(-1);
  }


  if(Print_Chol) {
    for(unsigned i=0;i<dim;++i) {
      for(unsigned j=0;j<dim;++j) {
	fprintf(stdout,"%f ",chol[i*dim+j]);
      }
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
  }

  double* randVec = new double[dim];
  for(unsigned sample=0; sample<cases; ++sample) {
    for(unsigned i=0; i<dim; ++i) {
      //randNormal = inverse_normal_func((double)rand()/RAND_MAX);
      randVec[i] = rnd.normal();
    }  
    if(Print_Rand_Vec) {
      fprintf(stdout,"( ");
      for(unsigned i=0; i<dim; ++i) {
	fprintf(stdout,"%.4f ", randVec[i]);
      }  
      fprintf(stdout,")\n");
    }
    vecMatProduct(randVec,chol,dim);  // randVec now contains the product
    for(unsigned i=0; i<dim; ++i) {
      randVec[i] += mean[i];
      fprintf(stdout,"%.4f ",randVec[i]);
    }
    fprintf(stdout,"\n");
  }

  

  delete [] chol;

  return 0;
}


#define END_INPUT 0
#define MORE_INPUT 1
#define NEW_ROW 2

unsigned parseMean(char* meanStr, double* mean) {

  char* pStr = meanStr;

  unsigned i=0;
  while(parseNumber(pStr,mean[i])==MORE_INPUT) i++;
  if(i==0) {
        fprintf(stderr,"Parse error: no number found in the mean vector\n"); exit(-1);
  }
  
  return i;  // dimension of mean vector

}



unsigned parseCov(char* covStr, double* cov) {
  char*    pStr = covStr;

  double tmp[MAX_DIM];

  unsigned numRows     = 0;
  unsigned prevNumCols = 0;
  unsigned numCols     = 0;
  unsigned ret         = NEW_ROW;
  while(ret==NEW_ROW) {
    prevNumCols=numCols;
    unsigned i=0;
    while((ret=parseNumber(pStr,tmp[i]))==MORE_INPUT) i++;
    if(i==0) {
    fprintf(stderr,"Parse error: no number found in the mean vector\n"); exit(-1);
    }
    numCols=i;
    if(numRows > 0 && prevNumCols != numCols) {
      fprintf(stderr,"Unequal number of entries in covariance matrix.\n"); exit(-1);
    }
    for(unsigned j=0; j< numCols; ++j) {
      cov[numRows*numCols + j]=tmp[j];
    }
    numRows++;
  }

  if(numRows != numCols) {
    fprintf(stderr,"Unequal num rows and columns in covariance matrix.\n"); exit(-1);
  }
  
  return numCols;  // dimension of cov
}


int parseNumber(char*& pStr, double& num) {
  char *pEnd;

  pStr += strspn(pStr," \t,"); // skip white space and commas
  if(*pStr == '\0') return END_INPUT;
  num = strtod(pStr,&pEnd);
  if(num==0 && pEnd==pStr) {  // no conversion was made
    if(*pStr==';') { pStr++; return NEW_ROW; }
    fprintf(stderr,"Parse error while reading number in the mean vector\n"); exit(-1);
  }

  pStr=pEnd;
  return MORE_INPUT;
}
