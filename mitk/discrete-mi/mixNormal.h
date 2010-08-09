/**
 *: mixNormal.h
*/


#ifndef MIX_NORMAL_H
#define MIX_NORMAL_H

#include <cstdio>
#include <fstream>
#include <string>
#include <map>
#include <vector>

#include "data_processing.h"
#include "config.h"
#include "error.h"

#ifndef DATA_TYPE_DEFINED
#error "DataType not defined"
#endif		/* ifndef DATA_TYPE_DEFINED */


#ifndef PROC_TYPE_DEFINED
#error "ProcType not defined"
#endif		/* ifndef PROC_TYPE_DEFINED */




//typedef Vector<ProcType> ProcVector;

using namespace std;

using std::ifstream;
using std::ofstream;


//////////////////// class MixNormal ////////////////////

// at the first stage, we just deal with diagonal covariance

/**
 * a gaussian mixture compoent
 */
class MixNormal {

public:
  // constructors and destructor
  MixNormal() { }
  MixNormal(unsigned numVars, unsigned dimX) { 
    _numVariables = numVars;
    _dX = dimX;
  }
  
  void setInit(unsigned numVars, unsigned dimX) { 
    _numVariables = numVars;
    _dX = dimX;
  }

  virtual ~MixNormal();

  void startEpoch(unsigned numVars, unsigned dimX);
  //bool addToEpoch(ProcVector *x, unsigned numData);
  bool addToEpoch(PtrArray<unsigned> *ptrArray, unsigned numData);
  double endEpoch(int numSamples);

private:
  
  map<vector<double>, int> xFreqTable;
  map<vector<double>, int> yFreqTable;
  map<vector<double>, int> xyFreqTable;
  /** the number of variables we are talking about */
  unsigned _numVariables;

  /** MI VARIABLES */
  /** size of X, the length of first split */
  unsigned _dX;

  double _pX, _pY, _pXY;

  /** information theoretic quantities */
  double _MI;

  /** number of samples used in LLN calculation by using "real" data */
  unsigned _nSamplesMI;
}; // end class MixNormal


#endif










