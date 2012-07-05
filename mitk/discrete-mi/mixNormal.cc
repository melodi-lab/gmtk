/**
 *: mixNormal.cc
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cfloat>

#include "rand.h"
#include "config.h"
#include "mixNormal.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

/**
 * methods for class MixNormal
 */


//////////////////// ~MixNormal ////////////////////

/**
 * destructor
 * cleanup the memory
 */
MixNormal::~MixNormal() {
  //cleanup();
} // end MixNormal


//////////////////// startEpoch ////////////////////

/**
 * add new epoch to the em training
 */
void MixNormal::startEpoch(unsigned numVars, unsigned dimX) {
  _numVariables = numVars;
  _dX = dimX;
} // end startEpoch


//////////////////// addToEpoch ////////////////////
/**
 *
 * @param data the Vector array of input data
 * @param numData the number of data
 * @return the learning status
 */
bool MixNormal::addToEpoch(PtrArray<unsigned> *ptrArray, unsigned numData) {

  vector<double> X(_dX), Y(_numVariables - _dX), XY(_numVariables); 
  //unsigned* data;
  for(unsigned i = 0; i < numData; i++) {  //go through all frames
    //data = ptrArray->toVec(i);
    for(unsigned j = 0; j < _dX; ++j) X[j] = XY[j] = ptrArray->getVal(i,j);
    for(unsigned j = _dX; j < _numVariables; ++j) 
      Y[j - _dX] = XY[j] = ptrArray->getVal(i,j);
    
    map<vector<double>,int>::iterator xfound  = xFreqTable.find(X);
    if(xfound != xFreqTable.end()) (*xfound).second++; 
    else xFreqTable[X] = 1;
    map<vector<double>,int>::iterator yfound  = yFreqTable.find(Y);
    if(yfound != yFreqTable.end()) (*yfound).second++; 
    else yFreqTable[Y] = 1;
    map<vector<double>,int>::iterator xyfound  = xyFreqTable.find(XY);
    if(xyfound != xyFreqTable.end()) (*xyfound).second++; 
    else xyFreqTable[XY] = 1;
  }

  return true;
} // end addToEpoch



//////////////////// endEpoch ////////////////////

/**
 * end session of an epoch. this will update all the parameters
 * normalize those alphas, copy the value to previous
 * and prepare for next epoch
 *
 * @return the learning status
 */
double MixNormal::endEpoch(int numSamples) {
  vector<double> X(_dX), Y(_numVariables - _dX), XY(_numVariables);
  double accumMI = 0.0, tmp;
  double invLog2 = 1.0 / log(2.0);
  for(map<vector<double>,int>::iterator it = xyFreqTable.begin();
      it != xyFreqTable.end();
      it++) 
    {
      //cout<<"(*it).second = "<<(*it).second<<endl;
      //cout<<"numSamples = "<<numSamples<<endl;
      _pXY = (double) (*it).second / numSamples;
      //cout<<"_pXY = "<<_pXY<<endl;
      XY = (*it).first;
      //cout<<"xyFreqTable = "<<xyFreqTable[XY]<<endl;
      for(unsigned j = 0; j < _dX; ++j) X[j] = XY[j];
      for(unsigned j = _dX; j < _numVariables; ++j) Y[j - _dX] = XY[j];
      _pX = (double) xFreqTable[X] / numSamples;
      //cout<<"xFreqTable = "<<xFreqTable[X]<<endl;
      //cout<<"_pX = "<<_pX<<endl;
      _pY = (double) yFreqTable[Y] / numSamples;
      //cout<<"_pY = "<<_pY<<endl;    
      tmp = log(_pXY) - log(_pX * _pY);
      tmp *= invLog2 * _pXY;
      accumMI += tmp;
    }
  
  //_MI = accumMI / numSamples;
  _MI = accumMI;
  cout<<"MI = "<<_MI<<endl;
  
  return _MI;

} // end endEpoch














