#define MAX_DIMENSIONALITY 20

#include "GMTK_ObservationMatrix.h"

template<class T>
class PtrArray {
 public:
  PtrArray() {
    
  }
  PtrArray(ObservationMatrix* O) {
    _obsMat = O;
    _stride = _obsMat->stride();
  }
  //return an array of elements pointed to by the array 
  //of pointers to the observation matrix
  T* toVec(unsigned frameno) {
    for(int i = 0; i < _size; ++i) {
      dataVec[i] = *(ptrVec[i] + _stride * frameno);
    }
    return dataVec;
  }

  T getVal(unsigned frameno, unsigned featPos) {
    return *(ptrVec[featPos] + _stride * frameno);
  }

  T getNextVal(unsigned featPos) {
    return *(currPtrVec[featPos] += _stride);
  }

  T* toCurrPtrVec(unsigned frameno) {
    for(int i = 0; i < _size; ++i) {
      currPtrVec[i] = (ptrVec[i] + _stride * frameno);
    }
    return currPtrVec;
  }
  
  T* nextPtrVec() {
    for(int i = 0; i < _size; ++i) {
      currPtrVec[i] = (currPtrVec[i] + _stride);
    }
    return currPtrVec;
  }
  
  void setPtrVec(unsigned i, T* pos) {
    ptrVec[i] = pos;
    currPtrVec[i] = pos; 
  }

  void setSize(int s) {
    _size = s;
  }

  T dataVec[MAX_DIMENSIONALITY];
  T* ptrVec[MAX_DIMENSIONALITY]; 
  T* currPtrVec[MAX_DIMENSIONALITY];
 private:
  int _size;
  unsigned _stride;
  ObservationMatrix* _obsMat;

};



