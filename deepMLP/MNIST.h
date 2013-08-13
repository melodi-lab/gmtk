#pragma once

#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>

#include "rand.h"

using namespace std;

class MNISTData {
  friend class Iterator;

  AllocatingMatrix _images;
  AllocatingMatrix _testImages;
  vector<unsigned char> _labels;
  int _numPixels, _numRows, _numCols, _numFeatures;

  AllocatingMatrix _trans, _invTrans;
  AllocatingVector _means;

  static const int MAX_PIXEL = 255;

  unsigned int SwapEndian(unsigned int x) {
    char * b = (char*) &x;
    unsigned int ret;
    char* br = (char*) &ret;
    br[0] = b[3];
    br[1] = b[2];
    br[2] = b[1];
    br[3] = b[0];
    return ret;
  }

public:
  int NumRows() const { return _numRows; }
  int NumCols() const { return _numCols; }

  MNISTData(const string & images, const string & labels) {

    ifstream imageStream(images, ios::in|ios::binary);
    if (!imageStream.is_open()) {
      cout << "Couldn't open feature file " << images.c_str() << endl;
      abort();
    }

    unsigned int header[4];
    imageStream.read((char*)header, 16);

    int magic = SwapEndian(header[0]);
    int numImages = SwapEndian(header[1]);
    _numRows = SwapEndian(header[2]);
    _numCols = SwapEndian(header[3]);

    _numPixels = _numRows * _numCols;

    _images.Resize(_numPixels, numImages);

    vector<unsigned char> charImage(_numPixels);
    for (int i = 0; i < numImages; ++i) {
      imageStream.read((char*)&charImage[0], _numPixels);
      MutableVector image = _images.GetCol(i);
      for (int j = 0; j < _numPixels; ++j) {
        double val01 = (double)charImage[j] / MAX_PIXEL;
        image[j] = 2 * val01 - 1.0;
      }
    }

#ifdef NDEBUG
    int reductionFactor = 1;
#else
    int reductionFactor = 1;
#endif

    if (reductionFactor > 1) {
      // make small
      cout << "WARNING!!! Reducing size for testing purposes!!!" << endl;
      int newNumRows = _numRows / reductionFactor;
      int newNumCols = _numCols / reductionFactor;
      int newNumPixels = newNumRows * newNumCols;
      AllocatingMatrix newImages(newNumPixels, numImages);
      for (int i = 0; i < numImages; ++i) {
        MutableVector newImage = newImages.GetCol(i);
        Vector oldImage = _images.GetCol(i);
        for (int r = 0; r < newNumRows; ++r) {
          for (int c = 0; c < newNumCols; ++c) {
            double newVal = 0;
            for (int oR = r * reductionFactor; oR < (r+1) * reductionFactor; ++oR) {
              for (int oC = c * reductionFactor; oC < (c+1) * reductionFactor; ++oC) {
                int idx = oR * _numRows + oC;
                newVal += oldImage[idx];
              }
            }
            newVal /= (reductionFactor * reductionFactor);
            int newIdx = r * newNumRows + c;
            newImage[newIdx] = newVal;
          }
        }
      }
      _numRows = newNumRows;
      _numCols = newNumCols;
      _numPixels = newNumPixels;
      _images.Swap(newImages);
    }

    ifstream labStream(labels, ios::in|ios::binary);
    if (!labStream.is_open()) {
      cout << "Couldn't open feature file " << labels.c_str() << endl;
      abort();
    }

    labStream.read((char*)header, 8);
    assert (SwapEndian(header[1]) == numImages);

    _labels.resize(numImages);
    labStream.read((char*)&_labels[0], numImages);

    // permute training instances
    AllocatingVector temp(_numPixels);
    for (int i = 1; i < numImages; ++i) {
      int j = (int)(rnd.drand48() * (i + 1));
      MutableVector imageI = _images.GetCol(i), imageJ = _images.GetCol(j);
      temp.CopyFrom(imageI);
      imageI.CopyFrom(imageJ);
      imageJ.CopyFrom(temp);
      swap(_labels[i], _labels[j]);
    }

    _numFeatures = -1;
  }

  int NumFeatures() const {
    return (_numFeatures == -1) ? _numPixels : _numFeatures;
  }

  int NumLabels() const { return 10; }

  int NumImages() const { return _labels.size(); }

  Vector GetImage(int i) const { return _images.GetCol(i); }

  Matrix GetImages(int start = 0, int end = -1) const { return _images.GetCols(start, end); }

  vector<unsigned char> GetLabels(int start = 0, int end = -1) const {
    if (end == -1) end = _labels.size();
    return vector<unsigned char>(_labels.begin() + start, _labels.begin() + end);
  }

  static void Print(Vector digit, int rows, int cols) {
    char chars[] = {' ', '.', '+', '*', '#'};
    for (int r = 0; r < rows; ++r) {
      for (int c = 0; c < cols; ++c) {
        double x = digit[r * cols + c];
        int val = (int)(5 * x);
        if (val < 0) val = 0; if (val >= 5) val = 4;
        cout << chars[val];
      }
      cout << endl;
    }
  }

  static void PrintSideBySide(Vector imageL, Vector imageR, int rows, int cols) {
    char chars[] = {' ', '.', '+', '*', '#'};
    for (int r = 0; r < rows; ++r) {
      for (int c = 0; c < cols; ++c) {
        double x = imageL[r * cols + c];
        int val = (int)(5 * x);
        if (val < 0) val = 0; if (val >= 5) val = 4;
        cout << chars[val];
      }
      cout << '|';
      for (int c = 0; c < cols; ++c) {
        double x = imageR[r * cols + c];
        int val = (int)(5 * x);
        if (val < 0) val = 0; if (val >= 5) val = 4;
        cout << chars[val];
      }
      cout << endl;
    }
  }

  void Print(int digitIdx) {
    Vector digit = _images.GetCol(digitIdx);

    Print(digit, _numRows, _numCols);
  }
};
