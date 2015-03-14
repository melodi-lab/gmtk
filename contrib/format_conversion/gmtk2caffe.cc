
/*
 Build with

	$(CXX) gmtk2caffe.cc -o gmtk2caffe -lhdf5_cpp -lhdf5

 The correct command should be something like

	h5c++ gmtk2caffe.cc -o gmtk2caffe

 where h5c++ is a wrapper around a C++ compiler that knows
 how to link with the desired HDF5 libarary installation.
 Unfortunately, the h5c++ on the MELODI machines seems to
 require staticly linked HDF5 libraries but only dynamically
 linked libraries appear to be installed.
*/


#include <stdio.h>
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
  using namespace H5;
#endif

int
main(int argc, char *argv[]) {

  if (argc != 3) {
    fprintf(stderr, "\ngmtk2caffe input output\n\n"
"input must be an HDF5 file containing an N x M 'continuous' dataset where\n"
"N is the number of M x 1 feature vectors and an N x 1 'discrete' dataset \n"
"containing the labels of the feature vectors. output will be an HDF5 file\n"
"containing an N x 1 dataset named 'label' (of type NATIVE_FLOAT) and an\n"
"N x 1 x 1 x M dataset named 'data' containing the feature vectors.\n\n");
    return 1;
  }
  const H5std_string inputFilename(argv[1]);
  const H5std_string discName("discrete");
  const H5std_string contName("continuous");

  H5File file(inputFilename, H5F_ACC_RDONLY);

  DataSet discDataset = file.openDataSet(discName);
  DataSet contDataset = file.openDataSet(contName);

  DataSpace discSpace = discDataset.getSpace();
  DataSpace contSpace = contDataset.getSpace();

  hsize_t discDims[2];
  hsize_t contDims[2];
  int ndims = discSpace.getSimpleExtentDims(discDims);
  ndims = contSpace.getSimpleExtentDims(contDims);


  const H5std_string outputFilename(argv[2]);
  H5File outputFile(outputFilename, H5F_ACC_TRUNC);

  hsize_t labelDims[2] = {discDims[0],1};
  DataSpace labelDataspace(2, labelDims);
  DataSet labelDataset = 
    outputFile.createDataSet("label", PredType::NATIVE_FLOAT, labelDataspace);

  hsize_t dataDims[4] = {contDims[0], 1, 1, contDims[1]};
  DataSpace dataDataspace(4, dataDims);
  DataSet dataDataset =
    outputFile.createDataSet("data", PredType::NATIVE_FLOAT, dataDataspace);
  
  hsize_t discStart[2]={0,0};
  hsize_t discLength[2]={1,1};
  hsize_t discMemdim[1] = {1};
  DataSpace discMemspace(1, discMemdim);

  hsize_t contStart[2]={0,0};
  hsize_t contLength[2]={1,contDims[1]};
  hsize_t contMemdim[2] = {1,contDims[1]};
  DataSpace contMemspace(2, contMemdim);

  hsize_t dataStart[4]={0,0,0,0};
  hsize_t dataLength[4]={1, 1, 1, contDims[1]};
  hsize_t dataMemdim[4] = {1, 1, 1, contDims[1]};
  DataSpace dataMemspace(4, dataMemdim);

  float featureVector[contDims[1]];
  float featureTensor[1][1][1][contDims[1]];

  for (unsigned i=0; i < discDims[0]; i+=1) {
    discStart[0] = i;
    discSpace.selectHyperslab(H5S_SELECT_SET, discLength, discStart);
    unsigned l; float lf;
    discDataset.read(&l, PredType::NATIVE_UINT32, discMemspace, discSpace);

    labelDataspace.selectHyperslab(H5S_SELECT_SET, discLength, discStart);
    lf = (float)l;
    labelDataset.write(&lf,PredType::NATIVE_FLOAT, discMemspace,labelDataspace);

    contStart[0]=i; contStart[1]=0;
    contSpace.selectHyperslab(H5S_SELECT_SET, contLength, contStart);
    contDataset.read(featureVector, PredType::NATIVE_FLOAT, contMemspace, contSpace);

    dataStart[0]=i;
    dataDataspace.selectHyperslab(H5S_SELECT_SET, dataLength, dataStart);
    for (unsigned j=0; j < contDims[1]; j+=1) {
      featureTensor[0][0][0][j] = featureVector[j];
    }
    dataDataset.write(featureTensor, PredType::NATIVE_FLOAT, dataMemspace, dataDataspace);
  }

  return 0;
}
