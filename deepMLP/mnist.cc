
/*
 * Written by Richard Rogers <rprogers@ee.washington.edu>
 *
 * Copyright (C) 2013 Jeff Bilmes
 * Licensed under the Open Software License version 3.0
 * See COPYING or http://opensource.org/licenses/OSL-3.0
 *
 */

#include "Matrix.h"
#include "MNIST.h"
#include "rand.h"

RAND rnd;

int
main(int argc, char *argv[]) {
  if (argc != 5) {
    fprintf(stderr, "%s trainImages trainLabels testImages testLabels\n\n", argv[0]);
    exit(1);
  }
  string trainImageFile(argv[1]);
  string trainLabelFile(argv[2]);
  string testImageFile(argv[3]);
  string testLabelFile(argv[4]);

  MNISTData mnist(trainImageFile, trainLabelFile);

  int inputSize = mnist.NumFeatures(), outputSize = 10, 
      numInstances = mnist.NumImages();

  AllocatingMatrix output(outputSize, numInstances);

  vector<unsigned char> labels = mnist.GetLabels();
  for (int i = 0; i < numInstances; ++i) {
    output.AtMutable(labels[i], i) = 1;
  }

  for (int i=0; i < numInstances; i+=1) {
    printf("0 %d", i);
    Vector image = mnist.GetImage(i);
    for (int j=0; j < image.Len(); j+=1) {
      printf(" %f", image[j]);
    }
    Vector label = output.GetCol(i);
    for (int j=0; j < label.Len(); j+=1) {
      printf(" %f", label[j]);
    }
    printf("\n");
  }
  return 0;
}
