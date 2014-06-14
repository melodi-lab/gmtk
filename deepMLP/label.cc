
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

  int //inputSize = mnist.NumFeatures(), outputSize = 10,  -- unused
      numInstances = mnist.NumImages();

  vector<unsigned char> labels = mnist.GetLabels();

  for (int i=0; i < numInstances; i+=1) {
    printf("0 %d %d\n", i, (int)labels[i]);
  }
  return 0;
}
