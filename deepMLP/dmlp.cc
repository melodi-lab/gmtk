#include "DBN.h"
#include "MNIST.h"

int
main(int argc, char *argv[]) {

  //   MNISTData mnist("C:\\svn\\data\\mnist\\train-images.idx3-ubyte", "C:\\svn\\data\\mnist\\train-labels.idx1-ubyte");
  //   MNISTData testMnist("C:\\svn\\data\\mnist\\t10k-images.idx3-ubyte", "C:\\svn\\data\\mnist\\t10k-labels.idx1-ubyte");
  if (argc != 5) {
    fprintf(stderr, "%s trainImages trainLabels testImages testLabels\n\n", argv[0]);
    exit(1);
  }
  string trainImageFile(argv[1]);
  string trainLabelFile(argv[2]);
  string testImageFile(argv[3]);
  string testLabelFile(argv[4]);

  MNISTData mnist(trainImageFile, trainLabelFile);

  int inputSize = mnist.NumFeatures(), outputSize = 10, hiddenSize = 200, numInstances = mnist.NumImages();
  int layers = 5;

  DBN dbn(layers, inputSize, hiddenSize, outputSize, Layer::ActFunc(Layer::ActFunc::TANH), Layer::ActFunc(Layer::ActFunc::CUBIC), Layer::ActFunc(Layer::ActFunc::LINEAR));

  AllocatingMatrix output(outputSize, numInstances);
  DBN::ObjectiveType objType = DBN::SOFT_MAX;

  vector<unsigned char> labels = mnist.GetLabels();
  for (int i = 0; i < numInstances; ++i) {
    output.At(labels[i], i) = 1;
  }

  vector<DBN::HyperParams> pretrainHyperParams(layers);
  for (int i = 0; i < layers; ++i) {
    pretrainHyperParams[i].pretrainType = DBN::AE;
    pretrainHyperParams[i].numUpdates = 25000;
    pretrainHyperParams[i].numAnnealUpdates = 10000;
  }

  DBN::HyperParams bpHyperParams;
  bpHyperParams.numUpdates = 200000;
  bpHyperParams.numAnnealUpdates = 20000;

  Random rand(2);
  dbn.Train(mnist.GetImages(), output, objType, rand, false, pretrainHyperParams, bpHyperParams);

  // really we shouldn't be using the test set like this, but this is just for debugging
  MNISTData testMnist(testImageFile, testLabelFile);
  Matrix testImages = testMnist.GetImages();
  vector<unsigned char> testLabels = testMnist.GetLabels();
  Matrix testScores = dbn.MapUp(testImages);

  int correct = 0;
  for (int i = 0; i < testScores.NumC(); ++i) {
    int answer = ArgMax(testScores.GetCol(i));
    if ((unsigned char)answer == testLabels[i]) correct++;
  }
  cout << "Accuracy: " << ((float)correct / testScores.NumC()) << endl;
  return 0;
}
