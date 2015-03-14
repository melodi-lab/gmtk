# Net Surgery as done in http://nbviewer.ipython.org/github/BVLC/caffe/blob/master/examples/net_surgery.ipynb

import numpy as np
import matplotlib.pyplot as plt


import sys
import caffe

# Set the right path to your model definition file, pretrained model weights

MODEL_FILE_DEPLOY = 'caffe/network_definition_deploy.prototxt'
PRETRAINED_NN = 'caffe/snapshot_adagrad/network_definition_adagrad_iter_2200000.caffemodel'


INPUT_LAYER_SIZE = 351 # Num of nodes in the 1st layer of NN


GMTK_WEIGHT_FILE = 'learnedParams/learned_dmlp'


FLAT_LENGTH = 1000

MAX_COUNT = 10000


# For 10 layer networks
params = ['fc1', 'fc2', 'fc3', 'fc4', 'fc5', 'fc6','fc7','fc8','fc9','fc10']


numMatrices = len(params)

net_dim_red = caffe.Classifier(MODEL_FILE_DEPLOY, PRETRAINED_NN, image_dims=(INPUT_LAYER_SIZE, 1))

fc_params = {pr: (net_dim_red.params[pr][0].data, net_dim_red.params[pr][1].data) for pr in params}

fdGMTKwt = open(GMTK_WEIGHT_FILE, "w")

fdGMTKwt.write(str(numMatrices)+"\n")

for fc in params:
    print '{} weights are {} dimensional and biases are {} dimensional'.format(fc, fc_params[fc][0].shape, fc_params[fc][1].shape)


for matId in range(0,numMatrices):
	fdGMTKwt.write(str(matId)+"\n")
	fdGMTKwt.write("g"+str(matId)+ " " + str(fc_params[params[matId]][0].shape[2]) + " " + str(fc_params[params[matId]][0].shape[3] + 1) + " ")
	
	for row in range(0,fc_params[params[matId]][0].shape[2]):
		for column in range (0,fc_params[params[matId]][0].shape[3]):
			fdGMTKwt.write(str(fc_params[params[matId]][0][0,0,row,column])+ " ")
		fdGMTKwt.write(str(fc_params[params[matId]][1][0,0,0,row])+ "\n")


	


