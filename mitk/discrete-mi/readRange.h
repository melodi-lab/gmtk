#ifndef READ_RANGE_H
#define READ_RANGE_H

#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include <assert.h>
#include <ctype.h>

#include<vector>

#define MAX_SET_SIZE 100
#define MAX_LINE_LEN 200

//---------A bunch of error definitions
#define COMMENT -1
#define SUCCESS 0
#define BRACKET_ERROR 1
#define UNRECOGNIZED_TOKEN 2
#define INVALID_FORMAT 3
#define EMPTY_SET 4
//------------------------------------

using namespace std;

typedef struct _feat_lag_pair {
    int feat, lag;
} feat_lag_pair;

class RangeSet {
    int size[2];

 public:
    RangeSet() {
	size[0] = size[1] = 0;
	min_lag = max_lag = 0;
    }

    int getSize(int setNum) {
	return size[setNum];
    }

    int createSets(char* line);
    int createOneSet(int setNum, char* &lp);

    feat_lag_pair set[2][MAX_SET_SIZE];
    int min_lag, max_lag;
};

class RangeSetCollection {

  int numSets;
  
public:
  std::vector<RangeSet> rs;

  RangeSetCollection() { numSets = 0; }
  int getSize() { return rs.size(); }
  void addToCollection(RangeSet r) {
    rs.push_back(r);
  }
  void readSets(const char* fileName);
  RangeSetCollection(const char* fileName) {
    readSets(fileName);
  }

};

#endif
