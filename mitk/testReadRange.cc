#include <stdio.h>
#include <stdlib.h>
#include "readRange.h"

#define FILE_NAME_LEN 25
#define DEFAULT_FILE_NAME "rangeFile.dat"

int main(int argc, char* argv[]) {
  
  char fileName[FILE_NAME_LEN];

  if(argc==2) strcpy(fileName, argv[1]);
  else strcpy(fileName, DEFAULT_FILE_NAME);

  RangeSetCollection* rsc = new RangeSetCollection(fileName);; 
  RangeSet X;
  for(int i=0; i<rsc->getSize(); i++) {
    X = rsc->rs[i];
    printf("A: ");
    for(int i=0;i<X.getSize(0);++i)     printf("(%d, %d) ", X.set[0][i].feat, X.set[0][i].lag);
    printf("\nB: ");
    for(int i=0;i<X.getSize(1);++i)     printf("(%d, %d) ", X.set[1][i].feat, X.set[1][i].lag);
    printf("\n");
  }

}
