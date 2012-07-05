#include "readRange.h"
#include "error.h"

#define FILE_NAME_LEN 25

//reads lines from the file fileName and creates sets A and B for each.
void RangeSetCollection::readSets(const char* fileName) {
  
  char line[MAX_LINE_LEN];
  FILE *ifp;

  if((ifp=fopen(fileName,"r")) == NULL) {
    printf("Couldn't open range input file for reading\n");
    exit(-1);
  }

  int lineno=0;
  while(fgets(line,MAX_LINE_LEN,ifp) != NULL) {
    RangeSet X;
    lineno++;
    int status = X.createSets(line);
    if(status == COMMENT) {
      //ignore empty and commented lines
      //printf("Comment line\n");
    }
    else if(status != SUCCESS) {
      error("Invalid line number %d (error code = %d)\n",lineno,status);
    }
    else {
      addToCollection(X);
#if 0
	printf("A: ");
	for(int i=0;i<X.getSize(0);++i)     printf("(%d, %d) ", X.set[0][i].feat, X.set[0][i].lag);
	printf("\nB: ");
	for(int i=0;i<X.getSize(1);++i)     printf("(%d, %d) ", X.set[1][i].feat, X.set[1][i].lag);
	printf("\n");
#endif
    }
  }
}


