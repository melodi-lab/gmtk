#include <stdlib.h>
#include <stdio.h>

int
main(int argc, char *argv[]) {
  int i;
  for (i=0; i < 20; i+=1) {
    printf("stdout %02d\n", i);
    fprintf(stderr,"stderr %02d\n", i);
  }
  exit(0);
}
