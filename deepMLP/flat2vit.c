
#include <stdlib.h>
#include <stdio.h>

/* flat2vit var1 ... */

int
main(int argc, char *argv[]) {
  int n = argc - 1;
  unsigned *vals = (unsigned *) malloc(n * sizeof(unsigned));
  unsigned seg, frame;

  while (fscanf(stdin, "%u %u", &seg, &frame) == 2) {
    int i;
    for (i=0; i < n; i+=1)
      fscanf(stdin, "%u", vals+i);
    for (i=0; i < n; i+=1)
      printf("seg %u: %s(%u)=%u\n", seg, argv[i+1], frame, vals[i]);
  }
  exit(0);
}
