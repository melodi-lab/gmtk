
#include "rand.h"
#include "range.h"

RAND rnd(true);

int
main(int argc, char *argv[]) {
  Range r(argv[1],0,1000);

  printf("'%s' length %u\n", argv[1], r.length());
  for (Range::iterator it=r.begin(); !it.at_end(); ++it) {
    printf("%3d ", *it);
  }
  printf("\n");
  for (Range::permuter p=r.permute(); !p.at_end(); ++p) {
    printf("%3d ", *p);
  }
  printf("\n");
  return 0;
}
