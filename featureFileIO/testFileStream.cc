
#include <stdlib.h>
#include <stdio.h>

#include "GMTK_FileStream.h"
#include "GMTK_PFileFile.h"
#include "GMTK_FileSource.h"

// testFileStream file nf ni
int
main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "testFileStream pfile nf ni\n");
    exit (1);
  }

  unsigned nf = atoi(argv[2]);
  unsigned ni = atoi(argv[3]);

  PFileFile *pf = new PFileFile(argv[1], nf, ni, 0, true);
  ObservationFile *ofs[1] = {pf};
  FileSource fsrc(1, ofs);
  FileStream fs(pf);

  unsigned seg = 0;
  unsigned frm = 0;

  Data32 const *frame;
  for (; !fs.EOS(); ) {
    frame = fs.getNextFrame();
    if (!frame) {
      printf("eos\n");
      seg += 1;
      frm = 0;
      continue;
    }
    printf("%03u %03u", seg, frm);
    for (unsigned i=0; i < nf; i+=1) {
      printf(" %f", ((float *)frame)[i]);
    }
    for (unsigned i=0; i < ni; i+=1) {
      printf(" %d", ((int *)frame)[fs.numContinuous() + i]);
    }
    printf("\n");
    frm += 1;
  }
  exit(0);
}

