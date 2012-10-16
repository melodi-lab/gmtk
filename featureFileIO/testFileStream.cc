
#include <stdlib.h>
#include <stdio.h>

#include "GMTK_FileStream.h"
#include "GMTK_PFileFile.h"
#include "GMTK_FileSource.h"


char    *ofs[1];
unsigned nfs[1];
unsigned nis[1];
const char   *fmts[1];
const char    *frs[1];
const char    *irs[1];
const char    *prefrs[1];
const char    *preirs[1];
const char    *sr[1];
char  *prepr[1];
char *postpr[1];
char *gpr_str;
unsigned justification;
bool iswp[1];
unsigned ifmts[1];
bool Cpp_If_Ascii;
char *cppCommandOptions;
char    *Per_Stream_Transforms[1];
char    *Post_Transforms;
int startSkip;
int endSkip;
unsigned    Action_If_Diff_Num_Frames[1];
unsigned    Action_If_Diff_Num_Sents[1];
unsigned Ftr_Combo;

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
  FileSource fsrc(pf);
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

