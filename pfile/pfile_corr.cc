#include <stdio.h>
#include <stddef.h>
#include <iostream.h>
#include <string.h>
#include <math.h>

#include "pfile.h"
#include "Range.H"
#include "error.h"

#define PROB_LOW 0.00001
#define PROB_MIN 0.0001

#define USAGESTRING "USAGE: pfile_confusion [options]\n
Option                                     default\n
-i1  input pfile1                            none\n
-i2  input pfile2                            none\n
-iswap1 do byteswapping on file 1            none\n
-iswap2 do byteswapping on file 2            none\n
-o   output confusion table                  none\n
-e   compute correlation between errors only off\n
-c   compute correlation between hits only   off\n"

float **
CreateMatrix(int size) {
  
  int i;
  float **m;

  
  m = new float*[size];

  for (i=0;i<size;i++) 
    m[i] = new float[size];
  return m;
}

void
ZeroMatrix(float ** m, int size) {
  
  int i,j;
  for (i=0;i<size;i++)
    for (j=0;j<size;j++)
      m[i][j] = 0;
  return;
}


FILE* 
iopen(const char *fname)
{
  FILE *file;
  char ch;

  if ((file=fopen(fname,"r")) != NULL) 
    {
      fprintf(stderr,"File %s exists, overwrite? (y/n) ",fname);
      if ((ch=getchar()) == 'y') 
	{
	  remove(fname);
	  file=fopen(fname,"wb");
	}
      else
	exit(03);
      rewind(stdin);
    }
  else 
    file=fopen(fname,"wb");
  return file;
}

int
GetMax(float *vec, int size) {

  float m = 0.0;
  int mInd = -1;

  for (int i=0;i<size;i++) {
    if (vec[i] > m) {
      m = vec[i];
      mInd = i;
    }
  }
  return mInd;
}

float
Entropy(float *vec, int size) {
  float out = 0.0;
  int i;

  for (i=0;i<size;i++) {
    if (vec[i] <= PROB_LOW)
      continue;
    else
      out += (vec[i] * log(vec[i]));
  }
  return out;
}


void 
pfile_corr(InFtrLabStream_PFile *data1,InFtrLabStream_PFile *data2,int err, int corr)
{

  const size_t n_sents1 = data1->num_segs();
  const size_t n_sents2 = data2->num_segs();
  const size_t n_ftrs1 =  data1->num_ftrs();
  const size_t n_ftrs2 =  data2->num_ftrs();
  const size_t n_labs1 = data1->num_labs();
  const size_t n_labs2 = data2->num_labs();
  size_t same_error = 0, diff_error = 0;
  size_t total_frames;
  size_t n_errors= 0;

  fprintf(stderr,"File 1: %li sents, %li ftrs, %li labs\n",
	  n_sents1,n_ftrs1,n_labs1);
  fprintf(stderr,"File 2: %li sents, %li ftrs, %li labs\n",
	  n_sents2,n_ftrs2,n_labs2);
  
  size_t n_frames = 0;
  float term1,term2,rho;
  float sum_xy = 0, sum_x2 = 0, sum_y2 = 0;
  float sum_x = 0, sum_y = 0;

  if (n_sents1 != n_sents2) {
    fprintf(stderr,"Number of sentences is different: %li vs. %li\n",
	    n_sents1,
	    n_sents2);
    exit(01);
  }
  if (n_ftrs1 != n_ftrs2) {
    fprintf(stderr,"Number of sentences is different: %li vs. %li\n",
	    n_ftrs1,
	    n_ftrs2);
    exit(01);
  }
  size_t buf_size = 300;
  float *ftr_buf1 = new float[n_ftrs1 * buf_size];
  float *ftr_buf2 = new float[n_ftrs2 * buf_size];
  UInt32* lab_buf1 = new UInt32[n_labs1 * buf_size];
  UInt32* lab_buf2 = new UInt32[n_labs2 * buf_size];

  for (int s=0; s <n_sents1;s++) {
    int pos = 0; 
    int lpos = 0;

    if ( (s % 100) == 0)
    printf("Processing sentence %i\n",s);

    int seg_frames1 = data1->num_frames(s);
    int seg_frames2 = data2->num_frames(s);
      
    if (seg_frames1 == SIZET_BAD) 
      {
	fprintf(stderr,
		"Couldn't find number of frames "
		"at sentence %lu in input pfile 1.\n",
		(unsigned long) s);
	error("Aborting.");
      }
    if (seg_frames2 == SIZET_BAD) 
      {
	fprintf(stderr,
		"Couldn't find number of frames "
		"at sentence %lu in input pfile 1.\n",
		(unsigned long) s);
	error("Aborting.");
      }
    const SegID seg_id1 = data1->set_pos(s,0);
    const SegID seg_id2 = data2->set_pos(s,0);
      
      if (seg_id1 == SEGID_BAD) {
	fprintf(stderr,
		"Couldn't seek to start of sentence %lu "
		"in input pfile 1.",
		(unsigned long) s);
	error("Aborting.");
      }
      if (seg_id2 == SEGID_BAD) {
	fprintf(stderr,
		"Couldn't seek to start of sentence %lu "
		"in input pfile 2.",
		(unsigned long) s);
	error("Aborting.");
      }
      if (seg_frames1 > buf_size) {
	delete ftr_buf1;
	delete ftr_buf2;
	delete lab_buf1;
	delete lab_buf2;
	buf_size = seg_frames1 * 2;
	ftr_buf1 = new float[buf_size * n_ftrs1];
	ftr_buf2 = new float[buf_size * n_ftrs2];
	lab_buf1 = new UInt32[buf_size * n_labs1];
	lab_buf2 = new UInt32[buf_size * n_labs2];
      }
      
      const size_t n_read1 = data1->read_ftrslabs(seg_frames1,ftr_buf1,lab_buf1);
      
      const size_t n_read2 = data2->read_ftrslabs(seg_frames2,ftr_buf2,lab_buf2);
      if ( n_read1 != seg_frames1) {
	fprintf(stderr, " At sentence %lu in input pfile 1, "
		"only read %lu frames when should have read %lu.\n",
		(unsigned long) s,
		(unsigned long) n_read1, (unsigned long) seg_frames1);
	error("Aborting.");
      }
      if ( n_read2 != seg_frames2) {
	fprintf(stderr, " At sentence %lu in input pfile 2, "
		"only read %lu frames when should have read %lu.\n",
		(unsigned long) s,
		(unsigned long) n_read2, (unsigned long) seg_frames2);
	error("Aborting.");
      }
      if (n_read1 != n_read2)
	error("Read different number of frames at sent %li\n",s);
      
      float *fp1 = ftr_buf1;
      float *fp2 = ftr_buf2;
      UInt32 *lp1 = lab_buf1;
      UInt32 *lp2 = lab_buf2;

      
      for (int f = 0; f < n_read1; f++) {
	if (*lp1 != *lp2)
	  error("Different labels in pfiles at sent %li frame %li!\n",s,f);

	// get maximum output unit
	int m = GetMax(fp1,n_ftrs1);
	int n = GetMax(fp2,n_ftrs2);
	if (m == -1 || n == -1)
	  error("Max output not found: file1: %i file2: %i at sent %li, frame %li\n",m,n,
		s,f);

	if (m != *lp1  && n != *lp1) {
	  n_errors++;
	  if (m == n)
	    same_error++;
	  else 
	    diff_error++;

	// only correlation among correct outputs

	  if (corr) {
	    fp1 += n_ftrs1;
	    fp2 += n_ftrs2;
	    lp1++; lp2++;
	    continue;
	  }

	}
	if (m == *lp1 || n == *lp1) {	
	  // only correlation among errors
	  if (err) {
	    fp1 += n_ftrs1;
            fp2 += n_ftrs2;
            lp1++; lp2++;
	    continue;
	  }
	}
	sum_xy += m*n;
	sum_x += m;
	sum_y += n;
	sum_x2 += m*m;
	sum_y2 += n*n;
	n_frames++;

	fp1 += n_ftrs1;
	fp2 += n_ftrs2;
	lp1++;
	lp2++;
      }
  }
  term1 = sum_xy - (sum_x * sum_y)/n_frames;
  term2 = sqrt((sum_x2 - (sum_x*sum_x)/n_frames) * (sum_y2 - (sum_y*sum_y)/n_frames));
  rho = term1/term2;
	       
  size_t total_error = same_error + diff_error;
  printf("%li frames processed\nCorrelation coefficient: %f\n",n_frames,rho);
  printf("Percentage of different errors: %.2f\n",diff_error/(float)n_errors * 100.0f);

  delete ftr_buf1;
  delete ftr_buf2;
  delete lab_buf1;
  delete lab_buf2;
  return;
}


void 
main(int argc, const char **argv) {

  const char *inname1=NULL,*inname2=NULL;      
  const char *program_name = *argv++;
  int debug_level = 0;
  bool iswap1 = false, iswap2 = false;
  int errors_only = 0, correct_only = 0;
  argc--;
  
  if (argc == 0)
    error(USAGESTRING);
  
  while (argc--)
    {
      const char *argp = *argv++;
      
      if (strcmp(argp, "-help")==0)
	{
	  printf("%s\n",USAGESTRING);
	}
      else if (strcmp(argp, "-i1")==0)
	{
	  if (argc>0)
	    {
	      inname1 = *argv++;
	      argc--;
	    }
	  else
	    error("No input pfile name specified\n");
	}
      else if (strcmp(argp,"-i2")==0)
	{
	  if (argc>0)
	    {
	      inname2 = *argv++;
	      argc--;
	    }
	  else
	    error("No input pfile name specified\n");
	}
      else if (strcmp(argp,"-e") == 0)
	{
	      errors_only = 1;
	}
      else if (strcmp(argp,"-c") == 0)
	{
	  correct_only = 1;
	}
       else if (strcmp(argp,"-iswap1") == 0)
         {
           iswap1 = true;
         }
      else if (strcmp(argp,"-iswap2") == 0)
         {
           iswap2 = true;
         }
      else {
	error("Unrecognized argument (%s).",argp);
      }
    }
  if (inname1 == 0)
    fprintf(stderr,"No input pfile name supplied.");
  if (inname2 == 0)
    fprintf(stderr,"No input pfile name supplied.");

  FILE *infile1 = fopen(inname1, "r");
  if (infile1 == NULL)
    error("Couldn't open input pfile for input\n");

  FILE *infile2 = fopen(inname2, "r");
  if (infile2 == NULL)
    error("Couldn't open input pfile for input\n");

  InFtrLabStream_PFile* input1
    = new InFtrLabStream_PFile(debug_level, "", infile1, 1,iswap1);
  
  InFtrLabStream_PFile* input2
    = new InFtrLabStream_PFile(debug_level, "", infile2, 1,iswap2);

  if (errors_only == 1 && correct_only == 1)
    error("Contradictory flags: -c and -e\n");
  if (errors_only)
    printf("Computing correlation among errors only\n");
  if (correct_only)
    printf("Computing correlation among correct outputs only\n");
  pfile_corr(input1,input2,errors_only,correct_only);
  
  delete input1;
  delete input2;
  fclose(infile1);
  fclose(infile2);
  return 0;
}
    





