#include <stdio.h>
#include <stddef.h>
#include <iostream.h>
#include <math.h>
#include <stdlib.h>
#include "pfile.h"
#include "error.h"


#define USAGESTRING "USAGE: pfile_combine [options]\n
Option                                 default\n
-i1   input pfile                       none\n
[-i2   second input pfile               none]\n
-o    output file                       none\n
-f    format for output file            none\n
-n    normalize                         off\n
-sel  selection                         off\n
-s    use sum rule                      off\n
-p   use product rule                   off\n
-max    use max rule                    off\n
-min    use min rule                    off\n
-smin    use smin rule                  off\n
-psmin  use psmin rule                  off\n
-esmin  use esmin rule                  off\n
-pesmin use pesmin rule                 off\n
-b    beta constant for smin rules      off\n
-iswap1 do byteswapping on file 1       off\n
-iswap2 do byteswapping on file 2       off\n
-oswap  do byteswapping on output       off\n
-v    verbose output
-debug  debug level                     0\n
"

#define FLOAT_MIN 1.0E-10
#define LOGMINARG 2.45E-308
#define EXPMINARG -708.3

// combination rules

typedef enum {
  NONE,
  MIN,
  MAX,
  PROD,
  SUM,
  SMIN,
  ESMIN,
  PSMIN,
  PESMIN,
  SELECT,
};
  

size_t num_zeros = 0;

int verbose = 0;

FILE * iopen(const char *fname)
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

float neglog(float f) {

  if (f < LOGMINARG) {
    num_zeros++;
    return FLOAT_MIN;
  }
  else
    return -log(f);
}

float mylog(float f) {

  if (f < LOGMINARG) {
    num_zeros++;
    return FLOAT_MIN;
  }
  else
    return log(f);
}


double myexp(double d) {

  if (d < EXPMINARG)
    return FLOAT_MIN;
  else
    return exp(d);
}

// hard min and max

double
minf(float f1,float f2) {

  return (double)(f1 > f2) ? f2 : f1;
}

double
maxf(float f1,float f2) {
  
  return (double) (f1 > f2) ? f1: f2;
}
  


// this is smin

double
sminf(float f1, float f2, float beta) {

  double d1,d2,tmp1,tmp2,b;
 
  d1 = (double)f1;
  d2 = (double)f2;
  b = -(double)beta;
  
  tmp1 = pow(d1,b) + pow(d2,b);
  tmp2 = pow(tmp1,1/b);
  if (verbose)
    printf("in: %e %e, tmp: %e %e, out: %f\n",
	   f1,f2,tmp1,tmp2,tmp2);
  return tmp2;
}

// this is psmin

double
psminf(float f1,float f2,float beta) {

  double d1,d2,tmp1,tmp2,tmp3,b;
  
  if ((f1 == 1.0 || f2 == 1.0) && beta < 1)
     return 1.0;
  if (f1 == 0.0 && f2 == 0.0)
    return 0.0;
    d1 = 1/(double)f1;
    d2 = 1/(double)f2;
    b = (double)beta;
    if (d1 == 0.0)
      tmp1 = pow((double)log(d2),beta);
    else if (d2 == 0.0)
      tmp1 =  pow((double)log(d1),beta);
    else
      tmp1 = pow((double)log(d1),beta) + pow((double)log(d2),beta);
    tmp2 = pow(tmp1,1/beta);
    tmp3 =  exp(-tmp2);
    if (verbose) {
   //      printf("psmin: in: %e %e, d: %e %e\nlog: %e %e\n%e %e\n tmp: %e %e, out: %f\n",
      //	     f1,f2,d1,d2,log(d1),log(d2),pow(log(d1),beta),pow(log(d2),beta),tmp1,tmp2,tmp3);
    }
  return tmp3;
}
  
// this is esmin


double
esminf(float f1, float f2, float beta) {

  double d1,d2,tmp1,tmp2,b,sum;
  d1 = (double)f1;
  d2 = (double)f2;
  b = -(double)beta;
  
  tmp1 = exp(b * d1);
  tmp2 = exp(b * d2);
  sum = tmp1+tmp2;
  if (verbose) {
    printf("%f %f %f\n",tmp1,tmp2,sum);
    printf("out: %f\n",((d1 *tmp1/sum) + (d2 *  tmp2/sum)));
  }
  return ((d1 *tmp1/sum) + (d2 *  tmp2/sum));
}


// and this is pesmin
double
pesminf(float f1, float f2, float beta) {

  double d1,d2,tmp1,tmp2,b,f_sum,s_sum;

  if (f1 == 0.0 && f2 == 0.0)
    return 0.0;

  b = (double)beta; 
  d1 = (1/(double)f1);
  d2 = (1/(double)f2);


  tmp1 = pow(d1,b);
  tmp2 = pow(d2,b);

  f_sum = tmp1 + tmp2;
  if (f1 == 0.0)
    s_sum = ((double)log(f2) * tmp2/f_sum);
  else if (f2 == 0.0)
    s_sum = ((double)log(f1) * tmp1/f_sum);
  else
    s_sum = ((double)log(f1) * tmp1/f_sum) + ((double)log(f2) * tmp2/f_sum);
  if (verbose) {
    //    printf("%e %e %e %e %e %e %e\n",f1,f2,tmp1,tmp2,f_sum,s_sum,exp(s_sum));
  }

  return exp(s_sum);
}


double
fselect(float f1, float f2, float end1, float end2) {

  if (end1 < end2)
    return (double) f1;
  else
    return (double) f2;
}
int
GetMax(float *vec, int size) {

  float *fp;
  float m;
  int mInd = -1;
  int i;

  for (i=0,fp = vec;i < size;i++,fp++) {
    if (i == 0) {
      m = *fp;
      mInd = 0;
    }
    else if (*fp > m) {
      m = *fp;
      mInd = i;
    }
  }
  return mInd;
}


void
combine(float *f1,float *f2,float *out,int size,int norm,int rule,float beta)

{

  float *fp1,*fp2,*outp;
  double sum=0.0;
  double *tmp,*tp;
  int i,b,c;
  tmp = new double[size];
  tp = tmp;
  fp1 = f1;
  fp2 = f2;
  outp = out;
  float *end1;
  float *end2;

  switch(rule) {
    // averaging
  case SUM:
    for (i = 0; i < size; i++) {
      *tp = (*fp1 + *fp2)/2;
      sum += *tp;
      tp++;
      fp1++;
      fp2++;
    }
    break;
    // product rule
  case PROD:
    for (i = 0; i < size; i++) {
      *tp = (*fp1 * *fp2);
      sum += *tp;
      tp++;
      fp1++;
      fp2++;
    }
    break;
    // hard min rule
  case MIN:
    for (i = 0; i < size; i++) {
      *tp = minf(*fp1,*fp2);
      sum += *tp;
      tp++;
      fp1++;
      fp2++;
    }
    break;
  case MAX:
    for (i = 0; i < size; i++) {
      *tp = maxf(*fp1,*fp2);
      sum += *tp;
      tp++;
      fp1++;
      fp2++;
    }
    break;
  case SMIN:
    for (i = 0; i < size; i++) {
      *tp = sminf(*fp1,*fp2,beta);
      sum += *tp;
      tp++;
      fp1++;
      fp2++;
    }
    break;
  case SELECT:
    end1 = fp1 + size;
    end2 = fp2 + size;
    for (i = 0; i < size-1; i++) {
      *tp = fselect(*fp1,*fp2,*end1,*end2);
      sum += *tp;
      tp++;
      fp1++;
      fp2++;
    }
    fp1++;
    fp2++;
    break;
  case ESMIN:
    for (i = 0; i < size; i++) {
      *tp = esminf(*fp1,*fp2,beta);
      sum += *tp;
      tp++;
      fp1++;
      fp2++;
    }
    break;
  case PSMIN:
    for (i = 0; i < size; i++) {
      *tp = psminf(*fp1,*fp2,beta);
      sum += *tp;
      tp++;
      fp1++;
      fp2++;
    }
    break;
  case PESMIN:
    for (i = 0; i < size; i++) {
      *tp = pesminf(*fp1,*fp2,beta);
      sum += *tp;
      tp++;
      fp1++;
      fp2++;
    }
    break;
  default:
   ;
  }
  tp = tmp;
  if (norm == 1) {
    for (i=0;i<size;i++) {
      *outp = (float)(*tp / sum);
      outp++;
      tp++;
    }
  }
  else {
    for (i=0;i<size;i++) {
      *outp = (float)*tp;
      tp++;
      outp++;
    }
  }
  return;
}

void 
pfile_combine(InFtrLabStream_PFile *data1,
	      InFtrLabStream_PFile *data2,
	      FILE *outfile,
	      int norm,
	      int format,
	      int rule,
	      float beta,
              int oswap)

{

  const int sep = -1;
  int num_phones = 0;
  const size_t n_sents1 = data1->num_segs();
  const size_t n_labs = data1->num_labs();
  const size_t n_ftrs1 =  data1->num_ftrs();
  const size_t n_frames1 = data1->num_frames();
 
  size_t n_sents2 = 0;
  size_t n_ftrs2 = 0;
  size_t n_frames2 = 0;
  int debug_level = 0;

  OutFtrLabStream_PFile* output=NULL;

  if (format == 0 && rule == SELECT) {
    output = new OutFtrLabStream_PFile(debug_level,"",outfile,
					  data1->num_ftrs()-1,
					  data1->num_labs(),1,oswap);
  }
  else if (format == 0) {
    output = new OutFtrLabStream_PFile(debug_level,"",outfile,
                                          data1->num_ftrs(),
                                          data1->num_labs(),1,oswap);
  }

  if (output == NULL) {
    error("Couldn't open output stream\n");
  }
  
  if (data2 != NULL) {
    n_sents2 = data2->num_segs();
    n_ftrs2 =  data2->num_ftrs();
    n_frames2 = data2->num_frames();

    if (n_sents1 != n_sents2)
      error("Different numbers of sentences in input files\n");
    if ( n_ftrs1 != n_ftrs2)
        error("Different numbers of features in input files\n");
    if (n_frames1 != n_frames2)
       error("Different numbers of frames in input files\n");
  }
  const size_t head = n_ftrs1;
  size_t buf_size = 2000;
  float *ftr_buf2;
  size_t n_read1;
  size_t n_read2;
  float *out_buf;

  if (rule == SELECT)
    out_buf = new float[(n_ftrs1-1) * buf_size];
  else
    out_buf = new float[n_ftrs1 * buf_size];

  if (format == 1)
    fwrite(&head,sizeof(head),1,outfile); // lna file header

  float *ftr_buf1 = new float[n_ftrs1 * buf_size];
  if (data2 != NULL)
    ftr_buf2 = new float[n_ftrs2 * buf_size];

  UInt32* lab_buf = new UInt32[n_labs * buf_size];

  //    cout << "Input file has " << n_ftrs << " features, " << n_sents <<      " sentences, " << n_frames << " frames\n";
  
  for (int r=0;r<data1->num_segs(); r++) {
    
    const size_t seg_frames = data1->num_frames(r);
    
    if ((r % 100) == 0)
      printf("Processing sentence %i\n",r);
    
    if (seg_frames == SIZET_BAD) 
      {
	error("Couldn't find number of frames "
	      "at sentence %lu in input pfile.\n",
	      (unsigned long) r);
      }
    const SegID seg_id = data1->set_pos(r, 0);
    if (data2 != NULL)
      data2->set_pos(r,0);
    
    if (seg_id == SEGID_BAD) {
      error("Couldn't seek to start of sentence %lu in input pfile.\n",r);
    }
    if (seg_frames > buf_size) {
      delete ftr_buf1;
      delete lab_buf;
      delete out_buf;
      buf_size = seg_frames * 2;
      ftr_buf1 = new float[buf_size * n_ftrs1];
      lab_buf = new UInt32[buf_size * n_labs];
      if (data2 != NULL) {
	delete ftr_buf2;
	ftr_buf2 = new float[buf_size * n_ftrs2];
      }
      if (rule == SELECT)
	out_buf = new float[buf_size * (n_ftrs1-1)];
      else
	out_buf = new float[buf_size * n_ftrs1];
    }
    n_read1 = data1->read_ftrslabs(seg_frames,ftr_buf1,lab_buf);
    if (data2 != NULL)
      n_read2 = data2->read_ftrs(seg_frames,ftr_buf2);
    
    if ( n_read1 != seg_frames) {
      error("At sentence %lu in input pfile, "
	      "only read %lu frames when should have read %lu.\n",
	      (unsigned long) r,
	    (unsigned long) n_read1, (unsigned long) n_frames1);
    }
    if (data2 != NULL && n_read1 != n_read2) 
      error("At sentence %lu in input pfile, "
	      "read %lu frames from file1 and  %lu from file2\n",
	      (unsigned long) r,
	      (unsigned long) n_read1, (unsigned long) n_read2);
    float *fp1  = ftr_buf1;
    float *fp2  = ftr_buf2;
    float *outp = out_buf;
    int s;

    // we do not combine, just change format
    if (data2 == NULL) {
      for ( s = 0; s < n_read1; s++) {
	if (format == 1)
	  fwrite(&s,sizeof(int),1,outfile);   // output framenumber for rapbin file
	if ((fwrite((float *)fp1,sizeof(float),n_ftrs1,outfile)) != n_ftrs1)
	  error("Fewer items written than requested\n");
	fp1 += n_ftrs1;
      }
    }
    // combination
    else if (data2 != NULL) {
      // for all sements
      for ( s = 0; s < n_read1; s++) {
	if (format == 1)
	  fwrite(&s,sizeof(int),1,outfile);
	combine(fp1,fp2,outp,n_ftrs1,norm,rule,beta);
	if (format == 1) {
	  if ((fwrite((float *)outp,sizeof(float),n_ftrs1,outfile))\
	      != n_ftrs1)
	    error("Fewer items written than requested\n");
	}
	else 
	  output->write_ftrslabs(1,outp,&lab_buf[s]);
	outp += n_ftrs1;	  
	fp1 += n_ftrs1;
	fp2 += n_ftrs1;
      }
    }
	if (format == 1)
	  fwrite((int *)&sep,sizeof(sep),1,outfile); 
	else
	  output->doneseg((SegID)r); 
      }
      delete ftr_buf1;
      if (data2 != NULL)
	delete ftr_buf2;
      delete lab_buf;
      if (output != NULL)
	delete output;
      
      return;
}


main(int argc, const char **argv) {

  const char *inname1 = NULL,*outname = NULL,*inname2 = NULL; 
  const char *priorsname1 = NULL, *priorsname2 = NULL;
  int debug_level = 0;
  FILE *infile1=NULL,*infile2=NULL,*priorsfile1 = NULL,*priorsfile2 = NULL;
  FILE *mapfile = NULL;
  int format = 0; // default: pfile
  int rule = NONE; // combination rule: default none
  const char *program_name = *argv++;
  const char *mapname = NULL;
  float beta=0.0;
  bool iswap1 = false,iswap2 = false, oswap = false;
  int norm = 0;
  int joint_prob = 0;
  const char *prob_table_fname = NULL;
  FILE *prob_table_fp = NULL;
  
  argc--;
  
  if (argc == 0)
    error(USAGESTRING);
  
  while (argc--)
    {
      char buf[BUFSIZ];
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
      else if (strcmp(argp,"-i2") == 0)
	{
	  if (argc>0)
	    {
	      inname2 = *argv++;
	      argc--;
	    }
	  else
	    error("Second input pfile name unspecified\n");
	}
      else if (strcmp(argp,"-o")==0)
	{
	  if (argc>0)
	    {
	      outname = *argv++;
	      argc--;
	    }
	  else
	    error("No output pfile name specified\n");
	}
      else if (strcmp(argp, "-debug")==0)
	{
	  // debug level
	  if (argc>0)
	    {
	      debug_level = atoi(*argv++);
	      argc--;
	    }
	  else
	    error("No debug level given\n");
	}
      else if (strcmp(argp,"-n") == 0)
	{
	  norm = 1;
	}
      else if (strcmp(argp,"-b") == 0)
	{
	  beta = atof(*argv++);
	  argc--;
	  printf("beta is %f\n",beta);
	}
      else if (strcmp(argp,"-p") == 0)
	{
	  rule = PROD;
	}
      else if (strcmp(argp,"-s") == 0)
	{
	  rule = SUM;
	}
      else if (strcmp(argp,"-max") == 0)
	{
	  rule = MAX;
	}
      else if (strcmp(argp,"-min") == 0)
	{
	  rule = MIN;
	}
      else if (strcmp(argp,"-smin") == 0)
	{
	  rule = SMIN;
	}
      else if (strcmp(argp,"-psmin") == 0)
	{
	  rule = PSMIN;
	}
      else if (strcmp(argp,"-esmin") == 0)
	{
	  rule = ESMIN;
	}
      else if (strcmp(argp,"-pesmin") == 0)
	{
	  rule = PESMIN;
	}
      else if (strcmp(argp,"-sel") == 0)
	{
	  rule = SELECT;
	  
	}	  
      else if (strcmp(argp,"-v") == 0)
	{
	  verbose = 1;
	}
      else if (strcmp(argp,"-f") == 0)
	{
	  if (argc>0) {
	    if (strcmp(*argv,"rapbin") == 0) 
	      format = 1;
	    else if (strcmp(*argv,"pfile") == 0)
	      format = 0;
	    else
	      error("Unrecognized output file format\n");
	    argc--;
	    argv++;
	  }
	  else
	    error("No argument for format specified\n");
	}
       else if (strcmp(argp,"-iswap1") == 0)
         {
           iswap1 = true;
          }
       else if (strcmp(argp,"-iswap2") == 0)
         {
           iswap2 = true;
          }

       else if (strcmp(argp,"-oswap") == 0)
         {
           oswap = true;
         }
       else {
	 error("Unrecognized argument (%s).",argp);
      }
    }
  
  if (rule == NONE)
    rule = SUM;
  if (inname1 == NULL) 
    error("No input pfile name supplied\n");

  if ((infile1 = fopen(inname1, "r")) == NULL)
    error("Couldn't open input pfile %s for reading.",inname1);

  if (inname2 != NULL) {
    if ((infile2 = fopen(inname2,"r")) == NULL)
      error("Couldn't open input pfile for reading: %s",inname2);
  }
  if (prob_table_fname != NULL)
    if ((prob_table_fp = fopen(prob_table_fname,"r")) == NULL)
      error("Couldn't open prob table for input\n");

  InFtrLabStream_PFile* input2;
  InFtrLabStream_PFile* input1
    = new InFtrLabStream_PFile(debug_level,"",infile1,1,iswap1);

  if (infile2 != NULL) {
   input2 = new InFtrLabStream_PFile(debug_level,"",infile2,1,iswap2);
  }
  else
    input2 = NULL;
  
  FILE *outfile = stdout;
  if (outname != NULL) {
    outfile = fopen(outname,"w");
    if (outfile == NULL)
      error("Couldn't open output pfile\n");
  }

  pfile_combine(input1,input2,outfile,norm,format,rule,beta,oswap);
  delete input1;

  if (input2 != NULL)
    delete input2;
  fclose(infile1);
  if (outname != NULL)
    fclose(outfile);
  if (infile2 != NULL)
    fclose(infile2);
  return 0;
}
    





