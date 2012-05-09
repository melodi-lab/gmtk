#include <stdio.h>
#include <stdlib.h>
#include "GMTK_ObservationMatrix.h"




#define USAGESTRING "  Usage: iotest   options       \n   options                    defaults \n   -n    number of files        0\n   -ifile list of input files  none\n   -iformat list of input formats none\n   -nf      list of number of floats none\n   -ni      list of number of ints none\n   -ext     list of output extensions  none\n   -icrng   list of input cont ranges  none\n   -idrng   list of input disc. ranges none\n   -ocrng   list of output cont ranges none\n   -odrng   list of output disc ranges none\n   -iswap   list of input bytes swap flags none\n   -oswap  list of output byte swap flags  none \n   -o num output files 0 \n   -oformat list of output formats none\n \n"



void
string2strings(char *s1, char **sarray, int maxn) {

  int nf = 0;
  char *cp = s1;
  sarray[nf] = s1;
  while (*cp != '\0') {
    if (*cp == '-') {
      *cp = '\0';
      if (++nf >= maxn) {
	fprintf(stderr,"More items than expected in string %s\n",s1);
	exit(-1);
      }
      sarray[nf] = ++cp;
    }
    else 
      cp++;
  }
}

void
string2ints(char *s1, unsigned *sarray, int maxn) {

  int nf = 0;
  char *cp = s1;
  char *tmp = s1;
  
  while (*cp != '\0') {
    if (*cp == '-') {
      *cp = '\0';
      sarray[nf] = atoi(tmp);
      if (++nf >= maxn) {
        fprintf(stderr,"More items than expected in string %s\n",s1);
        exit(-1);
      }
      tmp = ++cp;
    }
    else 
      cp++;
  }
  sarray[nf] = atoi(tmp);
}

int
main(int argc, const char **argv) {

  int num_files = 0;
  int num_outfiles = 0;

  ObservationMatrix buf;

  //  DataOutStream *out_str = NULL;

  const char **ap;

  const char *infile_str = NULL;
  const char *iformat_str = NULL;
  const char *iswap_str = NULL;
  const char *infloat_str = NULL;
  const char *inint_str = NULL;
  const char *icrange_str = NULL;
  const char *idrange_str = NULL;

  const char *ofile_str = NULL;
  const char *oformat_str = NULL;
  const char *oswap_str = NULL;
  const char *outext_str = NULL;
  const char *ocrange_str = NULL;
  const char *odrange_str = NULL;

  
  int n_args = argc;

  if (n_args == 1) {
    fprintf(stderr,USAGESTRING);
    exit(-1);
  }

  ap = argv;
    ap++;
    n_args--;

  while (n_args > 0) {
    if (strcmp(*ap,"-n") == 0) {
      num_files = atoi(*(++ap));
      n_args--;
    }
    else if (strcmp(*ap,"-ifile") == 0) {
      infile_str = *(++ap);
      n_args--;
    }
    else if (strcmp(*ap,"-ofile") == 0) {
      ofile_str = *(++ap);
      n_args--;
    }
    else if (strcmp(*ap,"-iformat") == 0) {
      iformat_str = *(++ap);
      n_args--;
    }
    else if (strcmp(*ap,"-nf") == 0) {
      infloat_str = *(++ap);
      n_args--;
    }
    else if (strcmp(*ap,"-ni") == 0) {
      inint_str = *(++ap);
      n_args--;
    }
    else if (strcmp(*ap,"-ext") == 0) {
      outext_str = *(++ap);
      n_args--;
    }
    else if (strcmp(*ap,"-ocrng") == 0) {
      ocrange_str = *(++ap);
      n_args--;
    }
    else if (strcmp(*ap,"-odrng") == 0) {
      odrange_str = *(++ap);
      n_args--;
    }
    else if (strcmp(*ap,"-icrng") == 0) {
      icrange_str = *(++ap);
      n_args--;
    }
    else if (strcmp(*ap,"-idrng") == 0) {
      idrange_str = *(++ap);
      n_args--;
    }
    else if (strcmp(*ap,"-iswap") == 0) {
      iswap_str = *(++ap);
      n_args--;
    }
    else if (strcmp(*ap,"-oswap") == 0) {
      oswap_str = *(++ap);
      n_args--;
    }
    else if (strcmp(*ap,"-o") == 0) {
      num_outfiles = atoi(*(++ap));
      n_args--;
    }
    else if (strcmp(*ap,"-oformat") == 0) {
      oformat_str = *(++ap);
      printf("oformat str: %s\n",oformat_str);
      n_args--;
    }
    else {
      fprintf(stderr,"Unsupported argument: %s\n",*ap);
      exit(-1);
    }
    ap++;
    n_args--;

  }

  if (infile_str == NULL) {
    fprintf(stderr,"File names must be specified\n");
    exit(-1);
  }
  if (infloat_str == NULL) {
    fprintf(stderr,"Number of floats must be specified\n");
    exit(-1);
  }
  if (inint_str == NULL) {
    fprintf(stderr,"Number of ints must be specified\n");
    exit(-1);
  }
  if (iformat_str == NULL) {
    fprintf(stderr,"File formats must be specified\n");
    exit(-1);
  }

  
  char **in_names = new char*[num_files];
  char **format_str = new char*[num_files];
  unsigned *iformats = new unsigned[num_files];
  unsigned *in_nfloats = new unsigned[num_files];
  unsigned *in_nints = new unsigned[num_files];
  char **icrange = new char*[num_files];
  char **idrange = new char*[num_files];
  bool *iswap = new bool[num_files];

  char **out_names = new char*[num_outfiles];
  char **out_formats = new char*[num_outfiles];
  unsigned *oformats = new unsigned[num_outfiles];
  unsigned *out_ints = new unsigned[num_outfiles];
  unsigned *out_floats = new unsigned[num_outfiles];
  char **ocrange = new char*[num_outfiles];
  char **odrange = new char*[num_outfiles];
  BP_Range **ocrng = new BP_Range*[num_outfiles];
  BP_Range **odrng = new BP_Range*[num_outfiles];
  bool *oswap = new bool[num_outfiles];
  char **out_ext = new char*[num_outfiles];


  if (num_outfiles > 0 && (ocrange_str == NULL || odrange_str == NULL)) {
    fprintf(stderr,"Output ranges must be specified\n");
    exit(-1);
  }
  
  /* parse infile names */

  string2strings((char *)infile_str,in_names,num_files);

  /* parse file format string */

  string2strings((char *)iformat_str,format_str,num_files);

  for (int i= 0; i < num_files; i++) {
    char *tmp = format_str[i];
    if (strcmp(tmp,"htk") == 0)
      iformats[i] = HTK;
   else if (strcmp(tmp,"binary") == 0)
      iformats[i] = RAWBIN;
    else if (strcmp(tmp,"ascii") == 0)
      iformats[i] = RAWASC;
    else if (strcmp(tmp,"pfile") == 0)
      iformats[i] = PFILE;
    else {
      printf("Unknown file format: '%s'\n",tmp);
      exit(-1); 
    }
  }
  
  
  /* parse float string */
  
  string2ints((char *)infloat_str,in_nfloats,num_files);
  
  /* parse int string */
  
  string2ints((char *)inint_str,in_nints,num_files);
  
  /* parse range strings */
  
  if (icrange_str != NULL) 
    string2strings((char *)icrange_str,icrange,num_files);
  else
    icrange = NULL;
  
  if (idrange_str != NULL)
    string2strings((char *)idrange_str,idrange,num_files);
  else
    idrange = NULL;
  
  if (ocrange_str != NULL) {
    printf("parsing ocrange string\n");
    string2strings((char *)ocrange_str,ocrange,num_outfiles);
  }
  else
    ocrange = NULL;
  
  if (odrange_str != NULL)
    string2strings((char *)odrange_str,odrange,num_outfiles);
  else
    odrange = NULL;
  
  /* parse byteswap strings */
  
  if (iswap_str != NULL)
    string2ints((char *)iswap_str,(unsigned *)iswap,num_files);
  else
    iswap = NULL;
  
  if (oswap_str != NULL)
    string2ints((char *)oswap_str,(unsigned *)oswap,num_outfiles);
  else
    oswap = NULL;

  if (outext_str != NULL)
    string2strings((char *)outext_str,out_ext,num_outfiles);
  else
    out_ext = NULL;

  if (oformat_str != NULL)
    string2strings((char *)oformat_str,out_formats,num_outfiles);
  else
    out_formats = NULL;

  
  printf("Command line summary for %i input files: \n",num_files);
  for (int f = 0; f < num_files; f++) {
    printf("file %s: nfloats: %i nints: %i format: %i ",
	   in_names[f],in_nfloats[f],in_nints[f],iformats[f]);
    if (icrange_str != NULL)
      printf("icrange: %s\n",icrange[f]);
    if (idrange_str != NULL)
      printf("idrange: %s\n",idrange[f]);
    if (iswap_str != NULL)
      printf("byteswap: %i",iswap[f]);
    printf("\n");
  }

  if (buf.active() == 0)
    printf("Observation matrix not active\n");
  
  buf.openFiles(num_files,
		(const char **)in_names,
		(const char **)icrange,
		(const char **)idrange,
		in_nfloats,
		in_nints,
		iformats,
		iswap);

  if (buf.active())
    printf("Observation matrix is active now\n");
  
  for (size_t f = 0; f < buf.numSegments() ; f++) {

    printf("Utterance %i\n",f);
    
    buf.loadSegment(f);

    buf.printSegmentInfo();
    
    // print first 3 frames 

    for (size_t a = 0; a < buf.numFrames; a++) 
      buf.printFrame(stdout,a);
    
  }


  delete [] in_names;
  delete [] format_str;
  delete [] iformats;
  delete [] in_nfloats;
  delete [] in_nints;
  delete [] icrange;
  delete [] idrange;
  delete [] iswap;

  delete [] out_names;
  delete [] out_formats;
  delete [] oformats;
  delete [] out_ints;
  delete [] out_floats;
  delete [] ocrange;
  delete [] odrange;
  delete [] ocrng;
  delete [] odrng;
  
  delete [] oswap;
  delete [] out_ext;
    
  return 0;
}



