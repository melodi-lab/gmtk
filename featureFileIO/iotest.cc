#include <stdio.h>
#include <stdlib.h>
#include "GMTK_io.h"

#define MAXINPUTFILES 10


main(int argc, const char **argv) {

  DataInStream *in_str;
  const char **ap;
  const char **names = new const char*[MAXINPUTFILES];
  const char **format_str = new const char*[MAXINPUTFILES];
  unsigned *formats = new unsigned[MAXINPUTFILES];
  unsigned *n_floats = new unsigned[MAXINPUTFILES];
  unsigned *n_ints = new unsigned[MAXINPUTFILES];
  char **crange_str = new char*[MAXINPUTFILES];
  char **drange_str = new char*[MAXINPUTFILES];
  bool *swaps = new bool[MAXINPUTFILES];
  char **out_names = new char*[MAXINPUTFILES];
  char **out_ext = new char*[MAXINPUTFILES];
  unsigned *out_floats  = new unsigned[MAXINPUTFILES];
  unsigned *out_ints = new unsigned[MAXINPUTFILES];
  unsigned *out_formats = new unsigned[MAXINPUTFILES];
  char **out_crange = new char*[MAXINPUTFILES];
  char **out_drange = new char*[MAXINPUTFILES];
  char *out_frange;
  DataOutStream *out_str;

  for (int m = 0; m < MAXINPUTFILES; m++) {
    out_names[m] = new char[MAXSTRLEN];
    out_ext[m] = new char[MAXSTRLEN];
    crange_str[m] = new char[MAXSTRLEN];
    drange_str[m] = new char[MAXSTRLEN];
  }

  ap = argv;
  ap++;
  int nf = 0;

   for (int a = 1; a < argc-1; a++,ap++) {
    printf("%s\n",*ap);
    names[nf] = *ap;

    n_floats[nf] = atoi(*(++ap));
    a++;

    if (n_floats[nf] > 0) {
      strcpy(crange_str[nf],*(++ap));
      a++;
   }
    else
      crange_str[nf] = "all";

    n_ints[nf] = atoi(*(++ap));
    a++;

    if (n_ints[nf] > 0) {
      strcpy(drange_str[nf],*(++ap));
      a++;
    }
    else
      drange_str[nf] = "all";

    printf("ranges: %s %s\n",crange_str[nf],crange_str[nf]);

   format_str[nf] = *(++ap);
   printf("%s\n",format_str[nf]);
   a++;
   if (strcmp(format_str[nf],"htk") == 0)
      formats[nf] = HTK;
    else if (strcmp(format_str[nf],"binary") == 0)
      formats[nf] = RAWBIN;
    else if (strcmp(format_str[nf],"ascii") == 0)
      formats[nf] = RAWASC;
    else if (strcmp(format_str[nf],"pfile") == 0)
      formats[nf] = PFILE;
    else {
      printf("Unknown file format: '%s'\n",format_str[nf]);
      exit(-1); 
    }
    swaps[nf] = atoi(*(++ap)); 
    a++;
	
    nf++;
  }
  printf("done\n");


  printf("found %i files\n",nf);
  

  out_formats[0] = RAWBIN;
  out_formats[1] = RAWASC;
  out_ext[0] = "outbin";
  out_ext[1] = "outasc";

  out_crange[0] = "0:12";
  //  out_drange[1] = "13:25";

  out_frange = "0:199";

   for (int j = 0; j < nf; j++)
     printf("%i %i %i %i %s\n",
	    n_floats[j],n_ints[j],formats[j],swaps[j],out_ext[j]);


      

  in_str = new DataInStream(nf,names,(const char **)crange_str,(const char **)drange_str,
                    	    n_floats,n_ints,formats,swaps);

  in_str->getNumFloats(out_floats);
  in_str->getNumInts(out_ints);

  out_str = new DataOutStream(nf,
			      out_floats,
			      out_ints,
			      out_formats,
			      NULL,
			      NULL,
			      NULL,
			      NULL);

  for (size_t f = 0; f < in_str->numSegs() ; f++) {

    printf("Utterance %li\n",f);

    size_t num_frames = in_str->readFile(f);

    printf("read  %li frames ok\n",num_frames);

    in_str->getFileNames(f,out_names);
    out_str->initFiles(out_names,out_ext,num_frames);
    out_str->writeData(f,out_crange,NULL,out_frange,in_str->getObs());
    out_str->closeFiles();
  }

  delete out_str;
  delete in_str;

  for (int i = 0; i < nf; i++) {
    delete [] out_names[i];
    delete [] out_ext[i];
    delete [] format_str[i];
    delete [] names[i];
    delete [] crange_str[i];
    delete [] drange_str[i];
  }

  delete [] ap;
  delete [] names;
  delete [] format_str;
  delete [] crange_str;
  delete [] drange_str;
  delete [] out_names;
  delete [] out_ext;
  delete [] out_crange;
  delete [] out_drange;
  delete [] out_formats;
  delete [] out_ints;
  delete [] out_floats;
  delete [] swaps;
  delete [] n_ints;
  delete [] n_floats;
  delete [] formats;
  
}



