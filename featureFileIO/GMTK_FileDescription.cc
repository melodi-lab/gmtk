/*
 *
 * Written by Katrin Kirchhoff <katrin@ee.washington.edu>
 *
 * Copyright (c) 2001
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any non-commercial purpose
 * and without fee is hereby granted, provided that the above copyright
 * notice appears in all copies.  The University of Washington,
 * Seattle make no representations about the suitability of this soft
 * for any purpose. It is provided "as is" without express or implied warranty.
 *
 */

#include <stdio.h>
#include <string.h>
#include "GMTK_FileDescription.h"

FileDescription::FileDescription(const char *name, const char *crng_str,
                                 const char *drng_str,
    	                         unsigned *nfloats, unsigned *nints, 
		                 unsigned *format, bool *swap, unsigned num) {


   if (name == NULL) 	
       error("FileDescription: File name is NULL for file %i\n",num);	

   fofName = new char[strlen(name)+1];
   strcpy(fofName,name);
   
   if (format == NULL)
     error("FileDescription:Data format unspecified for file %i\n",num);
   
   dataFormat = *format;
   
   if (nfloats == NULL)
     error("FileDescription: Number of floats unspecified for file %i\n",num);

   nFloats = *nfloats;

   if (nints == NULL) 
      error("FileDescription: Number of ints unspecified for file %i\n",num);

   nInts = *nints;

   cont_rng = new BP_Range(crng_str,0,nFloats);
   disc_rng = new BP_Range(drng_str,0,nInts);

   nFloatsUsed = cont_rng->length();
   nIntsUsed = disc_rng->length();

   if (dataFormat == PFILE) {

   if ((curDataFile = fopen(fofName,"rb")) == NULL)
	error("FileDescription: Can't open '%s' for input\n",
		   fofName);

   pfile_istr = new InFtrLabStream_PFile(0,fofName,curDataFile,1,*swap);

   fofSize = pfile_istr->num_segs();

   if (pfile_istr->num_ftrs() != nFloats) 
	error("FileDescription: File %s has %i floats, expected %i\n",
		   fofName,
		   pfile_istr->num_ftrs(),
		   nFloats);

    // in a pfile, only the labs can be ints

   if (pfile_istr->num_labs() != nInts)
	error("FileDescription: File %s has %i floats, expected %i\n",
		   fofName,
		   pfile_istr->num_labs(),
		   nInts);

    }

   // if not a pfile we are dealing with a file of file names

    else {
      if ((fofFile = fopen(fofName,"r")) == NULL)
	error("FileDescription: Can't open '%s' for input\n",fofName);
      
     size_t fsize = fileSize(fofFile);

     fofBuf = new char[fsize];
      
     size_t n_read = fread((char *)fofBuf,1,fsize,fofFile);

     if (n_read < fsize) 
	error("FileDescription: Only read %li bytes from file %s, 
                    expected %li for file %s\n", n_read,fofName,fsize);

     dataNames = new char*[MAXFILES];
     fofSize = readFof();

     pfile_istr = NULL;

     fclose(fofFile);
  }
  if (swap == NULL)
    bswap = false; //default
  else
    bswap = *swap; 


}


FileDescription::~FileDescription() {

  printf("deleting file description\n");

  delete [] fofName;
  delete [] fofBuf;
  delete [] dataNames;

  if (cont_rng != NULL)
    delete cont_rng;	
  if (disc_rng != NULL)
    delete disc_rng;

  
}
	
size_t
FileDescription::readFof() {

  size_t n_lines = 0;
  char *tmp2, *endp;
  char prev;
  int i;

  tmp2 = fofBuf;
  dataNames[n_lines] = tmp2;
  endp = fofBuf;

  while (*endp != '\0')
    endp++;


  prev = *tmp2;

  while (++tmp2 != endp) {
    
    if (isspace(*tmp2) && (!(isspace(prev)))) {
      prev = *tmp2;
      *tmp2 = '\0';
      printf("%s\n",dataNames[n_lines]);
    }
    else if (isspace(prev) && (!(isspace(*tmp2)))) {
      dataNames[++n_lines] = tmp2;
      prev = *tmp2;
    }
    else 
      prev = *tmp2;
  }
  
  // add 1 count if final newline if missing 
  
  printf("prev: %c last: '%c'\n",prev,*tmp2);
  
  
  return n_lines;
}




	








