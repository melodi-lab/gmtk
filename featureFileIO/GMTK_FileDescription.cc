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
#include "GMTK_io.h"
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

  if (swap == NULL)
    bswap = false; //default
  else
    bswap = *swap; 

   if (dataFormat == PFILE) {

     if ((curDataFile = fopen(fofName,"rb")) == NULL)
       error("FileDescription: Can't open '%s' for input\n",
	     fofName);

     pfile_istr = new InFtrLabStream_PFile(0,fofName,curDataFile,1,bswap);
     
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

   else {
     pfile_istr = NULL;

     if ((fofFile = fopen(fofName,"r")) == NULL)
       error("FileDescription: Can't open '%s' for input\n",fofName);

     fofSize = readFof(fofFile);

     fclose(fofFile);
   }
}


FileDescription::~FileDescription() {

  if (fofName != NULL)
    delete [] fofName;

  if (dataFormat != PFILE) {
    for (unsigned i = 0; i < fofSize; i++)
      delete [] dataNames[i];
    delete [] dataNames;
  }
  else 
    delete pfile_istr;

  if (cont_rng != NULL)
    delete cont_rng;   
  if (disc_rng != NULL)
    delete disc_rng;
}
	
size_t
FileDescription::readFof(FILE *f) {

  size_t n_lines = 0,maxlines = 0;
  char line[MAXSTRLEN];

  while (fgets(line,sizeof(line),f) != NULL)
    maxlines++;
  rewind(f);

  dataNames = new char*[maxlines];

  n_lines = 0;
  
  while (fgets(line,sizeof(line),f) != NULL) {
    int l = strlen(line);
    if (line[l-1] != '\n') {
      if (n_lines < maxlines-1) 
	error("GMTK_FileDescription::readFof: line %i too long in file '%s' - increase MAXSTRLEN\n",
	      n_lines+1,fofName);
    }
    else
      line[l-1] = '\0';
    
    dataNames[n_lines] = new char[l+1];

    strcpy(dataNames[n_lines],line);

    n_lines++;
  }
  return n_lines;
}




	








