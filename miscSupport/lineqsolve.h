//
// 
//  Copyright (C) 2001 Jeff Bilmes
//  Licensed under the Open Software License version 3.0
//
// $Header$
//
// Linear equation solver routines.
// 
//          Jeff Bilmes
//          bilmes@ee.washington.edu


extern void 
lubksb(float *a, 
       int n, 
       int *indx, 
       float *b);

extern void 
ludcmp(float *a,
       int n, 
       int *indx,
       float *d);

extern void
lineqsolve(const int n, const int nrhs,
	   float *a,
	   float *b);

extern void 
lubksb(double *a, 
       int n, 
       int *indx, 
       double *b);

extern void 
ludcmp(double *a,
       int n, 
       int *indx,
       double *d);

extern void
lineqsolve(const int n, const int nrhs,
	   double *a,
	   double *b);
