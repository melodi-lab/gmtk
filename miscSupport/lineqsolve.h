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
