/*
 * $Header$
 * Eigenanalysis
 *
 * Converted to "C" by
 *  Jeff Bilmes <bilmes@icsi.berkeley.edu>
 *
 */


#ifndef EIG_H
#define EIG_H

extern void
eigenanalyasis(int n,
	       double *cor, /* nxn real symmetric covariance matrix */
	       double *vals,   /* n-vector, space for eigenvalues */
	       double *vecs); 

#endif
