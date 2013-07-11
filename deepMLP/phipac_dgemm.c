/*
** phipac_dgemm.c:
**      BLAS compatible DGEMM interface for PHiPAC matrix-matrix
**      multiply code.
** 
** Written by:
**      CheeWhye Chin <cheewhye@icsi.berkeley.edu>
**
**
** "Copyright (c) 1997 The Regents of the University of California.  All
** rights reserved."  Permission to use, copy, modify, and distribute
** this software and its documentation for any purpose, without fee, and
** without written agreement is hereby granted, provided that the above
** copyright notice and the following two paragraphs appear in all copies
** of this software.
** 
** IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
** FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
** ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
** THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
** SUCH DAMAGE.
**  
** THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
** INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
** MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
** PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
** CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
** ENHANCEMENTS, OR MODIFICATIONS.
**
**/

#include <stdio.h>

#define max(a,b)     (((a) > (b)) ? (a) : (b))


/* function prototype for PHiPAC generated
 * matrix multiplication code
 */
extern void
mm_double_NN_1(int, int, int,
	       double*, double*, double*,
	       int, int, int,
	       double);

extern void
mm_double_NN_c(int, int, int, 
	       double*, double*, double*,
	       int, int, int,
	       double, double);

extern void
mm_double_NT_1(int, int, int,
	       double*, double*, double*,
	       int, int, int,
	       double);

extern void
mm_double_NT_c(int, int, int, 
	       double*, double*, double*,
	       int, int, int,
	       double, double);

extern void
mm_double_TN_1(int, int, int,
	       double*, double*, double*,
	       int, int, int,
	       double);

extern void
mm_double_TN_c(int, int, int, 
	       double*, double*, double*,
	       int, int, int,
	       double, double);

extern void
mm_double_TT_1(int, int, int,
	       double*, double*, double*,
	       int, int, int,
	       double);

extern void
mm_double_TT_c(int, int, int, 
	       double*, double*, double*,
	       int, int, int,
	       double, double);


/* error handling routine */
void
xerbla(char* srname, int* info)
{
  fprintf(stderr,
	  "** On entry to %6s parameter number %2u had an illegal value\n",
	  srname, *info);
}


void
dgemm(char* transA, char* transB,
      int* M, int* N, int* K,
      double* alpha,
      double* A, int* Astride,
      double* B, int* Bstride,
      double* beta,
      double* C, int* Cstride)
{
  int info = 0;

  /* error checking */

  if (*M < 0)
    {
      info = 3;
      xerbla ("DGEMM ", &info);
      return;
    }
  else if (*N < 0)
    {
      info = 4;
      xerbla ("DGEMM ", &info);
      return;
    }
  else if (*K < 0)
    {
      info = 5;
      xerbla ("DGEMM ", &info);
      return;
    }

  if (*transA == 'n' || *transA == 'N')
    {
      if (*transB == 'n' || *transB == 'N')
	{
	  /* error checking */
	  if (*Astride < max(1,*M))
	    {
	      info = 8;
	      xerbla ("DGEMM ", &info);
	    } 
	  else if (*Bstride < max(1,*K)) 
	    {
	      info = 10;
	      xerbla ("DGEMM ", &info);
	    } 
	  else if (*Cstride < max(1,*M))
	    {
	      info = 13;
	      xerbla ("DGEMM ", &info);
	    }
	  else if (*alpha == 1)
	    {
	      /* C = AB + beta*C */
	      mm_double_NN_1(*N, *K, *M,
			     B, A, C,
			     *Bstride, *Astride, *Cstride,
			     *beta);
	    }
	  else
	    {
	      /* C = alpha*AB + beta*C */
	      mm_double_NN_c(*N, *K, *M,
			     B, A, C, 
			     *Bstride, *Astride, *Cstride,
			     *alpha, *beta);
	    }
	} 
      else
	{
	  /* error checking */
	  if (*transB != 'c' && *transB != 'C' &&
	      *transB != 't' && *transB != 'T')
	    {
	      info = 2;
	      xerbla ("DGEMM ", &info);
	    } 
	  else if (*Astride < max(1,*M))
	    {
	      info = 8;
	      xerbla ("DGEMM ", &info);
	    }
	  else if (*Bstride < max(1,*N))
	    {
	      info = 10;
	      xerbla ("DGEMM ", &info);
	    }
	  else if (*Cstride < max(1,*M))
	    {
	      info = 13;
	      xerbla ("DGEMM ", &info);
	    }
	  else if (*alpha == 1)
	    {
	      /* C = AB' + beta*C */
	      mm_double_TN_1(*N, *K, *M,
			     B, A, C, 
			     *Bstride, *Astride, *Cstride,
			     *beta);
	    }
	  else
	    {
	      /* C = alpha*AB' + beta*C */
	      mm_double_TN_c(*N, *K, *M,
			     B, A, C,
			     *Bstride, *Astride, *Cstride,
			     *alpha, *beta);
	    }
	}
    }
  else
    {
      if (*transB == 'n' || *transB == 'N')
	{
	  /* error checking */
	  if (*transA != 'c' && *transA != 'C' &&
	      *transA != 't' && *transA != 'T')
	    {
	      info = 1;
	      xerbla ("DGEMM ", &info);
	    }
	  else if (*Astride < max(1,*K))
	    {
	      info = 8;
	      xerbla ("DGEMM ", &info);
	    }
	  else if (*Bstride < max(1,*K))
	    {
	      info = 10;
	      xerbla ("DGEMM ", &info);
	    } 
	  else if (*Cstride < max(1,*M))
	    {
	      info = 13;
	      xerbla ("DGEMM ", &info);
	    }
	  else if (*alpha == 1)
	    {
	      /* C = A'B + beta*C */
	      mm_double_NT_1(*N, *K, *M,
			     B, A, C,
			     *Bstride, *Astride, *Cstride, 
			     *beta);
	    }
	  else
	    {
	      /* C = alpha*A'B + beta*C */
	      mm_double_NT_c(*N, *K, *M, 
			     B, A, C,
			     *Bstride, *Astride, *Cstride,
			     *alpha, *beta);
	    }
	}
      else
	{
	  /* error checking */
	  if (*transA != 'c' && *transA != 'C' &&
	      *transA != 't' && *transA != 'T') 
	    {
	      info = 1;
	      xerbla ("DGEMM ", &info);
	    } 
	  else if (*transB != 'c' && *transB != 'C' &&
		   *transB != 't' && *transB != 'T')
	    {
	      info = 2;
	      xerbla ("DGEMM ", &info);
	    } 
	  else if (*Astride < max(1,*K)) 
	    {
	      info = 8;
	      xerbla ("DGEMM ", &info);
	    } 
	  else if (*Bstride < max(1,*N))
	    {
	      info = 10;
	      xerbla ("DGEMM ", &info);
	    } 
	  else if (*Cstride < max(1,*M))
	    {
	      info = 13;
	      xerbla ("DGEMM ", &info);
	    }
	  else if (*alpha == 1)
	    {
	      /* C = A'B' * beta*C */
	      mm_double_TT_1(*N, *K, *M,
			     B, A, C,
			     *Bstride, *Astride, *Cstride,
			     *beta);
	    }
	  else
	    {
	      /* C = alpha*A'B' * beta*C */
	      mm_double_TT_c(*N, *K, *M,
			     B, A, C,
			     *Bstride, *Astride, *Cstride,
			     *alpha, *beta);
	    }
	}
    }
}
