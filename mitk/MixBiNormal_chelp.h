/*
** return the value of P(X|l,G) without the normalizing
** constant, for each of the pairs 
** (x0[0],x1[0]),(x0[s],x1[s]), ... (x0[(n-1)*s],x1[(n-1)*s])
** where s is the stride.
** res{,a} and scratch must be length >= n
** This function corresponds to:
** void MixBiNormal::pg(double *x0,double *x1,int l,int n,double alpha,
**                      double *scratch,double *res,double *resa);
**
*/
void
pg_c1(double *x0,double *x1,
		int l,
		int n,
		double alpha,
		double *scratch,
		double *res,
		double *resa,
		const double u0,
		const double u1,
		const double s0,
		const double s2,
		const double ntw_s1,
		const double n_inv_det_2,
		const double sqrt_inv_det,
		const double norm
      );

/*
** This function corresponds to:
** void MixBiNormal::pg(double x0,double *x1,int l,int n,double alpha,
**                      double *scratch,double *res);
*/
void
pg_c2(double x0,double *x1,int l,int n,double alpha,
		   double *scratch,double *res,
		   const double u0,
		   const double u1,
		   const double s0,
		   const double s2,
		   const double ntw_s1,
		   const double d0,
		   const double n_inv_det_2,
		   const double sqrt_inv_det,
		   const double d0d0s2,
		   const double d0ntw_s1,
      const double norm);

/*
** This function corresponds to:
** void MixBiNormal::pg(double x0,double *x1,int l,int n,double alpha,
**    	                double *scratch1,double *scratch2,double *resa);
*/
void
pg_c3(double x0,
      double *x1,
      int l,
      int n,
      double alpha,
      double *scratch1,
      double *scratch2,
      double *resa,
      const double u0,
      const double u1,
      const double s0,
      const double s2,
      const double ntw_s1,
      const double d0,
      const double n_inv_det_2,
      const double sqrt_inv_det,
      const double d0d0s2,
      const double d0ntw_s1,
      const double norm);
