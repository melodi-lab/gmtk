//
// Log arithmetic .cc and test file.
//
//
// Written by: Jeff Bilmes
//             bilmes@icsi.berkeley.edu



#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "general.h"
VCID("$Header$");

#include "logp.h"

double logp_minLogExp = -log(-LZERO);

// These constants are used to determine when probabilities
// get to small to bother using them to multiply observations.
// We add an additional constant on to account for potentially
// small observations.
#define LOG_FLT_MIN_CNST 2.0
#define LOG_DBL_MIN_CNST 2.0
double log_FLT_MIN = log(FLT_MIN)+LOG_FLT_MIN_CNST;
double log_DBL_MIN = log(DBL_MIN)+LOG_DBL_MIN_CNST;

#ifdef _TABLE_
#include <iostream.h>
template <class FT, class iFT = double>
  FT logp<FT, iFT>::table[table_size];
template <class FT, class iFT = double>
  bool logp<FT, iFT>::initialized=false;
template <class FT, class iFT = double>
  double logp<FT, iFT>::inc;

template <class FT, class iFT>
void logp<FT, iFT>::table_init()
{
    inc = (table_size-1)/logp_minLogExp;
    for (int i=0; i<table_size; i++)
        table[i] = log(1.0 + exp(logp_minLogExp*i/table_size));
}
#endif

#ifdef MAIN

inline logpr foo(logpr p1, logpr p2)
{
  logpr tmp = p1*p2;
  return tmp;
}

int
main()
{
  {
  logp<float> a(0.5),b(0.25);
  printf("a=%f b=%f\n",a.unlog(),b.unlog());
  printf("a+b=%f\n",(a+b).unlog());
  printf("2*b=%f\n",(2.0*b).unlog());
  printf("b*2=%f\n",(b*2.0).unlog());
  printf("a+2*b=%f\n",(a+2.0*b).unlog());
  printf("a*b=%f\n",(a*b).unlog());
  printf("a-b=%f\n",(a-b).unlog());
  printf("a/b=%f\n",(a/b).unlog());
  printf("a==b=%d\n",(a==b));
  printf("a!=b=%d\n",(a!=b));
  printf("a<b=%d\n",(a<b));
  printf("a<=b=%d\n",(a<=b));
  printf("a>b=%d\n",(a>b));
  printf("a>=b=%d\n",(a>=b));
  }
  
  {
  logpr a(0.5),b(0.25);
  printf("a=%f b=%f\n",a.unlog(),b.unlog());
  printf("a+b=%f\n",(a+b).unlog());
  printf("2*b=%f\n",(2.0*b).unlog());
  printf("b*2=%f\n",(b*2.0).unlog());
  printf("a+2*b=%f\n",(a+2.0*b).unlog());
  printf("a*b=%f\n",(a*b).unlog());
  printf("a-b=%f\n",(a-b).unlog());
  printf("a/b=%f\n",(a/b).unlog());
  printf("a==b=%d\n",(a==b));
  printf("a!=b=%d\n",(a!=b));
  printf("a<b=%d\n",(a<b));
  printf("a<=b=%d\n",(a<=b));
  printf("a>b=%d\n",(a>b));
  printf("a>=b=%d\n",(a>=b));
  }
  
  logpr tmp(0,log(0.3));
  logp<float,float> tmp2(0,log(0.3));
  logp<double,float> tmp3(0,log(0.3));
  logp<float,double> tmp4(0,log(0.3));
  
  printf("tmp's = %f,%f,%f,%f\n",
	 tmp.unlog(),
	 tmp2.unlog(),
	 tmp3.unlog(),
	 tmp4.unlog());

  printf("sizeof(logp<float,float>) = %d\n",sizeof(logp<float,float>));
  printf("sizeof(logp<double,float>) = %d\n",sizeof(logp<double,float>));
  printf("sizeof(logp<float,double>) = %d\n",sizeof(logp<float,double>));
  printf("sizeof(logp<double,double>) = %d\n",sizeof(logp<double,double>));


  foo(tmp,tmp);

}


#endif
