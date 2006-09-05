/**
 *: config.h
 *
 *  Gang Ji, http://welcome.to/rainier
 */

#ifndef CONFIG_H
#define CONFIG_H


/*
#ifndef DEGUGGING
#define DEBUGGING
#endif
*/

#define USE_CHOLESKY

#define DOUBLE_PROCESSION_DEFINED
//define SINGLE_PROCESSION_DEFINED


//#define USE_FULL_COVARIANCE_MATRICES
#define USE_DIAGONAL_COVARIANCE_MATRICES


#ifndef DATA_TYPE_DEFINED
#define DATA_TYPE_DEFINED
typedef float DataType;
#endif


#ifndef PROC_TYPE_DEFINED
#define PROC_TYPE_DEFINED

#ifdef DOUBLE_PROCESSION_DEFINED
//typedef double ProcType;
#else
//typedef float ProcType;
#endif		/* ifdef DOUBLE_PROCESSION_DEFINED */
#endif		/* ifndef PROC_TYPE_DEFINED */

typedef unsigned ProcType;

#endif		/* #ifndef CONFIG_H */



