/*	Copyright (c) 1988 AT&T */
/*	  All Rights Reserved */

/*	THIS IS UNPUBLISHED PROPRIETARY SOURCE CODE OF AT&T */
/*	The copyright notice above does not evidence any */
/*	actual or intended publication of such source code. */

/*	Modifications	Copyright ICL 1988,1989 */
/*	ieeefp.h	ICL version  7.1	(89/03/29) */

#ifndef _IEEEFP_H
#define	_IEEEFP_H

#pragma ident	"@(#)ieeefp.h	1.14	94/07/29 SMI"	/* SVr4.0 1.18 */

#ifdef	__cplusplus
extern "C" {
#endif

/*
 * Floating point enviornment for machines that support
 * the IEEE 754 floating-point standard.  This file currently
 * supports the 80*87, and SPARC families.
 *
 * This header defines the following interfaces:
 *	1) Classes of floating point numbers
 *	2) Rounding Control
 *	3) Exception Control
 *	4) Exception Handling
 *	5) Utility Macros
 *	6) Full Exception Environment Control
 */

/*
 * Note that the following construct, "#machine()", is a non-standard
 * extension to ANSI-C.  It is maintained here to provide compatibility
 * for existing compilations systems, but should be viewed as transitional
 * and may be removed in a future release.  If it is required that this
 * file not contain this extension, edit this file to remove the offending
 * condition.
 */
#if defined(__STDC__)
#define	__sparc_conditional #machine(sparc) || defined(__sparc)
#else
#define	__sparc_conditional defined(sparc) || defined(__sparc)
#endif

#if defined(__STDC__)
#define	__i386_conditional #machine(i386) || defined(__i386)
#else
#define	__i386_conditional defined(i386) || defined(__i386)
#endif

#if defined(__STDC__)
#define	__ppc_conditional #machine(ppc) || defined(__ppc)
#else
#define	__ppc_conditional defined(__ppc)
#endif

/*
 * CLASSES of floating point numbers *************************
 * IEEE floating point values fall into 1 of the following 10
 * classes
 */
typedef	enum	fpclass_t {
	FP_SNAN = 0,	/* signaling NaN */
	FP_QNAN = 1,	/* quiet NaN */
	FP_NINF = 2,	/* negative infinity */
	FP_PINF = 3,	/* positive infinity */
	FP_NDENORM = 4, /* negative denormalized non-zero */
	FP_PDENORM = 5, /* positive denormalized non-zero */
	FP_NZERO = 6,	/* -0.0 */
	FP_PZERO = 7,   /* +0.0 */
	FP_NNORM = 8,	/* negative normalized non-zero */
	FP_PNORM = 9	/* positive normalized non-zero */
} fpclass_t;

#if defined(__STDC__)
extern fpclass_t fpclass(double);	/* get class of double value */
extern int	finite(double);
extern int	unordered(double, double);
#else
extern fpclass_t fpclass();	/* get class of double value */
#endif

/*
 * ROUNDING CONTROL ******************************************
 *
 * At all times, floating-point math is done using one of four
 * mutually-exclusive rounding modes.
 */

#if __i386_conditional

/*
 * NOTE: the values given are chosen to match those used by the
 * 80*87 rounding mode field in the control word.
 */
typedef	enum	fp_rnd {
	FP_RN = 0,	/* round to nearest representable number, tie -> even */
	FP_RM = 1,	/* round toward minus infinity */
	FP_RP = 2,	/* round toward plus infinity */
	FP_RZ = 3	/* round toward zero (truncate) */
} fp_rnd;

#endif

#if __sparc_conditional || __ppc_conditional

/*
 * NOTE: the values given are chosen to match those used by the
 * RD (Round Direction) field of the FSR (Floating Point State Register).
 */
typedef	enum	fp_rnd {
	FP_RN = 0,	/* round to nearest representable number, tie -> even */
	FP_RZ = 1,	/* round toward zero (truncate) */
	FP_RP = 2,	/* round toward plus infinity */
	FP_RM = 3	/* round toward minus infinity */
} fp_rnd;

#endif

#if defined(__STDC__)
extern fp_rnd	fpsetround(fp_rnd);	/* set rounding mode, return previous */
extern fp_rnd	fpgetround(void);	/* return current rounding mode */

#else
extern fp_rnd	fpsetround();		/* set rounding mode, return previous */
extern fp_rnd	fpgetround();		/* return current rounding mode */

#endif

/*
 * EXCEPTION CONTROL *****************************************
 *
 */

#define	fp_except	int

#define	FP_DISABLE	0	/* exception will be ignored */
#define	FP_ENABLE	1	/* exception will cause SIGFPE */
#define	FP_CLEAR	0	/* exception has not occurred */
#define	FP_SET		1	/* exception has occurred */

#if __i386_conditional

/*
 * There are six floating point exceptions, which can be individually
 * ENABLED (== 1) or DISABLED (== 0).  When an exception occurs
 * (ENABLED or not), the fact is noted by changing an associated
 * "sticky bit" from CLEAR (==0) to SET (==1).
 *
 * NOTE: the bit positions in fp_except are chosen to match those of
 * the 80*87 control word mask bits.  Although the 87 chips actually
 * ENABLE exceptions with a mask value of 0 (not 1, as on the 3b), it
 * is felt that switching these values may create more problems than
 * it solves.
 */

/* an fp_except can have the following (not exclusive) values: */
#define	FP_X_INV	0x01	/* invalid operation exception */
#define	FP_X_DNML	0x02	/* denormalization exception */
#define	FP_X_DZ		0x04	/* divide-by-zero exception */
#define	FP_X_OFL	0x08	/* overflow exception */
#define	FP_X_UFL	0x10	/* underflow exception */
#define	FP_X_IMP	0x20	/* imprecise (loss of precision) */

#endif

#if __sparc_conditional || __ppc_conditional

/*
 * There are five floating-point exceptions, which can be individually
 * ENABLED (== 1) or DISABLED (== 0).  When an exception occurs
 * (ENABLED or not), the fact is noted by changing an associated
 * "sticky bit" from CLEAR (==0) to SET (==1).
 *
 * NOTE: the bit positions in an fp_except are chosen to match that in
 * the Trap Enable Mask of the FSR (Floating Point State Register).
 */

/* an fp_except can have the following (not exclusive) values: */
#define	FP_X_INV	0x10	/* invalid operation exception */
#define	FP_X_OFL	0x08	/* overflow exception */
#define	FP_X_UFL	0x04	/* underflow exception */
#define	FP_X_DZ		0x02	/* divide-by-zero exception */
#define	FP_X_IMP	0x01	/* imprecise (loss of precision) */

#endif

#if defined(__STDC__)
extern fp_except fpgetmask(void);		/* current exception mask */
extern fp_except fpsetmask(fp_except);		/* set mask, return previous */
extern fp_except fpgetsticky(void);		/* return logged exceptions */
extern fp_except fpsetsticky(fp_except);	/* change logged exceptions */

#else
extern fp_except fpgetmask();	/* current exception mask */
extern fp_except fpsetmask();	/* set mask, return previous */
extern fp_except fpgetsticky();	/* return logged exceptions */
extern fp_except fpsetsticky();	/* change logged exceptions */

#endif

/*
 * UTILITY MACROS ********************************************
 */

#if defined(__STDC__)
extern int isnanf(float);
extern int isnand(double);

#else
extern int isnand();
#define	isnanf(x)	(((*(long *)&(x) & 0x7f800000L) == 0x7f800000L) && \
			    ((*(long *)&(x) & 0x007fffffL) != 0x00000000L))
#endif

#if __i386_conditional

#ifndef SS
#include <sys/reg.h>
#endif
/*
 * EXCEPTION HANDLING ****************************************
 *
 * When a signal handler catches an FPE, it will have a freshly initialized
 * coprocessor.  This allows signal handling routines to make use of
 * floating point arithmetic, if need be.  The previous state of the 87
 * chip is available, however.  There are two ways to get at this information,
 * depending on how the signal handler was set up.
 *
 * If the handler was set via signal() or sigset(), the old, SVR3, method
 * should be used: the signal handler assumes that it has a single parameter,
 * which is of type struct _fpstackframe, defined below.  By investigating
 * this parameter, the cause of the FPE may be determined.  By modifying it,
 * the state of the coprocessor can be changed upon return to the main task.
 * THIS METHOD IS OBSOLETE, AND MAY NOT BE SUPPORTED IN FUTURE RELEASES.
 *
 * If the handler was set via sigaction(), the new, SVR4, method should be
 * used: the third argument to the handler will be a pointer to a ucontext
 * structure (see sys/ucontext.h).  The uc_mcontext.fpregs member of the
 * ucontext structure holds the saved floating-point registers.  This can be
 * examined and/or modified.  By modifying it, the state of the coprocessor
 * can be changed upon return to the main task.
 */

struct _fpstackframe {		/* signal handler's argument */
	long signo;		/* signal number arg */
	long regs[SS+1];	/* all registers */
	struct _fpstate *fpsp;	/* address of saved 387 state */
	char *wsp;		/* address of saved Weitek state */
};

struct _fpreg {	/* structure of a temp real fp register */
	unsigned short significand[4];	/* 64 bit mantissa value */
	unsigned short exponent;	/* 15 bit exponent and sign bit */
};

struct _fpstate {		/* saved state info from an exception */
	unsigned long   cw,	/* cotrol word */
			sw,	/* status word after fnclex-not useful */
			tag,	/* tag word */
			ipoff,	/* %eip register */
			cssel,	/* code segment selector */
			dataoff, /* data operand address */
			datasel; /* data operand selector */
	struct _fpreg _st[8];	/* saved register stack */
	unsigned long status;	/* status word saved at exception */
};

/*
 * The structure of the 80*87 status and control words are given by the
 * following structs.
 */
struct _cw87 {
	unsigned    mask:   6,  /* exception masks */
		    res1:   2,  /* not used */
		    prec:   2,  /* precision control field */
		    rnd:    2,  /* rounding control field */
		    inf:    1,  /* infinity control (not on 387) */
		    res2:   3;  /* not used */
};

struct _sw87 {
	unsigned    excp:	6,	/* exception sticky bits */
		    res1:	1,	/* not used */
		    errs:	1,	/* error summary-set if unmasked excp */
		    c012:	3,	/* condition code bits 0..2 */
		    stkt:	3,	/* stack top pointer */
		    c3:		1,	/* condition code bit 3 */
		    busy:	1;	/* coprocessor busy */
};

#endif

#ifdef	__cplusplus
}
#endif

#endif	/* _IEEEFP_H */
