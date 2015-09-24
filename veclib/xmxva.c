/*****************************************************************************
 * xmxva() - Matrix - Vector Multiply w/skips.
 *
 * This following function computes the matrix-vector product C = A * B.
 *
 *      mxva(A,iac,iar,B,ib,C,ic,nra,nca)
 *
 *      A   ... double* ... matrix factor (input)
 *      iac ... int     ... increment in A between column elements
 *      iar ... int     ... increment in A between row elements
 *      B   ... double* ... vector factor (input)
 *      ib  ... int     ... increment in B between consecutive elements
 *      C   ... double* ... vector product (output)
 *      ic  ... int     ... increment in C between consecutive elements
 *      nra ... int     ... number of rows in A
 *      nca ... int     ... number of columns in A
 *
 * Consider BLAS2 xgemv as alternatives.
 *
 * $Id: xmxva.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <stdio.h>
#include <cfemdef.h>
#include <cveclib.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dmxva(double* A, int_t iac, int_t iar, double* B, int_t ib,
	   double* C, int_t ic,  int_t nra, int_t nca)
{
  register double  *a, *b,
                   *c = C;
  register double  sum;
  register int_t i, j;


  for (i = 0; i < nra; ++i) {
    sum = 0.0;
    a   =  A;
    A  += iac;
    b   =  B;
    for (j = 0; j < nca; ++j) {
      sum += (*a) * (*b); 
      a   += iar;
      b   += ib;
    }

    *c  = sum; 
     c += ic;
  }
}


void smxva (float* A, int_t iac, int_t iar, float* B, int_t ib,
	    float* C, int_t ic,  int_t nra, int_t nca)
{
  register float   *a, *b,
                   *c = C;
  register float   sum;
  register int_t i, j;

  for (i = 0; i < nra; ++i) {
    sum = 0.0F;
    a   = A;
    A  += iac;
    b   = B;
    for (j = 0; j < nca; ++j) {
      sum += (*a) * (*b); 
      a   += iar;
      b   += ib;
    }

    *c  = sum; 
     c += ic;
  }
}
