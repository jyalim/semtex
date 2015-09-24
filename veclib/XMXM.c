/*****************************************************************************
 * xmxm() - Matrix multiply.
 *
 * The following function computes the matrix product C = A * B.
 * The matrices are assumed to be in row-major order and must occupy
 * consecutive memory locations.  The input quantities are as follows:
 *
 *     mxm(A,nra,B,nca,C,ncb)
 *
 *     A    ... double* ... source vector one
 *     nra  ... int     ... number of rows in A
 *     B    ... double* ... source operand two
 *     nca  ... int     ... number of columns in A
 *     C    ... double* ... result vector
 *     ncb  ... int     ... number of columns in B
 *
 * $Id: XMXM.c,v 8.1 2015/04/20 11:14:19 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>


void dmxm (double* A, int_t nra,
	   double* B, int_t nca,
	   double* C, int_t ncb)
{
  register double  *a = A,
                   *b = B,
                   *c = C;
  register double  sum;
  register int_t i, j, k;

  for (i = 0; i < nra; i++) {
    for (j = 0; j < ncb; j++) {

      b   = B + j;                  /* Next column of B    */
      sum = 0.0;                    /* Clear sum           */

      for (k = 0; k < nca; k++)
	sum += a[k] * b[k*ncb];     /* Inner product loop  */

      *c++   = sum;                 /* Store and increment */
    }
    a += nca;                       /* Next row of A       */
  }
}


void smxm (float* A, int_t nra,
	   float* B, int_t nca,
	   float* C, int_t ncb)
{
  register float   *a = A,
                   *b = B,
                   *c = C;
  register float   sum;
  register int_t i, j, k;

  for (i = 0; i < nra; i++) {
    for (j = 0; j < ncb; j++) {

      b   = B + j;
      sum = 0.0F;

      for(k = 0; k < nca; k++)
	sum += a[k] * b[k*ncb];

      *c++   = sum;
    }
    a += nca;
  }
}
