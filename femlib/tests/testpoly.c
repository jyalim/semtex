/*****************************************************************************
 * TESTPOLY.C: test suite for polyops.c.
 *****************************************************************************/

static char
RCSid[] = "$Id: testpoly.c,v 8.1 2015/04/20 11:14:14 hmb Exp $";

#include <stdio.h>
#include <math.h>

#include <cveclib>
#include <cfemlib>


void polcoe(int n, double *x, double *y, double *c);
void polder(int n, double *c, double  x, double *poly, double *pder);

#define  VEC_MAX 64





int main()
/* ========================================================================= *
 * Driver routine.
 * ========================================================================= */
{
  double *c, *z, *y, *w, **D, **DT;
  double  x, xcheb, PI_n, poly, pder;
  int     i, j, n, np;


  c  = dvector(0, VEC_MAX-1);
  z  = dvector(0, VEC_MAX-1);
  y  = dvector(0, VEC_MAX-1);
  w  = dvector(0, VEC_MAX-1);
  D  = dmatrix(0, VEC_MAX-1, 0, VEC_MAX-1);
  DT = dmatrix(0, VEC_MAX-1, 0, VEC_MAX-1);

  /* ----------------------------------------------------------------------- *
   * Test routine that finds Gauss-Lobatto-Jacobi points (and, indirectly,   *
   * the one that evaluates Jacobi polynomials) by finding the Gauss-        *
   * Lobatto-Chebyshev points, which are available in closed form.           *
   * Routine is from Canuto, App. C.                                         *
   * ----------------------------------------------------------------------- */

  printf("\n---test value--- TEST JACGL ---true value---\n");

  for (i=2; i<=6; i++) {
    n = i << 1;
    np = n + 1;
    PI_n = M_PI / (double) n;

    jacgl(n, -0.5, -0.5, z);
    
    printf("\nChebyshev Points, N=%1d\n", n);
    for (j=0; j<np; j++) {
      xcheb = -cos(j*PI_n);
      printf("%20.17g\t%20.17g\n", z[j], xcheb);
    }
  }

  /* ----------------------------------------------------------------------- * 
   * Test the Gauss-Legendre points and weights: we have some closed-form    *
   * results available, e.g. Hughes \S 3.8. We will just try a 3-point rule. *
   * ----------------------------------------------------------------------- */

  printf("\n---test value--- TEST  ZWGL ---true value---\n");

  zwgl(z, w, 3);

  printf("\nPoints for NP=3\n");

  printf("%20.17g\t%20.17g\n", z[0], -sqrt(3.0/5.0));
  printf("%20.17g\t%20.17g\n", z[1], 0.0);
  printf("%20.17g\t%20.17g\n", z[2],  sqrt(3.0/5.0));

  printf("\nWeights for NP=3\n");

  printf("%20.17g\t%20.17g\n", w[0], 5.0/9.0);
  printf("%20.17g\t%20.17g\n", w[1], 8.0/9.0);
  printf("%20.17g\t%20.17g\n", w[2], 5.0/9.0);

  /* ----------------------------------------------------------------------- *
   * Test the production of derivative-operator matrices by comparing dgll,  *
   * dermat_k and dermat_g (which should return the same results, barring    *
   * rounding errors) to "analytic" derivative obtained by extracting the    *
   * coefficients of the interpolating polynomials & differentiating by      *
   * synthetic division.                                                     *
   * ----------------------------------------------------------------------- */

  printf("\nTESTING DERIVATIVE OPERATORS\n");
  jacgl(4, 0.0, 0.0, z);

  printf("*** analytic (5x5):\n");

  for (i=0; i<5; i++) {
    dzero(5, w, 1);
    w[i] = 1.0;
    polcoe(4, z, w, c); 
    for (j=0; j<5; j++)
      polder(4, c, z[j], &poly, DT[j] + i);
  }
  for (i=0; i<5; i++) {
    for (j=0; j<5; j++)
      printf("%23.17g", DT[i][j]);
    printf("\n");
  }

  
  printf("*** dgll (5x5):\n");

  dgll(5, z, D, DT);

  for (i=0; i<5; i++) {
    for (j=0; j<5; j++)
      printf("%23.17g", D[i][j]);
    printf("\n");
  }

  printf("*** dermat_k (5x5):\n");

  dermat_k(5, z, D, DT);

  for (i=0; i<5; i++) {
    for (j=0; j<5; j++)
      printf("%23.17g", D[i][j]);
    printf("\n");
  }

  printf("*** dermat_g (5x5):\n");

  dermat_g(5, z, 5, z, D, DT);

  for (i=0; i<5; i++) {
    for (j=0; j<5; j++)
      printf("%23.17g", D[i][j]);
    printf("\n");
  }

  /* ----------------------------------------------------------------------- *
   * Test intmat_g by comparing results to "analytic" evaluations.           *
   * Start with 5 points evenly-spaced on [-1, 1], interpolate to 4 internal *
   * Gauss-Legendre points.                                                  *
   * ----------------------------------------------------------------------- */

  printf("\nTESTING INTERPOLATION OPERATOR\n");

  z[0] = -(z[4] = 1.0);
  z[1] = -(z[3] = 0.5);
  z[2] = 0.0;

  jacg(3, 0.0, 0.0, w);


  printf("*** analytic: (4x5)\n");

  for (i=0; i<5; i++) {		/* 5 Lagrange interpolants. */
    dzero(5, y, 1);
    y[i] = 1.0;
    polcoe(4, z, y, c); 
    for (j=0; j<4; j++)		/* 4 points to interpolate at. */
      polder(4, c, w[j], DT[j]+i, &pder);
  }
  for (i=0; i<4; i++) {
    for (j=0; j<5; j++)
      printf("%23.17g", DT[i][j]);
    printf("\n");
  }

  
  printf("*** intmat_g (4x5):\n");

  intmat_g(5, z, 4, w, D, DT);

  for (i=0; i<4; i++) {
    for (j=0; j<5; j++)
      printf("%23.17g", D[i][j]);
    printf("\n");
  }

}





void polcoe(int n, double *x, double *y, double *c)
/* ========================================================================= *
 * polcoe(): calculate polynomial coefficients given knot points.            *
 *                                                                           *
 * Given arrays x[0..n] and y[0..n] containing tabulated function y_i=f(x_i) *
 * this routine returns an array of coefficients c[0..n], such that          *
 * y_i = c_j x_i^j.                                                          *
 *                                                                           *
 * Numerical Recipes /S 3.5.  See also the caveats there.                    *
 * ========================================================================= */
{
  int     i, j, k;
  double  phi, ff, b, *s;


  s = dvector(0, n);
  dzero(n+1, s, 1);
  dzero(n+1, c, 1);

  s[n] -= x[0];
  for (i=1; i<=n; i++) {
    for (j=n-i; j<=n-1; j++)
      s[j] -= x[i]*s[j+1];
    s[n] -= x[i];
  }

  for (j=0; j<=n; j++) {

    phi = n + 1;
    for (k=n; k>=1; k--)
      phi = k*s[k] + x[j]*phi;
    ff = y[j] / phi;

    b = 1.0;
    for (k=n; k>=0; k--) {
      c[k] += b * ff;
      b     = s[k] + x[j]*b;
    }

  }

  freeDvector(s, 0);
}





void polder(int n, double *c, double x, double *poly, double *pder)
/* ========================================================================= *
 * Given the coefficients of a polynomial c[0..n] and an abcissa x, return   *
 * the value of the polynomial as *poly and the value of the first derivat-  *
 * ive as *pder at location x.  Method generalizes to arbitrary derivatives: *
 * see Numerical Recipes, /S 5.3.   c[0] = constant term.                    *
 * ========================================================================= */
{
  int     j;
  double  p, dp;

  p  = c[j=n];
  dp = 0.0;

  while (j--) {
    dp = dp*x + p;
    p  =  p*x + c[j];
  }

  *poly = p;
  *pder = dp;
}

    
