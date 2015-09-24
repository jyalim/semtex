/*****************************************************************************
 * polyops.c:  Routines for manipulating polynomials.
 *
 * Copyright (c) 1994 <--> $Date: 2015/04/20 11:14:14 $, Hugh Blackburn
  *
 * This file is part of Semtex.
 * 
 * Semtex is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 * 
 * Semtex is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Semtex (see the file COPYING); if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 *
 *
 * Summary of routines:
 * --------------------
 * dermat_g: Derivative operator for Lagrange interpolant, arbitrary points.
 * dermat_k: Derivative operator for Lagrange interpolant, at nodes (knots).
 * intmat_g: Interpolation operator for Lagrange interpolant, arbitrary.
 * jacobf  : Jacobi polynomial and its derivative.
 * jacg    : Points for Gauss-Jacobi quadrature.
 * jacgr   : Points for Gauss-Radau-Jacobi quadrature.
 * jacgl   : Points for Gauss-Lobatto-Jacobi quadrature.
 * zwgl    : Points and weights for Gauss-Legendre quadrature.
 * zwgrl   : Points and weights for Gauss-Radau-Legendre quadrature.
 * zwgll   : Points and weights for Gauss-Lobatto-Legendre quadrature.
 * pnleg   : Evaluate Legendre polynomial.
 * pndleg  : Evaluate derivative of Legendre polynomial.
 * pnd2leg : Evaluate second derivative of Legendre polynomial.
 * dgll    : Derivative operator for Gauss-Lobatto-Legendre interpolant.
 * uniknot : Points uniformly distributed on [-1, 1].
 *
 * Routines that deal specifically with orthogonal polynomials come from
 * a library of spectral routines written in FORTRAN by Einar Ronquist, MIT.
 * Many of the formulae used may be found in Canuto, Hussaini, Quarteroni &
 * Zang, "Spectral Methods in Fluid Dynamics", Springer, 1988.
 * The jacobf routine comes from Funaro, as that has been verified to work
 * also for alpha, beta != 0.0, 0.5.
 *
 * Everything here is real_t (double) precision.
 *
 * $Id: polyops.c,v 8.1 2015/04/20 11:14:14 hmb Exp $
 *****************************************************************************/

#include <math.h>
#include <stdio.h>

#include <cfemdef.h>
#include <cfemlib.h>
#include <cveclib.h>

#define STOP 16
static real_t gammaF(real_t x);

/* Row-major access in an array with rows of length n. */

#define RMA(i,j,n) (j+i*n) 

void dermat_g (const int_t   K   ,
	       const real_t* zero, 
	       const int_t   I   ,
	       const real_t* x   ,
	       real_t*       D   ,
	       real_t*       DT  )
/* ------------------------------------------------------------------------- *
 * Given a set of zeros for a Lagrange interpolant order K-1 and a set of I
 * points x, return the derivative operator matrix D and its transpose DT.
 *
 * The matrix D is IxK and DT is KxI, both are zero-offset.  D is defined
 * so that
 *                           [D]{y} = {y'},
 * that is, when it premultiplies a vector of points (ordinates) {y},
 * located at (abscissae) {x}, it returns (ordinates) {y'}, the derivative
 * of the Lagrange interpolant with knots {zero}.
 *
 * The set of points {x} may be identical to the knot points {zero} in which
 * case the matrix [D] maps the values {y} onto the derivative of their
 * interpolating polynomial evaluated at the same locations.   In this case,
 * however, the routine dermat_k(), which should have lower rounding errors,
 * may be used instead.
 *
 * Reference: Abramowitz & Stegun 25.3.2.
 * ------------------------------------------------------------------------- */
{
  register int_t i, j, k, l;
  register real_t* a;
  register real_t  sum, prod;

  a = dvector (0, K-1);
  dfill (K, 1.0, a, 1);

  for (k=0; k<K; k++)
    for (l=0; l<K; l++)
      if (l != k) a[k] *= zero[k] - zero[l];

   for (j=0; j<I; j++)
     for (k=0; k<K; k++) {
       sum = 0.0;
       for (l=0; l<K; l++) {
	 if (l != k) {
	   prod = 1.0;
	   for (i=0; i<K; i++)
	     if (i != k && i != l) prod *= x[j] - zero[i];
	   sum += prod;
	 }
       }
       D[RMA(j,k,K)]  = DT[RMA(k,j,I)] = sum / a[k];
       /* D[j][k]  = DT[k][j] = sum / a[k]; */
     }

  freeDvector (a, 0);
}


void dermat_k (const int_t   K   ,
	       const real_t* zero,
	       real_t*       D   ,
	       real_t*       DT  )
/* ------------------------------------------------------------------------- *
 * Given a set of zeros for a Lagrange interpolant order K-1, return the
 * (collocation) derivative operator matrix D and its transpose DT.
 *
 * The matrix D is KxK and is zero-offset.  D is defined so that
 *                           [D]{y} = {y'},
 * that is, when it premultiplies a vector of points (ordinates) {y},
 * located at (abscissae) {zero}, it returns (ordinates) {y'}, the deriv-
 * ative of the Lagrange interpolant with those knot points (roots).
 *
 * The fact that the evaluation points lie at the zeros of the Lagrange
 * functions simplifies the operations that must be performed when compared
 * to the more general routine dermat_g().
 *
 * It may be more accurate to use closed-form results for any given class of
 * e.g. orthogonal Lagrange interpolants if available.
 *
 * Reference: Solomonoff, A. & Turkel, E.  1989.  "Global properties of
 *   pseudospectral methods".  JCP 81, 239--276
 * ------------------------------------------------------------------------- */
{
  register int_t   j, k, l;
  register real_t* a;
  register real_t  sum, prod;

  a = dvector (0, K-1);
  dfill (K, 1.0, a, 1);

  for (k = 0; k < K; k++)
    for (l = 0; l < K; l++)
      if (l != k) a[k] *= zero[k] - zero[l];

   for (j = 0; j < K; j++)
     for (k = 0; k < K; k++)
       if (j == k) {		    /* -- Diagonal term:     use eq (2.5). */
	 sum = 0.0;
	 for (l = 0; l < K; l++)
	   if (l != k) sum += 1.0 / (zero[k] - zero[l]);
	 D[RMA(k,k,K)] = DT[RMA(k,k,K)] = sum;
       } else {			    /* -- Off-diagonal term: use eq (2.7). */
	 prod = 1.0;
	 for (l = 0; l < K; l++) 
	   if (l != j) prod *= zero[j] - zero[l];
	 D[RMA(j,k,K)] = DT[RMA(k,j,K)] = prod / ((zero[j] - zero[k]) * a[k]);
       }

  freeDvector (a, 0);
}


void intmat_g (const int_t   K   ,
	       const real_t* zero,
	       const int_t   I   ,
	       const real_t* x   ,
	       real_t*       IN  ,
	       real_t*       IT  )
/* ------------------------------------------------------------------------- *
 * Given a set of zeros for a Lagrange interpolant order K-1 and a set of I
 * points x, return the interpolation matrix IN and its transpose IT.
 *
 * The matrix IN is IxK and IT is KxI, both are zero-offset.  IN is defined
 * so that
 *                           [IN]{y} = {z},
 * that is, it maps a set of K (ordinates) {y}, given at the Lagrange knot
 * points (abscissae, roots) onto a set of I ordinates {z}, given at arbitr-
 * ary locations {x}
 *
 * Reference: Abramowitz & Stegun 25.2.2.
 * ------------------------------------------------------------------------- */
{
  register int_t   i, j, k, l;
  register real_t* a;
  register real_t  prod;

  a = dvector (0, K-1);
  dfill (K, 1.0, a, 1);

  for (k = 0; k < K; k++)
    for (l = 0; l < K; l++)
      if (l != k) a[k] *= zero[k] - zero[l];

  for (j = 0; j < I; j++)
    for (k = 0; k < K; k++) {
      prod = 1.0;
      for (i = 0; i < K; i++)
	if (i != k) prod *= x[j] - zero[i];
      IN[RMA(j,k,K)] = IT[RMA(k,j,I)] = prod / a[k];
      /* IN[j][k] = IT[k][j] = prod / a[k]; */
    }

  freeDvector (a, 0);
}  


static void JACOBF (const int_t  n     ,
		    const real_t x     ,
		    const real_t alpha ,
		    const real_t beta  ,
		    real_t*      poly  ,
		    real_t*      pder  ,
		    real_t*      polym1,
		    real_t*      pderm1,
		    real_t*      polym2,
		    real_t*      pderm2)
/* ------------------------------------------------------------------------- *
 * Computes the Jacobi polynomial (poly) of degree n, and its derivative
 * (pder), at location x.  Values for lower degree are also returned.
 *
 * n          :  degree of approximation,
 * alpha, beta:  parameters in Jacobi weight function
 *                                   alpha           beta
 *                    w(x) = (1 - x)^      * (1 + x)^
 *               where special cases alpha = beta =  0.0 <==> G-L-Legendre
 *                                   alpha = beta = -0.5 <==> G-L-Chebyshev.
 *
 * References:
 *   "FORTRAN Routines for Spectral Methods", D. Funaro, 1993
 * ------------------------------------------------------------------------- */
{
  register int_t  i;
  register real_t y, ym, yp, dy, dym, dyp, ys, dys, apb;
  register real_t di, c0, c1, c2, c3, c4;

  *poly = y  = 1.0;
  *pder = dy = 0.0;
  if (n == 0) return;

  apb  = alpha + beta;
  *polym1 = y;
  *pderm1 = dy;
  *poly = y  = 0.5 * (apb + 2.0) * x + 0.5 * (alpha - beta);
  *pder = dy = 0.5 * (apb + 2.0);
  if (n == 1) return;

  yp   = 1.0;
  dyp  = 0.0;
  for (i=2; i <= n; i++) {
    di  = (real_t)(i);

    c0  = 2.0 * di + apb;
    c1  = 2.0 * di * (di + apb) * (c0 - 2.0);
    c2  = (c0 - 1.0) * (c0 - 2.0) * c0;
    c3  = (c0 - 1.0) * (alpha - beta) * apb;
    c4  = 2.0 * (di + alpha - 1.0) * c0 * (di + beta - 1.0);
    
    ys  = yp;
    dys = dyp;
    ym  = y;
    y   = ( (c2 * x + c3) * y - c4 * yp) / c1;
    yp  = ym;
    dym = dy;
    dy  = ( (c2 * x + c3) * dy - c4 * dyp + c2 * yp) / c1;
    dyp = dym;
  }

  *poly   = y;
  *pder   = dy;
  *polym1 = yp;
  *pderm1 = dyp;
  *polym2 = ys;
  *pderm2 = dys;
}


void JACG (const int_t  n    ,
	   const real_t alpha,
	   const real_t beta ,
	   real_t*      xjac )
/* ------------------------------------------------------------------------- *
 * Compute Gauss points xjac for a polynomial of degree n on range [-1, 1].
 * Points are returned in ascending order.  N + 1 points are found.
 *
 * Alpha, beta = 0.0 => Gauss-Legendre; = -0.5 => Chebyshev.
 * ------------------------------------------------------------------------- */
{
  register int_t  i, j, k, np, nh;
  register real_t dth, recsum, x, delx;
  real_t          poly, pder, polym1, pderm1, polym2, pderm2;

  np  = n + 1;
  nh  = np >> 1;
  dth = M_PI / (np << 1);

  for (j = 0; j < nh; j++) {
    x = cos (((j<<1) + 1) * dth);

    k = 0;
    do {
      JACOBF (np, x, alpha, beta, &poly,&pder,&polym1,&pderm1,&polym2,&pderm2);

      recsum = 0.0;
      for (i = 0; i < j-1; i++) recsum += 1.0 / (x - xjac[n - i]);

      delx = -poly / (pder - recsum*poly);
      x   += delx;
    } while (fabs (delx) > EPSDP && k++ < STOP);
    xjac [n - j] = x;
  }

  for (i = 0; i < nh; i++) xjac[i] = -xjac[n - i];

  if (np & 1) xjac[nh] = 0.0;
}


void JACGR (const int_t  n    ,
	    const real_t alpha,
	    const real_t beta ,
	    real_t*      xjac )
/* ------------------------------------------------------------------------- *
 * Computes the Gauss-Radau points xjac for polynomial of degree n on
 * [-1, 1].  Points are returned in ascending order.
 *
 * The GRL rule is asymmetric: this version assumes the desired
 * end-point is at z=-1. If you'd wanted the endpoint at +1, instead,
 * the ordering and sign of the points should be reversed.
 *
 * Alpha, beta = 0.0 => Gauss-Legendre; = -0.5 => Chebyshev.
 * ------------------------------------------------------------------------- */
{
  register int_t  np, i, j, k;
  register real_t x, delx, con, recsum;
  real_t          pn, pdn, pnp1, pdnp1, pnm1, pdnm1, func, funcd;

  np  = n + 1;
  con = 2.0 * M_PI / (n<<1 + 1);

  for (j = 0; j < np; j++) {
    x = -cos (con * j);

    k = 0;
    do {
      JACOBF (np, x, alpha, beta, &pn, &pdn, &pnp1, &pdnp1, &pnm1, &pdnm1);
      func  = pn  + pnp1;
      funcd = pdn + pdnp1;

      recsum = 0.0;
      for (i = 0; i < j-1; i++) recsum += 1.0 / (x - xjac[i]);

      delx  = -func  / (funcd - recsum*func);
      x    += delx;
    } while (fabs (delx) > EPSDP && k++ < STOP);

    xjac[j] = x;
  }
}


void JACGL (const int_t  n    ,
	    const real_t alpha,
	    const real_t beta ,
	    real_t*      xjac )
/* ------------------------------------------------------------------------- *
 * Computes the Gauss-Lobatto points xjac for polynomial of degree n on
 * [-1, 1].  Points are returned in ascending order.
 *
 * Alpha, beta = 0.0 => Gauss-Legendre, = -0.5 => Chebyshev.
 * ------------------------------------------------------------------------- */
{
  register int_t  i, j, k, np, jm, nh;
  register real_t a, b, det, con, rp, rm, x, delx, recsum;
  real_t          poly, pder, pnp1p, pdnp1p, pnp, pdnp, pnm1p, pdnm1;
  real_t          pnp1m, pdnp1m, pnm, pdnm, pnm1m;

  np      = n + 1;
  nh      = np >> 1;
  xjac[0] = -1.0;
  xjac[n] =  1.0;
  con     = M_PI / (real_t) n;

  JACOBF (np,  1.0, alpha, beta, &pnp1p, &pdnp1p, &pnp, &pdnp, &pnm1p, &pdnm1);
  JACOBF (np, -1.0, alpha, beta, &pnp1m, &pdnp1m, &pnm, &pdnm, &pnm1m, &pdnm1);
  det = pnp*pnm1m - pnm*pnm1p;
  rp  = -pnp1p;
  rm  = -pnp1m;
  a   = (rp*pnm1m - rm*pnm1p) / det;
  b   = (rm*pnp   -   rp*pnm) / det;

  for (j = 1; j < n; j++) {
    jm = j - 1;
    x  = cos (con * j);

    k  = 0;

    do {
      JACOBF (np, x, alpha,beta, &pnp1p, &pdnp1p, &pnp, &pdnp, &pnm1p, &pdnm1);
      poly = pnp1p  + a*pnp  + b*pnm1p;
      pder = pdnp1p + a*pdnp + b*pdnm1;
      
      recsum = 0.0;
      for (i = 0; i < jm; i++) recsum += 1.0 / (x - xjac[n - i]);

      delx = -poly / (pder - recsum*poly);
      x   += delx;
    } while (fabs (delx) > EPSDP && k++ < STOP);

    xjac[n - j] = x;
  }
}


void ZWGL (real_t*     z ,
	   real_t*     w ,
	   const int_t np)
/* ------------------------------------------------------------------------- *
 * Gauss-Legendre points and weights.
 *
 * Generate np G-L points (z) and weights (w) for integration over the
 * range [-1, 1] for a polynomial of degree n = np - 1.
 *
 * Reference: Canuto et al., eq (2.3.10).
 * ------------------------------------------------------------------------- */
{
  register int_t i, n;
  real_t         poly, pder, polym1, pderm1, polym2, pderm2;

  if (np < 2) { z[0] = 0.0;  w[0] = 2.0; return; }
  
  n  = np - 1;
  JACG (n, 0.0, 0.0, z);

  for (i = 0; i < np; i++) {
    JACOBF (np, z[i], 0.0, 0.0, &poly,&pder,&polym1,&pderm1,&polym2,&pderm2);
    w[i] = 2.0 / ( (1.0 - SQR(z[i])) * SQR(pder) );
  }
}
  

void ZWGRL (real_t*     z ,
	    real_t*     w ,
	    const int_t np)
/* ------------------------------------------------------------------------- *
 * Gauss-Radau-Legendre points and weights.
 *
 * Generate np G-R-L points (z) and weights (w) for integration over the
 * range [-1, 1] for a polynomial of degree n = np - 1.
 *
 * Reference: Canuto et al., eq (2.3.11).
 * ------------------------------------------------------------------------- */
{
  register int_t i, n;
  real_t         poly, pder, polym1, pderm1, polym2, pderm2, con;

  if (np < 2) { z[0] = 0.0; w[0] = 2.0; return; }

  n   = np - 1;
  con = 1.0 / SQR(np);

  JACGR (n, 0.0, 0.0, z);
  
  w[0] = 2.0 * con;

  for (i = 1; i < np; i++) {
    JACOBF (n, z[i], 0.0, 0.0, &poly,&pder,&polym1,&pderm1,&polym2,&pderm2);
    w[i] = con * (1.0 * z[i]) / SQR(poly);
  }
}
   

void ZWGLL (real_t*     z ,
	    real_t*     w ,
	    const int_t np)
/* ------------------------------------------------------------------------- *
 * Gauss-Lobatto-Legendre points and weights.
 *
 * Generate np G-L-L points (z) and weights (w) for integration over the
 * range [-1, 1] for a polynomial of degree n = np - 1.
 *
 * Reference: Canuto et al., eq (2.3.12).
 * ------------------------------------------------------------------------- */
{
  register int_t i, n;
  real_t         poly, pder, polym1, pderm1, polym2, pderm2, con;

  if (np  < 2) { z[0] = w[0] = 0.0; return; }

  if (np == 2) { z[0] = -(z[1] = w[0] = w[1] = 1.0); return; }

  n   = np - 1;
  con = 2.0 / (real_t) (n * np);
  
  JACGL (n, 0.0, 0.0, z);

  for (i = 0; i < np; i++) {
    JACOBF (n, z[i], 0.0, 0.0, &poly,&pder,&polym1,&pderm1,&polym2,&pderm2);
    w[i] = con / (SQR(poly));
  }
}
   

void ZWGLJ (real_t*      z ,
	    real_t*      w ,
	    const real_t alpha,
	    const real_t beta ,
	    const int_t  np)
/* ------------------------------------------------------------------------- *
 * Gauss-Lobatto-Jacobi points and weights, for Jacobi constants alpha & beta.
 *
 * Generate np G-L-J points (z) and weights (w) for integration over the
 * range [-1, 1] for a polynomial of degree n = np - 1.
 *
 * Reference: Canuto et al., eq (2.3.12).
 * ------------------------------------------------------------------------- */
{
  register int_t i, n;
  real_t         poly, pder, polym1, pderm1, polym2, pderm2, con;
  const real_t   apb = alpha + beta;

  if (np < 2) { z[0] = w[0] = 0.0; return; }

  n    = np - 1;
  con  = pow(2.0,apb + 1.0)*gammaF(alpha + np)*gammaF(beta + np);
  con /= n*gammaF(np)*gammaF(alpha + beta + np + 1.0);
  
  JACGL (n, alpha, beta, z);

  for (i = 0; i < np; i++) {
    JACOBF (n,z[i],alpha,beta,&poly,&pder,&polym1,&pderm1,&polym2,&pderm2);
    w[i] = con / (SQR(poly));
  }
  w[0] *= (beta  + 1.0);
  w[n] *= (alpha + 1.0);
}


real_t PNLEG (const real_t z,
	      const int_t  n)
/* ------------------------------------------------------------------------- *
 * Compute the value of the nth order Legendre polynomial at z, based on the
 * recursion formula for Legendre polynomials.
 * ------------------------------------------------------------------------- */
{
  register int_t k;
  register real_t  dk, p1, p2, p3;
 
  if (n < 1) return 1.0;

  p1 = 1.0;
  p3 = p2 = z;

  for (k = 1; k < n; k++) {
    dk = (real_t) k;
    p3 = ((2.0*dk + 1.0)*z*p2 - dk*p1) / (dk + 1.0);
    p1 = p2;
    p2 = p3;
  }
 
  return p3;
}


real_t PNDLEG (const real_t z,
	       const int_t  n)
/* ------------------------------------------------------------------------- *
 * Compute the value of the derivative of the nth order Legendre polynomial
 * at z, based on the recursion formula for Legendre polynomials.
 * ------------------------------------------------------------------------- */
{
  register int_t  k;
  register real_t dk, p1, p2, p3, p1d, p2d, p3d;

  if (n < 1) return 0.0;

  p2  = z;
  p1d = 0.0;
  p1  = p2d = p3d = 1.0;

  for (k = 1; k < n; k++) {
    dk  = (real_t) k;
    p3  = ((2.0*dk + 1.0)*z*p2 - dk*p1) / (dk + 1.0);
    p3d = ((2.0*dk + 1.0)*p2 + (2.0*dk + 1.0)*z*p2d - dk*p1d) / (dk + 1.0);
    p1  = p2;
    p2  = p3;
    p1d = p2d;
    p2d = p3d;
  }

  return p3d;
}


real_t PND2LEG (const real_t z,
		const int_t  n)
/* ------------------------------------------------------------------------- *
 * Compute the value of the second derivative of the nth order Legendre
 * polynomial at z, based on the definition of the singular Sturm-Liouville
 * problem that generates them (Canuto et al. eq 2.3.1):
 *
 *               (1 - z*z)L(z)'' - 2 z L' + n (n+1) L = 0.
 * ------------------------------------------------------------------------- */
{
  return (2.0 * z * PNDLEG (z, n) - n * (n+1) * PNLEG (z, n) / (1.0 - SQR(z)));
}


void DGLL (const int_t   nz,
	   const real_t* z ,
	   real_t**      D ,
	   real_t**      DT)
/* ------------------------------------------------------------------------- *
 * Compute the derivative operator matrix D and its transpose DT associated
 * with the nth order Lagrangian interpolants through the nz Gauss-Lobatto-
 * Legendre points z.
 *                     dU
 *                     --   = D_ij U_j evaluated at z = z_i
 *                     dz
 *
 * NB: D & DT are both nz x nz.  Canuto et al. eq (2.3.25).
 * ------------------------------------------------------------------------- */
{
  register int_t  i, j, n;
  register real_t d0;

  if (nz < 2) { D[0][0] = DT[0][0] = 0.0; return; }
  
  n  = nz - 1;
  d0 = n * (n + 1) * 0.25;

  D[0][0] = -d0;
  for (i = 1; i < n; i++) D[i][i] = 0.0;
  D[n][n] =  d0;

  for (i = 0; i < nz; i++)
    for (j=0; j<nz; j++) {
      if (i != j) D[i][j] = PNLEG (z[i], n) / (PNLEG(z[j], n) * (z[i] - z[j]));
      DT[j][i] = D[i][j];
    }
}


real_t PNMOD (const real_t z,
	      const int_t  n)
/* ------------------------------------------------------------------------- *
 * Compute the value of the nth order modal basis function at z.
 *
 * p_0(z) = 0.5  (1 + z)
 *
 * p_1(z) = 0.5  (1 - z)
 *                                1,1
 * p_n(z) = 0.25 (1 + z) (1 - z) J   (z)   n >= 2.
 *                                n-2
 *
 * where J is a Jacobi polynomial.
 *
 * Refs:
 *
 * eq.(2.16), R.D. Henderson, "Adaptive Spectral Element Methods", in
 * "High-Order Methods for Computational Physics", eds T.J. Barth &
 * H. Deconinck, Springer, 1999.
 *
 * eq.(2.40), G.E. Karniadakis & S.J. Sherwin, "Spectral/hp Element
 * Methods for CFD", Oxford, 1999.
 * ------------------------------------------------------------------------- */
{
  real_t poly, pder, polym1, pderm1, polym2, pderm2;
 
  if (n  < 1) return 0.5 * (1.0 + z);
  if (n == 1) return 0.5 * (1.0 - z);

  JACOBF (n-2, z, 1.0, 1.0, &poly, &pder, &polym1, &pderm1, &polym2, &pderm2);

  return 0.25 * (1.0 + z) * (1.0 - z) * poly;
}

					     
void uniknot (const int_t nk,
	      real_t*     k )
/* ------------------------------------------------------------------------- *
 * Return nk knot points with uniform spacing on [-1, 1].
 * ------------------------------------------------------------------------- */
{
  register int_t  i, nh;
  register real_t dx;

  if (nk < 2) { *k = 0.0; return; }

  nh = nk >> 1;

  dx      =  2.0 / (real_t) (nk - 1);
  k[0]    = -1.0;
  k[nk-1] =  1.0;
  for (i = 1; i < nh; i++) {
    k[i]      =  k[i-1] + dx;
    k[nk-1-i] = -k[i];
  }
  if (nk & 1) k[i] = 0.0;
}


static real_t gammaF (real_t x)
/* ------------------------------------------------------------------------- *
 * Gamma function for integer or semi-integer values of x.
 * ------------------------------------------------------------------------- */
{
  real_t gamma = 1.0;
  
  if      (x == -0.5) gamma = -2.0*sqrt(M_PI);
  else if (x ==  0.0) return gamma;
  else if ((x-(int_t)x) == 0.5) { 
    int_t  n = (int_t) x;
    real_t tmp = x;

    gamma = sqrt(M_PI);
    while(n--){
      tmp   -= 1.0;
      gamma *= tmp;
    }
  } else if ((x-(int_t)x) == 0.0){
    int_t  n = (int_t) x;
    real_t tmp = x;

    while(--n){
      tmp   -= 1.0;
      gamma *= tmp;
    }
  } else
    fprintf (stderr,"%lf is not of integer or half order\n",x);
  return gamma;
}
