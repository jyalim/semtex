/*****************************************************************************
 *                  3^i TENSOR & VECTOR CALCULATIONS
 *
 * Copyright (c) 1992 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn

 * NB: In all the following calculations, the 3x3 tensor which is used as
 * input is supplied as a 1-D array of double. It is assumed that the tensor
 * is supplied in column-major order (as it would be stored in FORTRAN):
 *
 *           --          --               --         --
 * TENSORij: | 11  12  13 |      STORAGE: | 0   3   6 |
 * i=row     | 21  22  23 |               | 1   4   7 |
 * j-col     | 31  32  33 |               | 2   5   8 |
 *           --          --               --         --
 *
 * There is a code in the naming of the tensors used below: if the operation
 * can be performed on a general tensor, then the tensor argument is T; if
 * the tensor is assumed symmetric, it is called S, A if it is antisymmetric
 * and if the operation has special significance for the velocity gradient
 * tensor, it is called VG.
 *
 * --
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
 *****************************************************************************/

static char
RCSid[] = "$Id: tensorcalcs.c,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <stdio.h>	
#include <math.h>
#include "tensorcalcs.h"

#define SWAP(a, b) {double temp = (a); (a) = (b); (b) = temp;}
#define MAX(a, b)  ((a) > (b) ? (a) : (b))
#define SGN(a)     ((a) >= 0.0 ? 1.0 : -1.0) 
#define SQR(a)     ((a) * (a))
#define CUBE(a)    ((a) * (a) * (a))
#define ATHIRD     0.333333333333333
#define SQRT3      1.73205080757
#define C1        -1.88988157484


void transpose(const double T[9], double O[9])
/* ------------------------------------------------------------------------- *
 * Return the transpose of 3x3 tensor T in O.
 * ------------------------------------------------------------------------- */
{
  O[0] = T[0]; O[4] = T[4]; O[8] = T[8];
  O[3] = T[1]; O[1] = T[3];
  O[6] = T[2]; O[2] = T[6];
  O[7] = T[5]; O[5] = T[7];
}
  

void transpose_in_place(double T[9])
/* ------------------------------------------------------------------------- *
 * C'est veut dire...
 * ------------------------------------------------------------------------- */
{
  SWAP(T[3], T[1]);
  SWAP(T[6], T[2]);
  SWAP(T[7], T[5]);
}

 
void symmetric_part(const double T[9], double S[9])
/* ------------------------------------------------------------------------- *
 * Return the symmetric part of tensor T in S.
 * ------------------------------------------------------------------------- */
{
  S[0] = T[0];
  S[4] = T[4];
  S[8] = T[8];
  S[3] = (S[1] = 0.5 * (T[1] + T[3]));
  S[6] = (S[2] = 0.5 * (T[2] + T[6]));
  S[7] = (S[5] = 0.5 * (T[5] + T[7]));
}


void anti_symmetric_part(const double T[9], double A[9])
/* ------------------------------------------------------------------------- *
 * Return antisymmetric part of T in A.
 * ------------------------------------------------------------------------- */
{
  A[0] =  (A[4] = (A[8] = 0.0));
  A[1] = -(A[3] = 0.5 * (T[3] - T[1]));
  A[2] = -(A[6] = 0.5 * (T[6] - T[2]));
  A[5] = -(A[7] = 0.5 * (T[7] - T[5]));
}


double contracted_product(const double T[9])
/* ------------------------------------------------------------------------- *
 * Calculate TijTji for a 3x3 tensor [=trace(T^2)].
 * ------------------------------------------------------------------------- */
{
  return(T[0]*T[0]+T[4]*T[4]+T[8]*T[8] + 2.0*(T[1]*T[3]+T[2]*T[6]+T[5]*T[7]));
}


void S2plusA2(const double V[9], double O[9])
/* ------------------------------------------------------------------------- *
 * Return in O the sum of the squares of the symmetric and
 * antisymmetric parts of V. The outcome is a symmetric ternsor.
 * ------------------------------------------------------------------------- */
{
  double W[9];

  symmetric_part(V, W);
  O[0] = SQR(W[0]) + SQR(W[3]) + SQR(W[6]);
  O[4] = SQR(W[1]) + SQR(W[4]) + SQR(W[7]);
  O[8] = SQR(W[2]) + SQR(W[5]) + SQR(W[8]);
  O[3] = O[1] = W[3]*(W[0]+W[4]) + W[5]*W[6];
  O[6] = O[2] = W[6]*(W[0]+W[8]) + W[3]*W[7];
  O[7] = O[5] = W[7]*(W[4]+W[8]) + W[1]*W[6];

  anti_symmetric_part(V, W);
  O[0] -= (SQR(W[3]) + SQR(W[6]));
  O[4] -= (SQR(W[3]) + SQR(W[7]));
  O[8] -= (SQR(W[6]) + SQR(W[7]));
  W[0]  = W[6]*W[5];
  W[4]  = W[3]*W[7];
  W[8]  = W[1]*W[6];
  O[3] += W[0]; O[1] += W[0];
  O[6] += W[4]; O[2] += W[4];
  O[7] += W[8]; O[5] += W[8];
}


double enstrophy(const double VG[9])
/* ------------------------------------------------------------------------- *
 * Calculate RijRji (enstrophy) from velocity gradient data.
 * ------------------------------------------------------------------------- */
{
  double A[9];

  anti_symmetric_part(VG, A);
  return(-2.0*(A[3]*A[1] + A[6]*A[2] + A[7]*A[5]));
}


double dissipation(const double VG[9])
/* ------------------------------------------------------------------------- *
 * Calculate SijSji (proportional to viscous dissipation) from vel grad.
 * ------------------------------------------------------------------------- */
{
  double S[9];

  symmetric_part(VG, S);
  return(contracted_product(S));
}


double strainrate(const double VG[9])
/* ------------------------------------------------------------------------- *
 * Calculate strain rate magnitude sqrt(2SijSji) from vel grad.
 * ------------------------------------------------------------------------- */
{
  double S[9];

  symmetric_part(VG, S);
  return(sqrt(2.0*contracted_product(S)));
}


double trace(const double T[9])
/* ------------------------------------------------------------------------- *
 * Return trace of 3x3 tensor.
 * ------------------------------------------------------------------------- */
{
  return(T[0] + T[4] + T[8]);
}


void remove_trace(double T[9])
/* ------------------------------------------------------------------------- *
 * Subtract Tii/3 from leading diagonal.
 * ------------------------------------------------------------------------- */
{
  double tr;

  tr    = ATHIRD * trace(T);
  T[0] -= tr;
  T[4] -= tr;
  T[8] -= tr;
}


double det(const double T[9])
/* ------------------------------------------------------------------------- *
 * Return determinant of 3x3 tensor.
 * ------------------------------------------------------------------------- */
{
  double d;

  d  = T[0] * (T[4]*T[8] - T[7]*T[5]);
  d -= T[3] * (T[1]*T[8] - T[7]*T[2]);
  d += T[6] * (T[1]*T[5] - T[4]*T[2]);

  return(d);
}


void invariants(const double T[9], double *I, double *II, double *III)
/* ------------------------------------------------------------------------- *
 * Calculate the invariants of a 3x3 tensor: the eigenvalues of T satisfy
 * x^3 + Ix^2 + IIx + III = 0.  See Aris Eq. (2.5.5).
 * ------------------------------------------------------------------------- */
{
  *I   = -trace(T);
  *II  = T[4]*T[8] - T[7]*T[5] + T[8]*T[0] - T[2]*T[6] + T[0]*T[4] - T[3]*T[1];
  *III = -det(T);
}


void real_eigenvalues(double  a2, double  a1, double  a0,
		      double *z1, double *z2, double *z3)
/* ------------------------------------------------------------------------- *
 * More generally, the real roots of a cubic (roots are guaranteed to be all
 * real if the coefficients ai are the invariants of a symmetric 3x3 matrix).
 * We find the roots of x^3 + a2.x^2 + a1.x + a0 = 0 and return them sorted
 * in descending order.  See Abramowitz & Stegun 3.8.
 * ------------------------------------------------------------------------- */
{
  double q, r, t, thetaon3, cosA, sinA;

  q = a1 / 3.0 - SQR(a2) / 9.0;
  r = (a1*a2 - 3.0*a0) / 6.0 - CUBE(a2) / 27.0;

  t = CUBE(q) + SQR(r);

  if (t > 0) {
#if 0 /* -- This is mostly caused by rounding errors. Ignore. */
    fprintf(stderr, "warning: t = %g: ==> complex conjugate roots\n", t);
#endif
    t = 0.0;
  }

  t = sqrt(-t);
  thetaon3 = atan2(t, r)/3.0;
  cosA = cos(thetaon3);
  sinA = sin(thetaon3);
  q = sqrt(-q);

  *z1 = 2.0*q*cosA - a2/3.0;
  *z2 = (*z3 = -q*cosA -a2/3.0);
  *z2 = *z2 - SQRT3*q*sinA;
  *z3 = *z3 + SQRT3*q*sinA;

  if (*z2 > *z1) SWAP(*z1, *z2);
  if (*z3 > *z2) SWAP(*z2, *z3);
  if (*z2 > *z1) SWAP(*z1, *z2);
}


void eigenvector(const double T[9], double z, double e[3])
/* ------------------------------------------------------------------------- *
 * Given z, a real eigenvalue of T 3x3, return the corresponding
 * eigenvector e, normalized with length = 1.
 * ------------------------------------------------------------------------- */
{
  e[1] = ((T[0] - z)*T[7] - T[1]*T[6]) / ((T[4] - z)*T[6] - T[3]*T[7]);
  e[2] = ((T[0] - z)*T[5] - T[2]*T[3]) / ((T[8] - z)*T[3] - T[6]*T[5]);
  e[0] = 1.0 / sqrt(1.0 + SQR(e[1]) + SQR(e[2]));
  e[1] *= e[0];
  e[2] *= e[0];
}


void vec(const double A[9], double V[3])
/* ------------------------------------------------------------------------- *
 * Return the vector equivalent of an antisymmetric 3x3 tensor A.
 * If A is the antisymmetric part of the velocity gradient tensor, the local 
 * angular velocity is -vec(A) and the vorticity is -2vec(A).
 * ------------------------------------------------------------------------- */
{
  V[0] = A[7];
  V[1] = A[2];
  V[2] = A[3];
}


void scale_vect(double v[3], double s)
/* ------------------------------------------------------------------------- *
 * Scale the components of vector v by factor s.
 * ------------------------------------------------------------------------- */
{
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}


void vorticity (const double VG[9], double w[3])
/* ------------------------------------------------------------------------- *
 * From A, antisymmetric part of VG, return the vorticity vector.
 * ------------------------------------------------------------------------- */
{
  double A[9];

  anti_symmetric_part(VG, A);
  vec(A, w);
  scale_vect(w, -2.0);
}


double dot(const double a[3], const double b[3])
/* ------------------------------------------------------------------------- *
 * Dot product of two vectors.
 * ------------------------------------------------------------------------- */
{
  return(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}


void cross(const double a[3], const double b[3], double c[3])
/* ------------------------------------------------------------------------- *
 * Return the cross product of two vectors.
 * ------------------------------------------------------------------------- */
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}


void normalize(double a[3])
/* ------------------------------------------------------------------------- *
 * Scale vector a to have length 1.
 * ------------------------------------------------------------------------- */
{
  double l;

  l = sqrt(SQR(a[0]) + SQR(a[1]) + SQR(a[2]));
  scale_vect(a, 1.0/l);
}


double classify_incompressible(const double VG[9], int *topology)
/* ------------------------------------------------------------------------- *
 * Classify velocity gradient tensor assumed trace-free as one of
 *     UFC (unstable focus/compressing  <--> blue)
 *     UNS (unstable node/saddle/saddle <--> green)
 *     SFS (stable focus/stretching     <--> yellow)
 *     SNS (stable node/saddle/saddle   <--> red)
 * The scalar magnitude given to each is the distance from the origin in q-r
 * space---a mapping of Q-R space in which the Q-R coordinates are first
 * mapped to a similarity coordinate system in which the topology-splitting
 * curve in the P=0 plane becomes invariant with change of scale, and then
 * sheared to produce quadrants.
 * ------------------------------------------------------------------------- */
{
  double P, Q, R, q, r;

  invariants(VG, &P, &Q, &R);

  r = SGN(R) * pow((SQR(R)), ATHIRD);
  q = Q - C1*fabs(r);

  if (r > 0.0) {
    if (q > 0.0)
      *topology = UFC;
    else
      *topology = UNS;
  } else {
    if (q > 0.0)
      *topology = SFS;
    else
      *topology = SNS;
  }

  return(sqrt(SQR(q) + SQR(r)));
}


double lambda2(const double VG[9])
/* ------------------------------------------------------------------------- *
 * Return the measure lambda2, the second eigenvalue of the sum of squares
 * of the symmetric and antisymmetric parts of VG, used by Jeong & Hussain
 * JFM 285 to identify vortex cores.
 * ------------------------------------------------------------------------- */
{
  double SS[9];
  double I, II, III, e1, e2, e3;

  S2plusA2(VG, SS);
  invariants(SS, &I, &II, &III);
  real_eigenvalues(I, II, III, &e1, &e2, &e3);

  return(e2);
}


double discrimi(const double VG[9])
/* ------------------------------------------------------------------------- *
 * Return the discriminant of VG, assuming the velocity field incompressible.
 * See Blackburn et at, JFM 310.
 * ------------------------------------------------------------------------- */
{
  double I, II, III, D;

  invariants(VG, &I, &II, &III);

  /* Assume (don't check) that I = 0. */

  D = 6.75*SQR(II) + CUBE(III);

  return(D);
}


void lambvector(const double VG[9], const double vel[3], double lamb[3])
/* ------------------------------------------------------------------------- *
 * Return the Lamb vector vorticity x velocity.
 * See Hamman et al, JFM 610.
 * ------------------------------------------------------------------------- */
{
  double vort[3];

  vorticity (VG, vort);
  
  lamb[0] = vort[1]*vel[2] - vort[2]*vel[1];
  lamb[1] = vort[2]*vel[0] - vort[0]*vel[2];
  lamb[2] = vort[0]*vel[1] - vort[1]*vel[0];
}


double helicity(const double VG[9], const double vel[3])
/* ------------------------------------------------------------------------- *
 * Return helicity 0.5*(vorticity . velocity).
 * ------------------------------------------------------------------------- */
{
  double vort[3];

  vorticity (VG, vort);
  
  return(0.5*(vort[0]*vel[0] + vort[1]*vel[1] + vort[2]*vel[2]));
}
