/*****************************************************************************
 * operators.c:  operators for mesh points, derivatives, quadratures.
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
 * All 2D matrices have row-major ordering (but we have moved to flat storage).
 * All routines are real_t (double) precision.
 * For definitions of constants (LL, GLL, etc) see cfemdef.h.
 *
 * $Id: operators.c,v 8.1 2015/04/20 11:14:14 hmb Exp $
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include <cfemdef.h>
#include <cveclib.h>
#include <cfemlib.h>
#include <polylib.h>


typedef struct quadop  { /* ------- quadrature operator information  ------- */
  char           rule  ; /* Quadrature rule: 'G', 'R', or 'L'                */
  int_t          np    ; /* Number of quadrature points                      */
  real_t         alpha ; /* Constants in definition of singular Sturm--      */
  real_t         beta  ; /*   Liouville problem, (Jacobi wt functions) >=-1  */
  real_t*        point ; /* Quadrature points on [-1, 1]                     */
  real_t*        weight; /* Quadrature weights                               */
  real_t*        deriv ; /* Derivative operator for Lagrangian interpolant   */
  real_t*        derivT; /* Transpose of the above (both in row-major order) */
  struct quadop* next  ; /* link to next                                     */
} QuadOp;		 /* ------------------------------------------------ */

static QuadOp* Qroot = 0;

void zquad (const real_t** point , /* Quadrature points.                     */
	    const real_t** weight, /* Quadrature weights.                    */
	    const real_t** dv    , /* Derivative operator at points.         */
	    const real_t** dt    , /* (Transpose) derivative operator.       */
	    const int_t    np    , /* Input: Number of quadrature points.    */
	    const char     rule  , /* Input: 'G'auss, 'R'adau, or 'L'obatto. */
	    const real_t   alpha , /* Input: Jacobi polynomial constant.     */
	    const real_t   beta  ) /* Input: Jacobi polynomial constant.     */
/* ------------------------------------------------------------------------- *
 * Maintain/return QUADRATURE operators for finite elements with
 * spectral basis functions, defined on the master interval [-1, +1].
 *
 * If the required operators are "in stock", they are returned from a list,
 * otherwise they are created and added to the list as a first step.
 * NULL input pointers are not assigned a value.
 *
 * NB: Gauss-Radau integration rules are asymmetric (they use a point
 * at one end of the interval, -1 or +1, as well as interior
 * points). Here the +1 endpoint is assumed.
 * ------------------------------------------------------------------------- */
{
  char    routine[] = "quad";
  int_t   found;
  QuadOp* p;

  for (found = 0, p = Qroot; p; p=p->next) {
    found = 
      p->rule  == rule  &&
      p->np    == np    &&
      p->alpha == alpha && 
      p->beta  == beta;
    if (found) break;
  }

  if (!found) {			/* -- Make new storage area, fill it. */

    p = (QuadOp *) calloc (1, sizeof (QuadOp));
    if (Qroot) p -> next = Qroot; Qroot = p;

    p->rule     = rule;
    p->np       = np;
    p->alpha    = alpha;
    p->beta     = beta;
    p -> point  = (real_t*) malloc (sizeof (real_t) * np);
    p -> weight = (real_t*) malloc (sizeof (real_t) * np);
    p -> deriv  = (real_t*) malloc (sizeof (real_t) * np*np);
    p -> derivT = (real_t*) malloc (sizeof (real_t) * np*np);

    if        (rule == 'L') {	/* Gauss-Lobatto-Jacobi -- most common case.*/
      zwglj   (p->point, p->weight, np, alpha, beta);
      Dglj    (p->deriv, p->derivT, p->point, np, alpha, beta);
    } else if (rule == 'R') {	/* Gauss-Radau-Jacobi   */
      zwgrjp  (p->point, p->weight, np, alpha, beta);
      Dgrjp   (p->deriv, p->derivT, p->point, np, alpha, beta);
    } else if (rule == 'G') {	/* Gauss-Jacobi         */
      zwgj    (p->point, p->weight, np, alpha, beta);
      Dgj     (p->deriv, p->derivT, p->point, np, alpha, beta);
    } else {
      char err[STR_MAX];
      sprintf (err, "unrecognized quadrature rule: %c", rule);
      message (routine, err, ERROR);
    }
  }

  /* -- p now points to valid storage: return requested operators. */

  if (point)  *point  = (const real_t*) p -> point;
  if (weight) *weight = (const real_t*) p -> weight;
  if (dv)     *dv     = (const real_t*) p -> deriv;
  if (dt)     *dt     = (const real_t*) p -> derivT;
}


typedef struct projop { /* ------- projection operator information  -------- */
  int_t          nfrom; /* Number of points projection is "from"             */
  char           rfrom; /* Mesh point definition: 'G', 'R', 'L', or 'U'      */
  real_t         afrom; /* Jacobi weight function powers of "from" mesh      */
  real_t         bfrom; /*   ignored for 'U' (uniform) mesh spacing          */
  int_t          nto  ; /* Number of points projection is "to".              */
  char           rto  ; /* Quadrature rule: 'G', 'R', or 'L'                 */
  real_t         ato  ; /* As for "from" mesh definitions.                   */
  real_t         bto  ; /*                                                   */
  real_t*        IN   ; /* Interpolant/projection matrix, row-major order    */
  real_t*        IT   ; /* Transpose of IN                                   */
  struct projop* next ; /* link to next                                      */
} ProjOp;		/* ------------------------------------------------- */

static ProjOp* Proot = 0;

void proj (const real_t** IN    , /* Interpolant operator matrix             */
	   const real_t** IT    , /* Transposed interpolant operator matrix  */
	   const int_t    nfrom , /* Input: Number of "from" points          */
	   const char     rulefr, /* Input: 'G', 'R', 'L', or 'U', "from"    */
	   const real_t   alphfr, /* Input: Jacobi constant, "from"          */
	   const real_t   betafr, /* Input: Jacobi constant, "from"          */
	   const int_t    nto   , /* Input: Number of "to" points            */
	   const char     ruleto, /* Input: 'G', 'R', 'L', or 'U', "to"      */
	   const real_t   alphto, /* Input: Jacobi constant, "to"            */
	   const real_t   betato) /* Input: Jacobi constant, "to"            */
/* ------------------------------------------------------------------------- *
 * Maintain/return operator matrices for projection *from* one mesh
 * *to* another.  Spectral basis functions, defined on the master
 * interval [-1, +1]. If the "from" or "to" meshes are "spectral"
 * (i.e. 'G', 'R', or 'L'), they are assumed to be defined by the
 * generating Gauss-Jacobi-type rule, and the alpha, beta
 * Jacobi-constant pair. Uniformly-spaced meshes are denoted by 'U',
 * in which case the constants alpha, beta, are ignored.
 *
 * If the required operators are "in stock", they are returned from a list,
 * otherwise they are created and added to the list as a first step.
 * NULL input pointers are not assigned a value.
 *
 * NB: Gauss-Radau integration rules are asymmetric (they use a point
 * at one end of the interval, -1 or +1, as well as interior
 * points). Here the +1 endpoint is assumed.
 * ------------------------------------------------------------------------- */
{
  char    routine[] = "proj";
  int_t   rules, constsfr, conststo;
  int_t   found;
  ProjOp* p;

  for (found = 0, p = Proot; p; p=p->next) {
    rules    = p->nfrom == nfrom && p->rfrom == rulefr &&
               p->nto   == nto   && p->rto   == ruleto ;
    constsfr = (rulefr == 'U') ? 1 : (p->afrom== alphfr && p->bfrom == betafr);
    conststo = (ruleto == 'U') ? 1 : (p->ato  == alphto && p->bto   == betato);
    found    = rules && constsfr && conststo;
    if (found) break;
  }

  if (!found) {			/* -- Make new storage area, fill it. */

    real_t  *zfrom, *zto, *wfrom, *wto, **IN, **IT;

    p = (ProjOp *) calloc (1, sizeof (ProjOp));
    if (Proot) p -> next = Proot; Proot = p;

    p->nfrom = nfrom;
    p->nto   = nto;
    p->rfrom = rulefr;
    p->rto   = ruleto;
    p->afrom = (p->rfrom == 'U') ? 0.0 : alphfr;
    p->bfrom = (p->rfrom == 'U') ? 0.0 : betafr;
    p->ato   = (p->rto   == 'U') ? 0.0 : alphto;
    p->bto   = (p->rto   == 'U') ? 0.0 : betato;
    p -> IN  = dvector (0, nto*nfrom-1);
    p -> IT  = dvector (0, nfrom*nto-1);

    zfrom = dvector (0, nfrom-1);
    wfrom = dvector (0, nfrom-1);
    zto   = dvector (0, nto-1);
    wto   = dvector (0, nto-1);
    IN    = dmatrix (0, nto-1,   0, nfrom-1);
    IT    = dmatrix (0, nfrom-1, 0, nto-1);

    if      (rulefr == 'G') zwgj    (zfrom, wfrom, nfrom, alphfr, betafr);
    else if (rulefr == 'R') zwgrjp  (zfrom, wfrom, nfrom, alphfr, betafr);
    else if (rulefr == 'L') zwglj   (zfrom, wfrom, nfrom, alphfr, betafr);
    else if (rulefr == 'U') uniknot (nfrom, zfrom);
    else {
      char err[STR_MAX];
      sprintf (err, "unrecognized input mesh rule: %c", rulefr);
      message (routine, err, ERROR);
    }

    if      (ruleto == 'G') zwgj    (zto, wto, nto, alphto, betato);
    else if (ruleto == 'R') zwgrjp  (zto, wto, nto, alphto, betato);
    else if (ruleto == 'L') zwglj   (zto, wto, nto, alphto, betato);
    else if (ruleto == 'U') uniknot (nto, zto);
    else {
      char err[STR_MAX];
      sprintf (err, "unrecognized output mesh rule: %c", ruleto);
      message (routine, err, ERROR);
    }

    intmat_g (nfrom, zfrom, nto, zto, *IN, *IT);
  
    dcopy (nto*nfrom, *IN, 1, p ->IN, 1);
    dcopy (nfrom*nto, *IT, 1, p ->IT, 1);

    freeDvector (zfrom, 0);
    freeDvector (wfrom, 0);
    freeDvector (zto,   0);
    freeDvector (wto,   0);
    freeDmatrix (IN, 0, 0);
    freeDmatrix (IT, 0, 0);
  }

  /* -- p now points to valid storage: return requested operators. */

  if (IN) *IN = (const real_t*) p->IN;
  if (IT) *IT = (const real_t*) p->IT;
}


void intp (real_t*      inr   ,  /* 1D shape function at r                   */
	   real_t*      ins   ,  /* 1D shape function at s                   */
	   real_t*      dvr   ,  /* 1D shape function derivative at r        */
	   real_t*      dvs   ,  /* 1D shape function derivative at s        */
	   const int_t  nr    ,
	   const char   ruler ,
	   const real_t alphar,
	   const real_t betar ,
	   const int_t  ns    ,
	   const char   rules ,
	   const real_t alphas,
	   const real_t betas ,
	   const real_t r     ,  /* location of r in [-1, 1]                 */
	   const real_t s     )  /* location of s in [-1, 1]                 */
/* ------------------------------------------------------------------------- *
 * Return the interpolation/derivative vectors for tensor-product shape
 * functions evaluated at canonical coordinates (r, s), each in [-1. 1].
 *
 * Take no action for empty vectors.
 * ------------------------------------------------------------------------- */
{
  const real_t *kr, *ks;
  real_t       x[1], *DVr, *DTr, *DVs, *DTs;

  /* -- NB: we should really pass this in as a work array. */

  DVr = malloc (2*nr*ns * sizeof (real_t));
  DTr = DVr + nr;
  DVs = DTr + nr;
  DTs = DVs + ns;

  zquad (&kr, NULL, NULL, NULL, nr, ruler, alphar, betar);
  zquad (&ks, NULL, NULL, NULL, ns, rules, alphas, betas);

  x[0] = r;
  if (inr) {
    intmat_g (nr, kr, 1, x, DVr, DTr);
    dcopy    (nr, DVr, 1, inr, 1);
  }
  if (dvr) {
    dermat_g (nr, kr, 1, x, DVr, DTr);
    dcopy    (nr, DVr, 1, dvr, 1);
  }

  x[0] = s;
  if (ins) {
    intmat_g (ns, ks, 1, x, DVs, DTs);
    dcopy    (ns, DVs, 1, ins, 1);
  }
  if (dvs) {
    dermat_g (ns, ks, 1, x, DVs, DTs);
    dcopy    (ns, DVs, 1, dvs, 1);
  }

  free (DVr);
}

 
typedef struct legcoef {	/* ---- Table for GLL Legendre transform --- */
  int_t           np  ;		/* Number of mesh points                     */
  real_t*         dtab;		/* (np+1)*np table of polynomials & weights  */
  struct legcoef* next;		/* link to next one                          */
} legCoef;			/* ----------------------------------------- */

static legCoef* lChead = 0;


void dglldpc (const int_t    np,
	      const real_t** cd)
/* ------------------------------------------------------------------------- *
 * Return pointers to look-up tables of Legendre polynomials and
 * weights, based on the Gauss--Lobatto nodes.  The tables can be used
 * in calculating discrete Legendre transforms.
 *
 * Tables contain values of Legendre polynomials evaluated at GLL
 * quadrature points, and weights.  Storage is row-major: the first np
 * points contain values of the zeroth Legendre polynomial at the GLL
 * points (these are identically 1.0), the second np points the values
 * of the first polynomial, etc.  So rows index polynomial, while
 * columns index spatial location. 
 *
 * NB: this is the opposite ordering used for dglmdpc, below.
 *
 * The final np points (the np-th row of the table, for 0-based
 * indexing) contains the weights 1/gamma_k (see Canuto et al.,
 * 2.3.13).
 * -------------------------------------------------------------------------
 * */
{
  register int_t    found = 0;
  register legCoef* p;

  for (p = lChead; p; p = p->next) {
    found = p -> np == np;
    if (found) break;
  }

  if (!found) {		/* -- Make more storage and operators. */
    register int_t i, j, k;
    const int_t    nm = np - 1;
    const real_t*  z;

    p = (legCoef*) calloc (1, sizeof (legCoef));
    if (lChead) p -> next = lChead;
    lChead = p;

    p -> np   = np;
    p -> dtab = dvector (0, np * (np + 1) - 1);

    zquad (&z, NULL, NULL, NULL, np, 'L', 0.0, 0.0);

    for (i = 0; i < np; i++) {
      k = i + np * np;
      p -> dtab[k] = (i < nm) ?  0.5*(i+i+1) : 0.5*nm;
      for (j = 0; j < np; j++) {
	k = j + i * np;
	p -> dtab[k] = PNLEG (z[j], i);
      }
    }
  }

  /* p now points to valid storage: return requested operators. */

  if (cd) *cd = (const real_t*) p -> dtab;
}


typedef struct legtran {	/* ---- GLL Legendre transform matrices  --- */
  int_t           np  ;		/* Number of mesh points                     */
  real_t*         FW  ;		/* np*np forward Legendre transform matrix.  */
  real_t*         FT  ;		/* Transpose of FW.                          */
  real_t*         BW  ;		/* np*np inverse Legendre transform matrix.  */
  real_t*         BT  ;		/* Transpose of BW.                          */
  real_t*         UF  ;		/* np^2*np^2 forward transform matrix.       */
  real_t*         UB  ;		/* np^2*np^2 inverse transform matrix.       */
  struct legtran* next;		/* link to next one                          */
} legTran;			/* ----------------------------------------- */

static legTran* lThead = 0;


void dglldpt (const int_t    np,
	      const real_t** fw,
	      const real_t** ft,
	      const real_t** bw,
	      const real_t** bt,
	      const real_t** fu,
	      const real_t** bu)
/* ------------------------------------------------------------------------- *
 * Return pointers to matrices for 1D and 2D Legendre polynomial
 * transforms of np data, based on the Gauss--Lobatto nodes.  See the
 * definitions of forward and inverse polynomial transforms in Canuto
 * et al., Sections 2.2.3, 2.2.13.
 *
 * The 1D discrete transform is defined as (summation implied)
 *
 *   Ak = Ck Wj Uk Lkj,
 *
 * with inverse
 *
 *   Uj = Ak Lkj,
 *
 * where Lij is the ith Legendre polynomial evaluated at the jth
 * quadrature point, Ck correspond to 1/\gamma_k in Canuto et al., and
 * Wj are the quadrature weights for the corresponding points.
 *
 * The equivalent 2D tensor-product (or sum-factorisation) relationships are
 *                                            t
 *   Aij = Ci Cj Wp Lip Wq Upq Ljq = Rip Upq Rqj
 *                                    t
 *   Uij = Lpi Apq Lqj             = Lip Apq Lqj
 *
 * where Rip = Ci Wp Lip.
 *
 * Instead of these tensor-product forms, we can also "unroll" the
 * matrix multiplies above to be a single matrix-vector product in
 * each case, so that
 *
 *   Aij = Fijpq Upq,  Uij = Bijpq Apq.
 *
 * This routine supplies the transform matrices.  If any pointer is
 * passed in as NULL, we don't return the corresponding value.
 * Matrices are supplied 1D, with row-major ordering.
 * ------------------------------------------------------------------------- */
{
  register int_t    found = 0;
  register legTran* p;

  for (p = lThead; p; p = p->next) {
    found = p -> np == np;
    if (found) break;
  }

  if (!found) {		/* -- Make more storage and operators. */
    register int_t i, j, k, l, r, s;
    const int_t    np2 = np * np;
    const real_t   *tab, *w;
    real_t         ci;

    p = (legTran*) calloc (1, sizeof (legTran));
    if (lThead) p -> next = lThead;
    lThead = p;

    p -> np = np;
    p -> FW = dvector (0, 4 * np2 + 2 * np2*np2 - 1);
    p -> FT = p -> FW + np2;
    p -> BW = p -> FT + np2;
    p -> BT = p -> BW + np2;
    p -> UF = p -> BT + np2;
    p -> UB = p -> UF + np2*np2;

    dglldpc  (np, &tab);
    zquad (NULL, &w, NULL, NULL, np, 'L', 0.0, 0.0);

    /* -- Create forward & inverse 1D transform matrices. */

    for (i = 0; i < np; i++) {
      ci = tab[np2 + i];
      for (j = 0; j < np; j++) {
	p->FW[i*np + j] = ci * w[j] * tab[i*np +j];
	p->BW[i*np + j] = tab[j*np+ i];
      }
    }

    /* -- And their transposes. */

    for (i = 0; i < np; i++)
      for (j = 0; j < np; j++) {
	p->FT[i*np + j] = p->FW[j*np + i];
	p->BT[i*np + j] = p->BW[j*np + i];
      }

    /* -- Manufacture 2D forward, inverse DLT matrices. */

    for (k = 0, i = 0; i < np; i++)
      for (j = 0; j < np; j++, k++)
	for (l = 0, r = 0; r < np; r++)
	  for (s = 0; s < np; s++, l++) {
	    p->UF[k*np2 + l] = p->FW[i*np + r] * p->FT[s*np + j];
	    p->UB[k*np2 + l] = p->BW[i*np + r] * p->BT[s*np + j];
	  }
  }

  /* p now points to valid storage: return requested operators. */

  if (fw) *fw = (const real_t*) p -> FW;
  if (ft) *ft = (const real_t*) p -> FT;
  if (bw) *bw = (const real_t*) p -> BW;
  if (bt) *bt = (const real_t*) p -> BT;
  if (fu) *fu = (const real_t*) p -> UF;
  if (bu) *bu = (const real_t*) p -> UB;
}


typedef struct modcoef {	/* -- Table for modal expansion transform -- */
  int_t           np  ;		/* Number of mesh points                     */
  real_t*         dtab;		/* np*np table of basis function values      */
  struct modcoef* next;		/* link to next one                          */
} modCoef;			/* ----------------------------------------- */

static modCoef* mChead = 0;


void dglmdpc (const int_t    np,
	      const real_t** cd)
/* ------------------------------------------------------------------------- *
 * Return pointers to look-up tables of modal expansion functions and
 * weights, based on the Gauss--Lobatto nodes.  The tables can be used
 * in calculating discrete Legendre transforms.
 *
 * Tables contain values of modal expansion functions evaluated at GLL
 * quadrature points, and the associated quadrature weights.  Storage
 * is column-major: the first np points contain values of the all the
 * modal basis function at the first GLL point, the second np points
 * the values of the all the basis functions at the second point, etc.
 * So rows index quadrature point, while columns index polynomial.
 *
 * NB: this ordering is THE OPPOSITE of that used for dglldpc, but it
 * conforms to that used by Karniadakis & Sherwin.
 *
 * The modal expansion functions are a set of hierarchical functions
 * based on the Jacobi polynomials.
 *
 * p_0(z) = 0.5  (1 + z)
 *
 * p_1(z) = 0.5  (1 - z)
 *                                1,1
 * p_n(z) = 0.25 (1 + z) (1 - z) J   (z)   n >= 2.
 *                                n-2
 *
 * where J is a Jacobi polynomial.
 * ------------------------------------------------------------------------- */
{
  register int_t    found = 0;
  register modCoef* p;

  for (p = mChead; p; p = p->next) {
    found = p -> np == np;
    if (found) break;
  }

  if (!found) {		/* -- Make more storage and operators. */
    register int_t i, j, k;
    const real_t   *z;

    p = (modCoef*) calloc (1, sizeof (modCoef));
    if (mChead) p -> next = mChead;
    mChead = p;

    p -> np   = np;
    p -> dtab = dvector (0, np * np - 1);

    zquad (&z, NULL, NULL, NULL, np, 'L', 0.0, 0.0);

    for (i = 0; i < np; i++)
      for (j = 0; j < np; j++) {
	k = i + j * np;
	p -> dtab[k] = PNMOD (z[j], i);
      }
  }

  /* -- p now points to valid storage: return requested operators. */

  if (cd) *cd = (const real_t*) p -> dtab;
}


typedef struct modtran {	/* - GL modal expansion transform matrices - */
  int_t           np  ;		/* Number of mesh points                     */
  real_t*         FW  ;		/* np*np forward modal transform matrix.     */
  real_t*         FT  ;		/* Transpose of FW.                          */
  real_t*         BW  ;		/* np*np inverse modal transform matrix.     */
  real_t*         BT  ;		/* Transpose of BW.                          */
  real_t*         UF  ;		/* np^2*np^2 forward transform matrix.       */
  real_t*         UB  ;		/* np^2*np^2 inverse transform matrix.       */
  struct modtran* next;		/* link to next one                          */
} modTran;			/* ----------------------------------------- */

static modTran* mThead = 0;


void dglmdpt (const int_t    np,
	      const real_t** fw,
	      const real_t** ft,
	      const real_t** bw,
	      const real_t** bt,
	      const real_t** fu,
	      const real_t** bu)
/* ------------------------------------------------------------------------- *
 * Return pointers to matrices for 1D and 2D modal expansion
 * transforms of np data, based on the Gauss--Lobatto nodes.  The
 * "modal" expansion functions are a set of hierarchical basis
 * functions closely associated with the GLL Lagrange basis. The
 * transform equations are just the standard "Normal form" projections
 * of one basis onto another.  (For a Legendre tranform -- see dglldpt
 * -- the equations are simplified by the discrete orthogonality of
 * all the basis functions.)  Here, not all (but most) of the basis
 * functions are orthogonal.  The treatment and intended use of the
 * transform matrices is similar to that for dglldpt.
 *
 * The 1D discrete forward transform is defined as
 *
 *   ^     t     -1 t                                  ^
 *   u = [B W B ]  B  W u  = F u,  with inverse  u = B u
 *   ~                  ~      ~                 ~     ~
 *
 * where u_i is a vector of function values evaluated at the ith
 * quadrature point, B_ij is the jth basis function evaluated at the
 * ith quadrature point and W_i are the associated quadrature weights.
 * The vector of coefficients of the new basis functions is u^.
 *
 * The equivalent 2D tensor-product (sum-factorisation) relationships are
 *
 *   ^              t                 ^    t
 *   uij = Fip upq Fqj,     uij = Bip upq Bqj
 *
 * Instead of these tensor-product forms, we can also "unroll" the
 * matrix multiplies above to be a single matrix-vector product in
 * each case, so that
 *   ^                                  ^
 *   uij = Fijpq upq,       uij = Bijpq upq.
 *
 * This form, although involving more operations, can actually be
 * faster on vector architectures (depending on np).
 *
 * This routine supplies the transform matrices.  If any pointer is
 * passed in as NULL, we don't return the corresponding value.
 * Matrices are supplied 1D, with row-major ordering.
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
  const char routine[] = "dglmdpt";

  register int_t    found = 0;
  register modTran* p;

  for (p = mThead; p; p = p->next) {
    found = p -> np == np;
    if (found) break;
  }

  if (!found) {		/* -- Make more storage and operators. */
    register int_t i, j, k, l, r, s;
    const int_t    np2 = np * np;
    const real_t   *B, *W;
    real_t         *work, *BtW, *BtWB, *rwrk;
    int_t          *iwrk, info;

    p = (modTran*) calloc (1, sizeof (modTran));
    if (mThead) p -> next = mThead;
    mThead = p;

    p -> np = np;
    p -> FW = dvector (0, 4 * np2 + 2 * np2*np2 - 1);
    p -> FT = p -> FW + np2;
    p -> BW = p -> FT + np2;
    p -> BT = p -> BW + np2;
    p -> UF = p -> BT + np2;
    p -> UB = p -> UF + np2*np2;

    dglmdpc  (np, &B);
    zquad (NULL, &W, NULL, NULL, np, 'L', 0.0, 0.0);

    iwrk = ivector (0, np - 1);
    work = dvector (0, 4 * np2 - 1);
    BtW  = work + np2;
    BtWB = BtW  + np2;
    rwrk = BtWB + np2;

    /* -- Create matrices BtW & (symmetric) BtWB. */
    
    for (i = 0; i < np; i++)
      dsmul (np, W[i], B + i * np, 1, BtW + i, np);

    dmxm (BtW, np, (real_t*) B, np, BtWB, np);

    /* -- Invert BtWB. */

    dgetrf (np, np, BtWB, np, iwrk, &info);
    if (info) message (routine, "matrix BtWB has singular factor", ERROR);

    dgetri (np, BtWB, np, iwrk, rwrk, 2*np2, &info);
    if (info) message (routine, "matrix BtWB is singular", ERROR);

    /* -- Create forward (FW) & inverse (BW) 1D transform matrices. */

    dmxm  (BtWB, np, BtW, np, p->FW, np);
    dcopy (np2, B, 1,         p->BW,  1);

    /* -- And their transposes. */

    for (i = 0; i < np; i++)
      for (j = 0; j < np; j++) {
	p->FT[i*np + j] = p->FW[j*np + i];
	p->BT[i*np + j] = p->BW[j*np + i];
      }

    /* -- Manufacture 2D forward, inverse DPT matrices. */

    for (k = 0, i = 0; i < np; i++)
      for (j = 0; j < np; j++, k++)
	for (l = 0, r = 0; r < np; r++)
	  for (s = 0; s < np; s++, l++) {
	    p->UF[k*np2 + l] = p->FW[i*np + r] * p->FT[s*np + j];
	    p->UB[k*np2 + l] = p->BW[i*np + r] * p->BT[s*np + j];
	  }

    freeIvector (iwrk, 0);
    freeDvector (work, 0);
  }

  /* p now points to valid storage: return requested operators. */

  if (fw) *fw = (const real_t*) p -> FW;
  if (ft) *ft = (const real_t*) p -> FT;
  if (bw) *bw = (const real_t*) p -> BW;
  if (bt) *bt = (const real_t*) p -> BT;
  if (fu) *fu = (const real_t*) p -> UF;
  if (bu) *bu = (const real_t*) p -> UB;
}
