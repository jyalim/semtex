#ifndef C_FEMLIB_H
#define C_FEMLIB_H
/*****************************************************************************
 *         FUNCTION PROTOTYPES FOR ROUTINES IN LIBRARY LIBFEM.A
 *
 * $Id: cfemlib.h,v 8.1 2015/04/20 11:14:14 hmb Exp $
 *****************************************************************************/

#include "cfemdef.h"

/* -- Routines from initial.c: */

void    yy_initialize (void);
real_t  yy_interpret  (const char*);

void    yy_vec_init   (const char*, const char*);
void    yy_vec_interp (const int_t, ...);

void    yy_help       (void);
void    yy_show       (void);

/* -- Routines from polyops.c: */

void   dermat_g (const int_t, const real_t*, const int_t,
		 const real_t*, real_t*, real_t*);
void   intmat_g (const int_t, const real_t*, const int_t,
		 const real_t*, real_t*, real_t*);
void   dermat_k (const int_t, const real_t*, real_t*, real_t*);

void   JACG     (const int_t, const real_t, const real_t, real_t*);
void   JACGR    (const int_t, const real_t, const real_t, real_t*);
void   JACGL    (const int_t, const real_t, const real_t, real_t*);

void   ZWGL     (real_t*, real_t*, const int_t);
void   ZWGRL    (real_t*, real_t*, const int_t);
void   ZWGLL    (real_t*, real_t*, const int_t);
void   ZWGLJ    (real_t*, real_t*, const real_t, const real_t, const int_t);

void   DGLL     (const int_t, const real_t*, real_t**, real_t**);

real_t PNLEG    (const real_t, const int_t);
real_t PNDLEG   (const real_t, const int_t);
real_t PND2LEG  (const real_t, const int_t);
real_t PNMOD    (const real_t, const int_t);

void   uniknot  (const int_t, real_t*);

/* -- Routines from operators.c: */

void zquad (const real_t** point ,  /* Quadrature points.                    */
	    const real_t** weight,  /* Quadrature weights.                   */
	    const real_t** dv    ,  /* Derivative operator at points.        */
	    const real_t** dt    ,  /* (Transpose) derivative operator.      */
	    const int_t    np    ,  /* Input: Number of quadrature points.   */
	    const char     rule  ,  /* Input: 'G'auss, 'R'adau, or 'L'obatto.*/
	    const real_t   alpha ,  /* Input: Jacobi polynomial constant.    */
	    const real_t   beta  ); /* Input: Jacobi polynomial constant.    */

void proj  (const real_t** IN    ,  /* Interpolant operator matrix           */
	    const real_t** IT    ,  /* Transposed interpolant operator matrix*/
	    const int_t    nfrom ,  /* Input: Number of "from" points        */
	    const char     rulefr,  /* Input: 'G', 'R', 'L', or 'U', "from"  */
	    const real_t   alphfr,  /* Input: Jacobi constant, "from"        */
	    const real_t   betafr,  /* Input: Jacobi constant, "from"        */
	    const int_t    nto   ,  /* Input: Number of "to" points          */
	    const char     ruleto,  /* Input: 'G', 'R', 'L', or 'U', "to"    */
	    const real_t   alphto,  /* Input: Jacobi constant, "to"          */
	    const real_t   betato); /* Input: Jacobi constant, "to"          */

void intp  (real_t*        inr   ,  /* 1D shape function at r                */
	    real_t*        ins   ,  /* 1D shape function at s                */
	    real_t*        dvr   ,  /* 1D shape function derivative at r     */
	    real_t*        dvs   ,  /* 1D shape function derivative at s     */
	    const int_t    nr    ,  /* Input: number of points, r-direction  */
	    const char     basisr,  /* Input: 'G', 'R', 'L', r-direction     */
	    const real_t   alphar,  /* Input: Jacobi constant, r-direction   */
	    const real_t   betar ,  /* Input: Jacobi constant, r-direction   */
	    const int_t    ns    ,  /* Input: number of points, s-direction  */
	    const char     basiss,  /* Input: 'G', 'R', 'L', s-direction     */
	    const real_t   alphas,  /* Input: Jacobi constant, s-direction   */
	    const real_t   betas ,  /* Input: Jacobi constant, s-direction   */
	    const real_t   r     ,  /* Input: location of r in [-1, 1]       */
	    const real_t   s     ); /* Input: location of s in [-1, 1]       */

void dglldpc (const int_t    np,    /* input: number of points for Leg polys */
	      const real_t** cd);   /* output: pointer to table of coeffs.   */

void dglldpt (const int_t    np,    /* input:  number of points for DLT      */
	      const real_t** fw,    /* output: 1D forward transform matrix   */
	      const real_t** ft,    /* output: transpose of fw               */
	      const real_t** bw,    /* output: 1D inverse transform matrix   */
	      const real_t** bt,    /* output: transpose of bw               */
	      const real_t** fu,    /* output: 2D forward transform matrix   */
	      const real_t** bu);   /* output: 2D inverse transform matrix   */

/* -- Routines from mapping.c: */

void edgemaps (const int_t np, const int_t dim,
	       int_t** emap, int_t** pmap);

/* -- Routines from family.c: */

void iadopt   (const int_t, int_t**);
void dadopt   (const int_t, real_t** );
void sadopt   (const int_t, float**  );

void iabandon (int_t**);
void dabandon (real_t**);
void sabandon (float**);

int_t FamilySize (int_t*, int_t*, int_t*);

/* -- Routines from RCM.f: */

void genrcm_ (int_t*, int_t*, int_t*, int_t*, int_t*, int_t*);
#define genrcm(neqns, xadj, adjncy, perm, mask, xls) \
(_vecIreg[0] = neqns, genrcm_(_vecIreg, xadj, adjncy, perm, mask, xls))

void F77NAME(fnroot) (int_t*, int_t*, int_t*, int_t*,
	      int_t*, int_t*, int_t*);
#define fnroot(root, xadj, adncy, mask, nlvl, xls, ls) \
(_vecIreg[0] = root, _vecIreg[1] = nlvl,                   \
 fnroot_(_vecIreg, xadj, adjncy, mask, _vecIreg + 1, xls, ls))

void rcm_    (int_t*, int_t*, int_t*, int_t*,
	      int_t*, int_t*, int_t*);
#define rcm(root, xadj, adjncy, mask, perm, ccsize, deg)  \
(_vecIreg[0] = root, rcm_(_vecIreg, xadj, adjncy, mask, perm, ccsize, deg))

/* -- Routines from fftpack.f (NETLIB/FFTPACK): */

void drffti_ (int_t*, real_t*, int_t*);
void drfftf_ (int_t*, real_t*, real_t*, real_t*, int_t*);
void drfftb_ (int_t*, real_t*, real_t*, real_t*, int_t*);

#define drffti(n,wa,ifac) \
  (_vecIreg[0]=n, drffti_ (_vecIreg, wa, ifac))
#define drfftf(n,c,ch,wa,ifac) \
  (_vecIreg[0]=n, drfftf_ (_vecIreg, c, ch, wa, ifac))
#define drfftb(n,c,ch,wa,ifac) \
  (_vecIreg[0]=n, drfftb_ (_vecIreg, c, ch, wa, ifac))

/* -- Routines from canfft.f (Canuto FFT routines): */

void factor_ (int_t*, int_t*, int_t*);

#define factor(n, nfac, ifac) (_vecIreg[0]=n, factor_ (_vecIreg, nfac, ifac))

#define dpreft(n,nfac,ifac,trig)  \
(_vecIreg[0]=n, dpreft_ (_vecIreg, nfac, ifac, trig))
#define dmrcft(v, np, nz, w, nfac, ifac, trig, sign)  \
(_vecIreg[0]=np, _vecIreg[1]=nz, _vecIreg[2]=nfac, _vecIreg[3]=sign,  \
dmrcft_(v, _vecIreg, _vecIreg+1, w, _vecIreg+2, ifac, trig, _vecIreg+3))

/* -- Routines from temfftd.f (Temperton FFT routines): */

void prf235_ (int_t*, int_t*, int_t*, int_t*, int_t*);
#define prf235(n,ip,iq,ir,ipqr2) (prf235_(n,ip,iq,ir,ipqr2))

void dsetpf_ (real_t*, int_t*, int_t*, int_t*, int_t*);
void dmpfft_ (real_t*, real_t*,  int_t*, int_t*, int_t*,
	      int_t*, int_t*, real_t*, int_t*);

#define dsetpf(trig,n,ip,iq,ir)                               \
(_vecIreg[0]=n,_vecIreg[1]=ip,_vecIreg[2]=iq, _vecIreg[3]=ir, \
dsetpf_(trig,_vecIreg,_vecIreg+1,_vecIreg+2,_vecIreg+3))
#define dmpfft(v,w,np,nz,ip,iq,ir,trig,isign)                 \
(_vecIreg[0]=np,_vecIreg[1]=nz,_vecIreg[2]=ip,                \
_vecIreg[3]=iq, _vecIreg[4]=ir,_vecIreg[5]=isign,             \
dmpfft_(v,w,_vecIreg,_vecIreg+1,_vecIreg+2,                   \
_vecIreg+3,_vecIreg+4,trig,_vecIreg+5))

/* -- Routines from matops.F: */

void dgrad2_ (real_t*,real_t*,real_t*,real_t*,real_t*,real_t*,int_t*,int_t*);

#define dgrad2(u,v,ur,vs,dv,dt,np,nel)                          \
(_vecIreg[0]=np,_vecIreg[1]=nel,dgrad2_(u,v,ur,vs,dv,dt,_vecIreg,_vecIreg+1))

#if defined(SX)

/* -- Routines from NEC FFT library: floating precision depends on library. */

void rftfax_ (int_t*,int_t*,real_t*);
void rfft_   (real_t*,real_t*,real_t*,int_t*,int_t*,int_t*,real_t*);

#define rftfax(n,ifax,trigs)                        \
(_vecIreg[0]=n,rftfax_(_vecIreg,ifax,trigs))
#define rfft(r,w,trigs,ifax,n,l,xnorm)              \
(_vecIreg[0]=n,_vecIreg[1]=l,_vecDreg[0]=xnorm,     \
rfft_(r,w,trigs,ifax,_vecIreg,_vecIreg+1,_vecDreg))

#endif

/* -- Routines from fourier.c */

void dDFTr (real_t*, const int_t, const int_t, const int_t);

/* -- Routines from filter.c */

void bvdFilter (const int_t,const int_t,const int_t, const real_t, real_t*);

/* -- Routines from message.c: */

void message_init      (int*, char***);
void message_stop      ();
void message_sync      ();
void message_dsend     (real_t*, const int_t, const int_t);
void message_drecv     (real_t*, const int_t, const int_t);
void message_ssend     (float*,  const int_t, const int_t);
void message_srecv     (float*,  const int_t, const int_t);
void message_isend     (int_t*,  const int_t, const int_t);
void message_irecv     (int_t*,  const int_t, const int_t);
void message_dexchange (real_t*, const int_t, const int_t,const int_t);
void message_sexchange (float*,  const int_t, const int_t,const int_t);
void message_iexchange (int_t*,  const int_t, const int_t,const int_t);

#endif
