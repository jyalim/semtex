#ifndef FEMLIB_H
#define FEMLIB_H
///////////////////////////////////////////////////////////////////////////////
// C++ function prototypes for routines in library libfem.a
//
// $Id: femlib.h,v 8.1 2015/04/20 11:14:14 hmb Exp $
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cmath>

#include "cfemdef.h"
#include "polylib.h"

extern "C" {

extern char buf[STR_MAX];

// -- Routines from initial.y:

void   yy_initialize (void);
real_t yy_interpret  (const char*);

void   yy_vec_init   (const char*, const char*);
void   yy_vec_interp (const int_t, ...);

void   yy_help       (void);
void   yy_show       (void);

// -- Routines from polyops.c:

void   dermat_g (const int_t, const real_t*, const int_t,
		 const real_t*, real_t*, real_t*);
void   intmat_g (const int_t, const real_t*, const int_t,
		 const real_t*, real_t*, real_t*);
void   uniknot  (const int_t, double*);

real_t PNLEG    (const real_t, const int_t);
real_t PNMOD    (const real_t, const int_t);	
void   ZWGLL    (real_t*, real_t*, const int_t);

// -- Routines from operators.c:

void zquad (const real_t** point ,  // Quadrature points.
	    const real_t** weight,  // Quadrature weights.
	    const real_t** dv    ,  // Derivative operator at points.
	    const real_t** dt    ,  // (Transpose) derivative operator.
	    const int_t    np    ,  // Input: Number of quadrature points.
	    const char     rule  ,  // Input: 'G'auss, 'R'adau, or 'L'obatto.
	    const real_t   alpha ,  // Input: Jacobi polynomial constant.
	    const real_t   beta  ); // Input: Jacobi polynomial constant.

void proj  (const real_t** IN    ,  // Interpolant operator matrix
	    const real_t** IT    ,  // Transposed interpolant operator matrix
	    const int_t    nfrom ,  // Input: Number of "from" points
	    const char     rulefr,  // Input: 'G', 'R', 'L', or 'U', "from"
	    const real_t   alphfr,  // Input: Jacobi constant, "from"
	    const real_t   betafr,  // Input: Jacobi constant, "from"
	    const int_t    nto   ,  // Input: Number of "to" points
	    const char     ruleto,  // Input: 'G', 'R', 'L', or 'U', "to"
	    const real_t   alphto,  // Input: Jacobi constant, "to"
	    const real_t   betato); // Input: Jacobi constant, "to"

void intp  (real_t*        inr   ,  // 1D shape function at r
	    real_t*        ins   ,  // 1D shape function at s
	    real_t*        dvr   ,  // 1D shape function derivative at r
	    real_t*        dvs   ,  // 1D shape function derivative at s
	    const int_t    nr    ,  // Input: number of points, r-direction
	    const char     basisr,  // Input: 'G', 'R', 'L', r-direction
	    const real_t   alphar,  // Input: Jacobi constant, r-direction
	    const real_t   betar ,  // Input: Jacobi constant, r-direction
	    const int_t    ns    ,  // Input: number of points, s-direction
	    const char     basiss,  // Input: 'G', 'R', 'L', s-direction
	    const real_t   alphas,  // Input: Jacobi constant, s-direction
	    const real_t   betas ,  // Input: Jacobi constant, s-direction
	    const real_t   r     ,  // Input: location of r in [-1, 1]
	    const real_t   s     ); // Input: location of s in [-1, 1]

void dglldpc (const int_t    np,     // input:  number of points for Leg polys
	      const real_t** cd);    // output: pointer to table of coeffs

void dglldpt (const int_t    np,     // input:  number of points for DLT
	      const real_t** fw,     // output: 1D forward transform matrix
	      const real_t** ft,     // output: transpose of fw
	      const real_t** bw,     // output: 1D inverse transform matrix
	      const real_t** bt,     // output: transpose of bw
	      const real_t** fu,     // output: 2D forward transform matrix
	      const real_t** bu);    // output: 2D inverse transform matrix

void dglmdpc (const int_t    np,     // input:  no. of points for expansions.
	      const real_t** cd);    // output: pointer to table of coeffs

void dglmdpt (const int_t    np,     // input:  number of points for DPT
	      const real_t** fw,     // output: 1D forward transform matrix
	      const real_t** ft,     // output: transpose of fw
	      const real_t** bw,     // output: 1D inverse transform matrix
	      const real_t** bt,     // output: transpose of bw
	      const real_t** fu,     // output: 2D forward transform matrix
	      const real_t** bu);    // output: 2D inverse transform matrix

// -- Routines from mapping.c:

void edgemaps (const int_t np, const int_t dim, int_t** emap, int_t** pmap);

// -- Routines from family.c:

void iadopt   (const int_t, int_t**);
void dadopt   (const int_t, real_t**);
void sadopt   (const int_t, float**);

void iabandon (int_t**);
void dabandon (real_t**);
void sabandon (float**);

int_t  FamilySize (int_t*, int_t*, int_t*);

// -- Routines from fourier.c:

void dDFTr (real_t*, const int_t, const int_t, const int_t);

// -- Routines from filter.c

void bvdFilter (const int_t,const real_t,const real_t, const real_t, real_t*);

// -- Routines from message.c:

void message_init      (int*, char***);
void message_stop      ();
void message_sync      ();
void message_dsend     (real_t*, const int_t, const int_t);
void message_drecv     (real_t*, const int_t, const int_t);
void message_ssend     (float*,  const int_t, const int_t);
void message_srecv     (float*,  const int_t, const int_t);
void message_isend     (int_t*,  const int_t, const int_t);
void message_irecv     (int_t*,  const int_t, const int_t);
void message_dexchange (real_t*, const int_t, const int_t, const int_t);
void message_sexchange (float*,  const int_t, const int_t, const int_t);
void message_iexchange (int_t*,  const int_t, const int_t, const int_t);

// -- FORTRAN
// -- Routines from NETLIB.f:

void F77NAME(genrcm) (const int_t&,int_t*,int_t*,int_t*,int_t*,int_t*);
void F77NAME(fnroot) (int_t&,int_t*,int_t*,int_t*,int_t&,int_t*,int_t*);
void F77NAME(rcm)    (const int_t&,int_t*,int_t*,int_t*,int_t*,int_t&,int_t*);

void F77NAME(drffti) (const int_t&,real_t*,int_t*);
void F77NAME(drfftf) (const int_t&,real_t*,real_t*,const real_t*,const int_t*);
void F77NAME(drfftb) (const int_t&,real_t*,real_t*,const real_t*,const int_t*);

double F77NAME(fbrent) (const real_t&, const real_t&,
			real_t(*)(const real_t&), const real_t&);
void   F77NAME(braket) (real_t&, real_t&, real_t&, real_t&, real_t&, real_t&,
			real_t(*)(const real_t&));

// -- Routines from canfft.f:

void F77NAME(factor) (const int_t&, int_t&, int_t*);

void F77NAME(dpreft) (const int_t&, int_t&, int_t*, real_t*);
void F77NAME(dmrcft) (real_t*, const int_t&, const int_t&, real_t*,
		      const int_t&, int_t*, real_t*, const int_t&);
void F77NAME(dfft1)  (real_t*, real_t*, const int_t&, const int_t&,
		      const int_t*, const int_t&, const real_t*,
		      const int_t&);

// -- Routines from temfft.f:

void F77NAME(prf235) (int_t&, int_t&, int_t&, int_t&, int_t&);

void F77NAME(dsetpf) (const real_t*, const int_t&, const int_t&,
		      const int_t&, const int_t&);
void F77NAME(dmpfft) (real_t*, real_t*, const int_t&, const int_t&,
		      const int_t&, const int_t&, const int_t&, 
		      const real_t*, const int_t&);
void F77NAME(dgpfa)  (real_t*, real_t*, const real_t*, const int_t&,
		      const int_t&, const int_t&, const int_t&,
		      const int_t&, const int_t&, const int_t&,
		      const int_t&);

// -- Routines from matops.F:

void F77NAME(dgrad2) (const real_t*, const real_t*, real_t*, real_t*,
		      const real_t*, const real_t*, const int_t&,
		      const int_t&, const int_t&);
void F77NAME(dtpr2d) (const real_t*, real_t*, real_t*, const real_t*,
		      const real_t*, const int_t&, const int_t&, const int_t&);

}

class Femlib {
public:

  static void initialize (int* argc, char*** argv)
    { message_init (argc, argv); }
  static void prepVec (const char* v, const char* f)
    { yy_vec_init (v, f); }
  static void finalize () 
    { message_stop (); }
  static void synchronize ()
    { message_sync(); }

#if 0
  // -- This form seems to cause a problem, hence workaround:
  static void parseVec (const int_t n ... )
    { yy_vec_interp (n ... ); }
#else
  #define Femlib__parseVec yy_vec_interp
#endif

  static void value (const char* s, const real_t p)
    { sprintf (buf, "%s = %.17g", s, p); yy_interpret (buf); }
  static void ivalue (const char* s, const int_t p)
    { sprintf (buf, "%s = %1d", s, p); rint (yy_interpret (buf)); }
  static real_t value (const char* s)
    { return yy_interpret (s); }
  static int_t ivalue (const char* s)
    { return rint (yy_interpret (s)); }
  
  static void equispacedMesh (const int_t np, real_t* z)
    { uniknot (np, z); }

  static void LagrangeInt (const int_t N,  const real_t* zero,
			   const int_t I,  const real_t* x   ,
			   real_t*     IN, real_t*       IT  )
    { intmat_g (N, zero, I, x, IN, IT); }
  static void LagrangeDer (const int_t N,  const real_t* zero,
			   const int_t I,  const real_t* x   ,
			   real_t*     IN, real_t*       IT  )
    { dermat_g (N, zero, I, x, IN, IT); }

  static void   GLLzw (const int_t N,real_t* z,real_t* w)  {ZWGLL(z,w,N);}
  static real_t LegendreVal (const int_t N,const real_t x) {return PNLEG(x,N);}
  static real_t ModalVal    (const int_t N,const real_t x) {return PNMOD(x,N);}

  static void quadrature (const real_t** point ,
			  const real_t** weight,
			  const real_t** dv    ,
			  const real_t** dt    ,
			  const int_t    np    ,
			  const char     rule  ,
			  const real_t   alpha ,
			  const real_t   beta  )
  { zquad (point, weight, dv, dt, np, rule, alpha, beta); }

  static void projection (const real_t** IN    ,
			  const real_t** IT    ,
			  const int_t    nfrom ,
			  const char     rulefr,
			  const real_t   alphfr,
			  const real_t   betafr,
			  const int_t    nto   ,
			  const char     ruleto,
			  const real_t   alphto,
			  const real_t   betato)
  { proj (IN, IT, nfrom,rulefr,alphfr,betafr, nto,ruleto,alphto,betato); }

  static void interpolation (real_t*      inr   ,
			     real_t*      ins   ,
			     real_t*      dvr   ,
			     real_t*      dvs   ,
			     const int_t  nr    ,
			     const char   basisr,
			     const real_t alphar,
			     const real_t betar ,
			     const int_t  ns    ,
			     const char   basiss,
			     const real_t alphas,
			     const real_t betas ,
			     const real_t r     ,
			     const real_t s     )
  { intp (inr,ins,dvr,dvs,nr,basisr,alphar,betar,ns,basiss,alphas,betas,r,s); }

  static void legCoef (const int_t n,
		       const real_t**  c)
    { dglldpc (n, c); }

  static void legTran (const int_t n,
		       const real_t** f1, const real_t** ft,
		       const real_t** i1, const real_t** it,
		       const real_t** f2, const real_t** i2)
    { dglldpt (n, f1, ft, i1, it, f2, i2); }

  static void modCoef (const int_t    n,
		       const real_t** c)
    { dglmdpc (n, c); }

  static void modTran (const int_t n,
		       const real_t** f1, const real_t** ft,
		       const real_t** i1, const real_t** it,
		       const real_t** f2, const real_t** i2)
    { dglmdpt (n, f1, ft, i1, it, f2, i2); }

  static void buildMaps (const int_t np, const int_t dim,
			 int_t** e, int_t** i)
    { edgemaps (np, dim, e, i); }

  static void adopt (const int_t np, int_t** v)
    { iadopt (np, v); }
  static void adopt (const int_t np, real_t** v)
    { dadopt (np, v); }
  static void adopt (const int_t np, float** v)
    { sadopt (np, v); }
  static void abandon (int_t** v)
    { iabandon (v); }
  static void abandon (real_t** v)
    { dabandon (v); }
  static void abandon (float** v)
    { sabandon (v); }
  static int_t  fwords (int_t* ni, int_t* nd, int_t* ns)
    { return FamilySize (ni, nd, ns); }

  static void genrcm (const int_t& n, int_t* x, int_t* a,
		      int_t* p, int_t* m, int_t* l) 
    { F77NAME(genrcm) (n, x, a, p, m, l); }
  static void fnroot (int_t& r, int_t* x, int_t* a, int_t* m,
		      int_t& n, int_t* l, int_t* p) 
    { F77NAME(fnroot) (r, x, a, m, n, l, p); }
  static void rcm (const int_t& n, int_t* x, int_t* a, int_t* m,
		   int_t* p, int_t& c, int_t* l)
    { F77NAME(rcm) (n, x, a, m, p, c, l); }

  static void rffti (const int_t& n , real_t* w, int_t* i)
    { F77NAME(drffti) (n, w, i); }
  static void rfftf (const int_t& n, real_t* c, real_t* ch,
	const real_t* wa, const int_t* ifac)
    { F77NAME(drfftf) (n, c, ch, wa, ifac); }
  static void rfftb (const int_t& n, real_t* c, real_t* ch,
  	const real_t* wa, const int_t* ifac) 
    { F77NAME(drfftb) (n, c, ch, wa, ifac); }

  static real_t brent  (const real_t& ax, const real_t& bx,
			real_t(*f)(const real_t&),
			const real_t& tol) {
    return F77NAME(fbrent) (ax, bx, f, tol);
  }
  static void bracket   (real_t& ax, real_t& bx, real_t& cx,
			 real_t& fa, real_t& fb, real_t& fc,
			 real_t(*func)(const real_t&)) {
    F77NAME(braket) (ax, bx, cx, fa, fb, fc, func);
  }
  static void DFTr (real_t*  data, const int_t nz, const int_t np,
		    const int_t sign)
    { dDFTr (data, nz, np, sign); }

  static void preft  (const int_t& n, int_t& nfax, int_t* ifax,
		      real_t* trig) 
    { F77NAME(dpreft) (n, nfax, ifax, trig); }
  static void mrcft  (real_t* v, const int_t& np, const int_t& nz,
		      real_t* w, const int_t& nfx, int_t* ifx,
		      real_t* trig, const int_t& sign)
    { F77NAME(dmrcft) (v, np, nz, w, nfx, ifx, trig, sign); }
  static void primes23 (const int_t& n, int_t& np, int_t* fact)
    { F77NAME(factor) (n, np, fact); }
  static void fft1 (real_t* a, real_t* c, const int_t& n,
		    const int_t& nfax,
                    const int_t* ifax, const int_t& isign,
		    const real_t* trig, const int_t& len)
    { F77NAME(dfft1) (a,c,n,nfax,ifax,isign,trig,len); }  

  static void primes235 (int_t& n,
			 int_t&ip, int_t&iq, int_t& ir, int_t& ipqr2)
    { F77NAME(prf235) (n, ip, iq, ir, ipqr2); }
  static void setpf (const real_t* t, const int_t& n, const int_t& ip,
		     const int_t& iq, const int_t& ir)
    { F77NAME(dsetpf) (t, n, ip, iq, ir); }
  static void mpfft (real_t* v, real_t* w, const int_t& np,
		     const int_t& nz,
		     const int_t& ip, const int_t& iq, const int_t& ir,
		     const real_t* trig, const int_t& isign)
    { F77NAME(dmpfft) (v, w, np, nz, ip, iq, ir, trig, isign); }
  static void gpfa (real_t* a, real_t* b, const real_t* trig,
		    const int_t& inc, const int_t& jump, const int_t& n,
		    const int_t& ip, const int_t& iq, const int_t& ir,
		    const int_t& lot, const int_t& isign)
    { F77NAME(dgpfa) (a, b, trig, inc, jump, n, ip, iq, ir, lot, isign); }

  static void send  (real_t* data, const int_t N, const int_t tgt)
    { message_dsend (data, N, tgt); }
  static void recv  (real_t* data, const int_t N, const int_t src)
    { message_drecv (data, N, src); }
  static void send  (float*  data, const int_t N, const int_t tgt)
    { message_ssend (data, N, tgt); }
  static void recv  (float*  data, const int_t N, const int_t src)
    { message_srecv (data, N, src); }
  static void send  (int_t*  data, const int_t N, const int_t tgt)
    { message_isend (data, N, tgt); }
  static void recv  (int_t*  data, const int_t N, const int_t src)
    { message_irecv (data, N, src); }
  static void exchange  (real_t* data, const int_t nZ,
			 const int_t nP, const int_t sign)
    { message_dexchange (data, nZ, nP, sign); }
  static void exchange  (float* data, const int_t nZ,
			 const int_t nP, const int_t sign)
    { message_sexchange (data, nZ, nP, sign); }
  static void exchange  (int_t* data,  const int_t nZ,
			 const int_t nP, const int_t sign)
    { message_iexchange (data, nZ, nP, sign); }

  static void grad2 (const real_t* x, const real_t* y, real_t* xr, real_t* ys,
		     const real_t* dv, const real_t* dt,
		     const int_t& nr, const int_t& ns, const int_t& nel)
    { F77NAME(dgrad2) (x, y, xr, ys, dv, dt, nr, ns, nel); }

  static void tpr2d (const real_t* x, real_t* y, real_t* t,
		     const real_t* dv, const real_t* dt,
		     const int_t& nr, const int_t& ns, const int_t& nel)
    { F77NAME(dtpr2d) (x, y, t, dv, dt, nr, ns, nel); }

  static void erfcFilter (const int_t N, const real_t p, const real_t s,
			  const real_t a, real_t* f)
    { bvdFilter (N, p, s, a, f); }

};

#endif
