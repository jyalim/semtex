#ifndef VECLIB_H
#define VECLIB_H
///////////////////////////////////////////////////////////////////////////////
// Veclib.h:  C++ wrappers for veclib subroutine calls.
//
// Veclib is described in the iPSC/2 Programmer's Reference Manual.
//
// $Id: veclib.h,v 8.1 2015/04/20 11:14:19 hmb Exp $
///////////////////////////////////////////////////////////////////////////////

#include <cfemdef.h>

extern "C" {

  /* -- MATHEMATICAL PRIMITIVES */

  void  dcopy  (int_t n,
		const double* x, int_t incx, double* y, int_t incy);
  void  icopy  (int_t n,
		const int_t*  x, int_t incx, int_t*  y, int_t incy);
  void  scopy  (int_t n,
		const float*  x, int_t incx, float*  y, int_t incy);
  
  void  dfill  (int_t n, double alpha, double* x, int_t incx);
  void  ifill  (int_t n, int_t  alpha, int_t*  x, int_t incx);
  void  sfill  (int_t n, float  alpha, float*  x, int_t incx);
  
  void  dznan  (int_t n, double* x, int_t incx);
  void  sznan  (int_t n, float*  x, int_t incx);
  
  void  dneg   (int_t n, double* x, int_t incx);
  void  ineg   (int_t n, int_t*  x, int_t incx);
  void  sneg   (int_t n, float*  x, int_t incx);
  
  void  dvneg  (int_t n,
		const double* x, int_t incx, double* y, int_t incy);
  void  ivneg  (int_t n,
		const int_t*  x, int_t incx, int_t*  y, int_t incy);
  void  svneg  (int_t n,
		const float*  x, int_t incx, float*  y, int_t incy);
  
  void  dvsgn  (int_t n,
		const double* x, int_t incx, double* y, int_t incy);
  void  ivsgn  (int_t n,
		const int_t*  x, int_t incx, int_t*  y, int_t incy);
  void  svsgn  (int_t n,
		const float*  x, int_t incx, float*  y, int_t incy);
  
  void  dsadd  (int_t n, double  alpha,
		const double*  x, int_t incx, double*  y, int_t incy);
  void  isadd  (int_t n, int_t alpha,
		const int_t* x, int_t incx, int_t* y, int_t incy);
  void  ssadd  (int_t n, float  alpha,
		const float*   x, int_t incx, float*   y, int_t incy);
  
  void  dspow  (const int_t n, const double alpha,
		const double* x, int_t incx, double* y, int_t incy);
  void  sspow  (const int_t n, const float  alpha,
		const float*  x, int_t incx, float*  y, int_t incy);
  
  void  dvadd  (int_t n, const double* x, int_t incx,
		const double* y, int_t incy, double* z, int_t incz);
  void  ivadd  (int_t n, const int_t*  x, int_t incx, 
		const int_t*  y, int_t incy, int_t*  z, int_t incz);
  void  svadd  (int_t n, const float*  x, int_t incx,
		const float*  y, int_t incy, float*  z, int_t incz);
  
  void  dssub  (int_t n, double alpha,
		const double* x, int_t incx, double* y, int_t incy);
  void  issub  (int_t n, int_t alpha,
		const int_t*  x, int_t incx, int_t*  y, int_t incy);
  void  sssub  (int_t n, float   alpha,
		const float*  x, int_t incx, float*  y, int_t incy);
  
  void  dvsub  (int_t n, const double* x, int_t incx,
		const double* y, int_t incy, double* z, int_t incz);
  void  ivsub  (int_t n, const int_t*  x, int_t incx,
		const int_t*  y, int_t incy, int_t*  z, int_t incz);
  void  svsub  (int_t n, const float*  x, int_t incx,
		const float*  y, int_t incy, float*  z, int_t incz);
  
  void  dsmul  (int_t n, double alpha,
		const double* x, int_t incx, double* y, int_t incy);
  void  ismul  (int_t n, int_t alpha,
		const int_t*  x, int_t incx, int_t*  y, int_t incy);
  void  ssmul  (int_t n, float  alpha,
		const float*  x, int_t incx, float*  y, int_t incy);
  
  void  dvmul  (int_t n, const double* x, int_t incx,
		const double* y, int_t incy, double* z, int_t incz);
  void  ivmul  (int_t n, const int_t*  x, int_t incx,
		const int_t*  y, int_t incy, int_t*  z, int_t incz);
  void  svmul  (int_t n, const float*  x, int_t incx,
		const float*  y, int_t incy, float*  z, int_t incz);
  
  void  dsdiv  (int_t n, double alpha,
		const double* x, int_t incx, double* y, int_t incy);
  void  isdiv  (int_t n, int_t  alpha,
		const int_t*  x, int_t incx, int_t*  y, int_t incy);
  void  ssdiv  (int_t n, float  alpha,
		const float*  x, int_t incx, float*  y, int_t incy);
  
  void  dvrecp (int_t n,
		const double* x, int_t incx, double* y, int_t incy);
  void  svrecp (int_t n,
		const float*  x, int_t incx, float*  y, int_t incy);
  
  void  dvdiv  (int_t n, const double* x, int_t incx,
		const double* y, int_t incy, double* z, int_t incz);
  void  svdiv  (int_t n, const float*  x, int_t incx,
		const float*  y, int_t incy, float*  z, int_t incz);
  
  void  dzero  (int_t n, double* x, int_t incx);
  void  izero  (int_t n, int_t*  x, int_t incx);
  void  szero  (int_t n, float*  x, int_t incx);
  
  /* -- OTHER MATHEMATICAL FUNCTIONS */
  
  void    dvabs    (int_t n,
		    const double*  x, int_t incx, double*  y, int_t incy);
  void    ivabs    (int_t n,
		    const int_t* x, int_t incx, int_t* y, int_t incy);
  void    svabs    (int_t n,
		    const float*   x, int_t incx, float*   y, int_t incy);
  
  void    dvamax   (int_t n, const double*  x, int_t incx,
		    const double*  y, int_t incy, double*  z, int_t incz);
  void    ivamax   (int_t n, const int_t* x, int_t incx,
		    const int_t* y, int_t incy, int_t* z, int_t incz);
  void    svamax   (int_t n, const float*   x, int_t incx,
		    const float*   y, int_t incy, float*   z, int_t incz);
  
  void    dvexp    (int_t n,
		    const double* x, int_t incx, double* y, int_t incy);
  void    svexp    (int_t n,
		    const float*  x, int_t incx, float*  y, int_t incy);
  
  void    dvlg10   (int_t n,
		    const double* x, int_t incx, double* y, int_t incy);
  void    svlg10   (int_t n,
		    const float*  x, int_t incx, float*  y, int_t incy);
  
  void    dvlog    (int_t n,
		    const double* x, int_t incx, double* y, int_t incy);
  void    svlog    (int_t n,
		    const float*  x, int_t incx, float*  y, int_t incy);
  
  void    dvatan   (int_t n,
		    const double* x, int_t incx, double* y, int_t incy);
  void    svatan   (int_t n,
		    const float*  x, int_t incx, float*  y, int_t incy);
  
  void    dvatn2   (int_t n, const double* x, int_t incx,
		    const double* y, int_t incy, double* z, int_t incz);
  void    svatn2   (int_t n, const float*  x, int_t incx,
		    const float*  y, int_t incy, float*  z, int_t incz);
  
  void    dvcos    (int_t n,
		    const double* x, int_t incx, double* y, int_t incy);
  void    svcos    (int_t n,
		    const float*  x, int_t incx, float*  y, int_t incy);
  
  void    dvsin    (int_t n,
		    const double* x, int_t incx, double* y, int_t incy);
  void    svsin    (int_t n,
		    const float*  x, int_t incx, float*  y, int_t incy);
  
  void    dvsqrt   (int_t n,
		    const double* x, int_t incx, double* y, int_t incy);
  void    svsqrt   (int_t n,
		    const float*  x, int_t incx, float*  y, int_t incy);
  
  void    dvtanh   (int_t n,
		    const double* x, int_t incx, double* y, int_t incy);
  void    svtanh   (int_t n,
		    const float*  x, int_t incx, float*  y, int_t incy);
  
  void    raninit  (int_t flag);
  double  dranu    (void);
  float   sranu    (void);
  double  drang    (void);
  float   srang    (void);
  double  dnormal  (double mean, double sdev);
  float   snormal  (float  mean, float  sdev);
  void    dvrandom (int_t  n, double* x, int_t incx);
  void    svrandom (int_t  n, float*  x, int_t incx);
  void    dvgauss  (int_t  n, double* x, int_t incx);
  void    svgauss  (int_t  n, float*  x, int_t incx);
  void    dvnormal (int_t  n, double mean, double sdev,
		    double* x, int_t incx);
  void    svnormal (int_t  n, float  mean, float  sdev, 
		    float*  x, int_t incx);

  void    dvhypot  (int_t n, const double* x, int_t incx,
		    const double* y, int_t incy, double* z, int_t incz);
  void    svhypot  (int_t n, const float*  x, int_t incx,
		    const float*  y, int_t incy, float*  z, int_t incz);
  void    dvmag    (int_t n, const double* w, int_t incw, 
		    const double* x, int_t incx, const double* y,
		    int_t incy, double* z, int_t incz);
  void    svmag    (int_t n, const float* w, int_t incw, 
		    const float*  x, int_t incx, const float*  y,
		    int_t incy, float*  z, int_t incz);

  void    dvpow    (int_t n, const double* x, int_t incx,
		    const double* y, int_t incy, double* z, int_t incz);
  void    svpow    (int_t n, const float*  x, int_t incx,
		    const float*  y, int_t incy, float*  z, int_t incz);

  
  /* -- TRIAD OPERATIONS */

  void   dsvmvt (int_t n, double alpha, const double* x, int_t incx,
		 const double* y, int_t incy, double* z, int_t incz);
  void   ssvmvt (int_t n, float  alpha, const float*  x, int_t incx,
		 const float*  y, int_t incy, float*  z, int_t incz);
  
  void   dsvpvt (int_t n, double alpha, const double* x, int_t incx,
		 const double* y, int_t incy, double* z, int_t incz);
  void   ssvpvt (int_t n, float  alpha, const float*  x, int_t incx,
		 const float*  y, int_t incy, float*  z, int_t incz);
  
  void   dsvtsp (int_t n, double alpha, double beta,
		 const double* x, int_t incx, double* y, int_t incy);
  void   ssvtsp (int_t n, float  alpha, float  beta,
		 const float*  x, int_t incx, float*  y, int_t incy);
  
  void   dsvtvm (int_t n, double alpha, const double* x, int_t incx,
		 const double* y, int_t incy, double* z, int_t incz);
  void   ssvtvm (int_t n, float  alpha, const float*  x, int_t incx,
		 const float*  y, int_t incy, float*  z, int_t incz);
  
  void   dsvtvp (int_t n, double alpha, const double* x, int_t incx,
		 const double* y, int_t incy, double* z, int_t incz);
  void   ssvtvp (int_t n, float  alpha, const float*  x, int_t incx,
		 const float*  y, int_t incy, float*  z, int_t incz);
  
  void   dsvvmt (int_t n, double alpha, const double* x, int_t incx,
		 const double* y, int_t incy, double* z, int_t incz);
  void   ssvvmt (int_t n, float  alpha, const float*  x, int_t incx,
		 const float*  y, int_t incy, float*  z, int_t incz);
  
  void   dsvvpt (int_t n, double alpha, const double* x, int_t incx,
		 const double* y, int_t incy, double* z, int_t incz);
  void   ssvvpt (int_t n, float  alpha, const float*  x, int_t incx,
		 const float*  y, int_t incy, float*  z, int_t incz);
  
  void   dsvvtm (int_t n, double alpha, const double* x, int_t incx,
		 const double* y, int_t incy, double* z, int_t incz);
  void   ssvvtm (int_t n, float  alpha, const float*  x, int_t incx,
		 const float*  y, int_t incy, float*  z, int_t incz);
  
  void   dsvvtp (int_t n, double alpha, const double* x, int_t incx,
		 const double* y, int_t incy, double* z, int_t incz);
  void   ssvvtp (int_t n, float  alpha, const float*  x, int_t incx,
		 const float*  y, int_t incy, float*  z, int_t incz);

  void   dsvvtt (int_t n, double alpha, const double* x, int_t incx,
		 const double* y, int_t incy, double* z, int_t incz);
  void   ssvvtt (int_t n, float  alpha, const float*  x, int_t incx,
		 const float*  y, int_t incy, float*  z, int_t incz);
  
  void   dvvmvt (int_t n, const double* w, int_t incw, const double* x,
		 int_t incx, const double* y, int_t incy,
		 double* z, int_t incz);
  void   svvmvt (int_t n, const float*  w, int_t incw, const float*  x,
		 int_t incx, const float*  y, int_t incy,
		 float*  z, int_t incz);
  
  void   dvvpvt (int_t n, const double* w, int_t incw, const double* x,
		 int_t incx, const double* y, int_t incy,
		 double* z, int_t incz);
  void   svvpvt (int_t n, const float*  w, int_t incw, const float*  x,
		 int_t incx, const float*  y, int_t incy,
		 float*  z, int_t incz);

  void   dvvtvm (int_t n, const double* w, int_t incw, const double* x,
		 int_t incx, const double* y, int_t incy,
		 double* z, int_t incz);
  void   svvtvm (int_t n, const float*  w, int_t incw, const float*  x,
		 int_t incx, const float*  y, int_t incy,
		 float*  z, int_t incz);
  
  void   dvvtvp (int_t n, const double* w, int_t incw, const double* x,
		 int_t incx, const double* y, int_t incy,
		 double* z, int_t incz);
  void   svvtvp (int_t n, const float*  w, int_t incw, const float*  x,
		 int_t incx, const float*  y, int_t incy,
		 float*  z, int_t incz);
  
  void   dvvvtm (int_t n, const double* w, int_t incw, const double* x,
		 int_t incx, const double* y, int_t incy,
		 double* z, int_t incz);
  void   svvvtm (int_t n, const float*  w, int_t incw, const float*  x,
		 int_t incx, const float*  y, int_t incy,
		 float*  z, int_t incz);

  void   dvvvtt (int_t n,    const double* w, int_t incw, const double* x,
		 int_t incx, const double* y, int_t incy,
		 double* z,  int_t incz);
  void   svvvtt (int_t n,    const float*  w, int_t incw, const float*  x,
		 int_t incx, const float*  y, int_t incy,
		 float*  z,  int_t incz);
  
  /* -- RELATIONAL PRIMITIVE OPERATIONS */

  void   iseq (int_t n, int_t  alpha,
	       const  int_t*  x, int_t incx, int_t *y, int_t incy);
  void   dsge (int_t n, double alpha,
	       const  double* x, int_t incx, int_t *y, int_t incy);
  void   isge (int_t n, int_t  alpha,
	       const  int_t*  x, int_t incx, int_t *y, int_t incy);
  void   ssge (int_t n, float  alpha,
	       const  float*  x, int_t incx, int_t *y, int_t incy);
  void   dsle (int_t n, double alpha,
	       const  double* x, int_t incx, int_t *y, int_t incy);
  void   isle (int_t n, int_t  alpha,
	       const  int_t*  x, int_t incx, int_t *y, int_t incy);
  void   ssle (int_t n, float  alpha,
	       const  float*  x, int_t incx, int_t *y, int_t incy);
  void   dslt (int_t n, double alpha,
	       const  double* x, int_t incx, int_t *y, int_t incy);
  void   islt (int_t n, int_t  alpha,
	       const  int_t*  x, int_t incx, int_t *y, int_t incy);
  void   sslt (int_t n, float  alpha,
	       const  float*  x, int_t incx, int_t *y, int_t incy);
  void   dsne (int_t n, double alpha,
	       const  double* x, int_t incx, int_t *y, int_t incy);
  void   isne (int_t n, int_t  alpha,
	       const  int_t*  x, int_t incx, int_t *y, int_t incy);
  void   ssne (int_t n, float  alpha,
	       const  float*  x, int_t incx, int_t *y, int_t incy);
  
  /* -- REDUCTION FUNCTIONS */

  double dsum   (int_t n, const double* x, int_t incx);
  int_t  isum   (int_t n, const int_t*  x, int_t incx);
  float  ssum   (int_t n, const float*  x, int_t incx);
  int_t  idmax  (int_t n, const double* x, int_t incx);
  int_t  iimax  (int_t n, const int_t*  x, int_t incx);
  int_t  ismax  (int_t n, const float*  x, int_t incx);
  int_t  idmin  (int_t n, const double* x, int_t incx);
  int_t  iimin  (int_t n, const int_t*  x, int_t incx);
  int_t  ismin  (int_t n, const float*  x, int_t incx);
  int_t  icount (int_t n, const int_t*  x, int_t incx);
  int_t  ifirst (int_t n, const int_t*  x, int_t incx);
  int_t  lany   (int_t n, const int_t*  x, int_t incx);  
  int_t  lisame (int_t n, const int_t*  x, int_t incx,
		          const int_t*  y, int_t incy);
  int_t  ldsame (int_t n, const double* x, int_t incx,
		          const double* y, int_t incy);
  int_t  lssame (int_t n, const float*  x, int_t incx,
		          const float*  y, int_t incy);

  /* -- CONVERSION PRIMITIVES */
  
  void   vdble   (int_t n,
		  const float*  x, int_t incx, double* y, int_t incy);
  void   vsngl   (int_t n,
		  const double* x, int_t incx, float*  y, int_t incy);
  void   dvfloa  (int_t n,
		  const int_t*  x, int_t incx, double* y, int_t incy);
  void   svfloa  (int_t n,
		  const int_t*  x, int_t incx, float*  y, int_t incy);


  int_t iformat (void);
  void    format  (char*);
  void    dbrev   (int_t n,
		   const double* x, int_t incx, double* y, int_t incy);
  void    ibrev   (int_t n,
		   const int_t*  x, int_t incx, int_t*  y, int_t incy);
  void    sbrev   (int_t n,
		   const float*  x, int_t incx, float*  y, int_t incy);
  
  /* -- MISCELLANEOUS FUNCTIONS */

  double dclock  (void);
  float  sclock  (void);
  
  void   dscatr  (int_t n, const double* x, const int_t* y, double* z);
  void   iscatr  (int_t n, const int_t*  x, const int_t* y, int_t*  z);
  void   sscatr  (int_t n, const float*  x, const int_t* y, float*  z);
  
  void   dgathr  (int_t n, const double* x, const int_t* y, double* z);
  void   igathr  (int_t n, const int_t*  x, const int_t* y, int_t*  z);
  void   sgathr  (int_t n, const float*  x, const int_t* y, float*  z);

  void   dgathr_scatr (int_t n, const double*  w, 
		       const int_t* x, const int_t* y, double* z);
  void   igathr_scatr (int_t n, const int_t* w,
		       const int_t* x, const int_t* y, int_t*  z);
  void   sgathr_scatr (int_t n, const float*   w,
		       const int_t* x, const int_t* y,  float* z);
  
  void   dscatr_sum (int_t n,
		     const double* x, const int_t *y, double* z);
  void   iscatr_sum (int_t n,
		     const int_t*  x, const int_t* y, int_t*  z);
  void   sscatr_sum (int_t n,
		     const float*  x, const int_t *y, float*  z);
  
  void   dgathr_sum (int_t n,
		     const double* x, const int_t* y, double* z);
  void   igathr_sum (int_t n,
		     const int_t*  x, const int_t* y, int_t*  z);
  void   sgathr_sum (int_t n,
		     const float*  x, const int_t* y, float*  z);

  void   dgathr_scatr_sum (int_t n, const double*  w, 
			   const int_t* x, const int_t* y, double*  z);
  void   igathr_scatr_sum (int_t n, const int_t* w,
			   const int_t* x, const int_t* y, int_t* z);
  void   sgathr_scatr_sum (int_t n, const float*   w,
			   const int_t* x, const int_t* y, float*   z);
  
  void   dramp   (int_t n, double  alpha, double  beta,
		  double*  x, int_t incx);
  void   iramp   (int_t n, int_t alpha, int_t beta,
		  int_t* x, int_t incx);
  void   sramp   (int_t n, float   alpha, float   beta,
		  float*   x, int_t incx);
  
  void   dcndst  (int_t n, const  double*  x, int_t incx,
		  const int_t* y, int_t incy, double*  z, int_t incz);
  void   icndst  (int_t n, const  int_t* x, int_t incx,
		  const int_t* y, int_t incy, int_t* z, int_t incz);
  void   scndst  (int_t n, const  float*   x, int_t incx,
		  const int_t* y, int_t incy, float*   z, int_t incz);
  
  void   dmask   (int_t n,
		  const double* w, int_t incw,
		  const double* x, int_t incx,
		  const int_t*  y, int_t incy,
		        double* z, int_t incz);
  void   imask   (int_t n,
		  const int_t*  w, int_t incw,
		  const int_t*  x, int_t incx,
		  const int_t*  y, int_t incy,
		        int_t*  z, int_t incz);
  void   smask   (int_t n,
		  const float*  w, int_t incw,
		  const float*  x, int_t incx,
		  const int_t*  y, int_t incy,
		        float*  z, int_t incz);

  void dclip (int_t n, const  double alpha,  const double   beta,
	      const double* x,  int_t incx,  const double*  y, int_t incy);
  void iclip (int_t n, const  int_t alpha, const int_t  beta,
	      const int_t* x, int_t incx,  const int_t* y, int_t incy);
  void sclip (int_t n, const  float alpha,   const float    beta,
	      const float* x,  int_t incx,   const float*   y, int_t incy);

  void dclipup (int_t n, const double alpha, const double* x,
		int_t incx,  const double* y, int_t incy);
  void iclipup (int_t n, const int_t alpha, const int_t* x,
		int_t incx, const int_t* y, int_t incy);
  void sclipup (int_t n, const float alpha, const float* x,
		int_t incx, const float* y, int_t incy);

  void dclipdn (int_t n, const double alpha, const double* x,
		int_t incx,  const double* y, int_t incy);
  void iclipdn (int_t n, const int_t alpha, const int_t* x,
		int_t incx, const int_t* y, int_t incy);
  void sclipdn (int_t n, const float alpha, const float* x,
		int_t incx, const float* y, int_t incy);

  void diclip (int_t n, const  double alpha,  const double   beta,
	       const double* x,  int_t incx, const double*  y, int_t incy);
  void iiclip (int_t n, const  int_t alpha, const int_t  beta,
	       const int_t* x, int_t incx, const int_t* y, int_t incy);
  void siclip (int_t n, const  float alpha,   const float    beta,
	       const float* x,  int_t incx,  const float*   y, int_t incy);
  
  void   dvpoly  (int_t n, const double* x, int_t incx, int_t m,
		           const double* c, int_t incc, 
		                 double* y, int_t incy);
  void   svpoly  (int_t n, const float*  x, int_t incx, int_t m,
		           const float*  c, int_t incc, 
		                 float*  y, int_t incy);
  
  double dpoly   (int_t n, double x, const double* xp, const double* yp);
  float  spoly   (int_t n, float  x, const float*  xp, const float*  yp);
  void   dpolint (const double* xa, const double* ya, int_t n,
		  double x,        double* y,  double* dy);
  void   spolint (const float*  xa, const float*  ya, int_t n,
		  float  x,        float*  y, float*  dy);
  
  void   dspline  (int_t n, double yp1, double ypn,
		   const double* x, const double* y,  double* y2);
  void   sspline  (int_t n, float yp1, float ypn,
		   const float*  x, const float*  y,  float*  y2);
  double dsplint  (int_t n, double x,
		   const double* xa, const double* ya, const double* y2a);
  float  ssplint  (int_t n, float  x,
		   const float*  xa, const float*  ya,  const float*  y2a);
  double dsplquad (const double* xa, const double* ya,  const double* y2a,
		   const int_t n,    const double xmin, const double  xmax);
  float  ssplquad (const float*  xa, const float*  ya,  const float*  y2a,
		   const int_t   n,  const float  xmin, const float   xmax);

  void   dvvtvvtp (int_t n, const double* v, int_t incv,
          	            const double* w, int_t incw,
	                    const double* x, int_t incx,
	                    const double* y, int_t incy,
	                          double* z, int_t incz);
  void   svvtvvtp (int_t n, const float*  v, int_t incv,
	                    const float*  w, int_t incw,
	                    const float*  x, int_t incx,
	                    const float*  y, int_t incy,
	                          float*  z, int_t incz); 

  void   dvvtvvtm (int_t n, const double* v, int_t incv,
          	            const double* w, int_t incw,
	                    const double* x, int_t incx,
	                    const double* y, int_t incy,
	                          double* z, int_t incz);
  void   svvtvvtm (int_t n, const float*  v, int_t incv,
	                    const float*  w, int_t incw,
	                    const float*  x, int_t incx,
	                    const float*  y, int_t incy,
	                          float*  z, int_t incz); 

  void   dsvvttvp (int_t n, const double  alpha,
          	            const double* w, int_t incw,
	                    const double* x, int_t incx,
	                    const double* y, int_t incy,
	                          double* z, int_t incz);
  void   ssvvttvp (int_t n, const float   alpha,
	                    const float*  w, int_t incw,
	                    const float*  x, int_t incx,
	                    const float*  y, int_t incy,
	                          float*  z, int_t incz);  
}


class Veclib {
 public:

  // -- ZERO-OFFSET 2D MATRIX ADDRESSES:

  static inline int_t row_major (int_t i, int_t j, int_t n)
  { return j + i * n; }
  static inline int_t col_major (int_t i, int_t j, int_t n)
  { return i + j * n; }


  // -- MATHEMATICAL PRIMITIVES:

  static void copy (int_t n, const double* x, int_t incx,
                                   double* y, int_t incy)
  { dcopy (n, x, incx, y, incy); }
  static void copy (int_t n, const int_t*  x, int_t incx,
                                   int_t*  y, int_t incy)
  { icopy (n, x, incx, y, incy); }
  static void copy (int_t n, const float*  x, int_t incx,
                                   float*  y, int_t incy)
  { scopy (n, x, incx, y, incy); } 

  
  static void fill (int_t n, double alpha, double* x, int_t incx)
  { dfill (n, alpha, x, incx); }
  static void fill (int_t n, int_t  alpha, int_t*  x, int_t incx)
  { ifill (n, alpha, x, incx); }
  static void fill (int_t n, float  alpha, float*  x, int_t incx)
  { sfill (n, alpha, x, incx); }


  static void znan (int_t n, double* x, int_t incx)
  { dznan (n, x, incx); }
  static void znan (int_t n, float*  x, int_t incx)
  { sznan (n, x, incx); }

  static void neg (int_t n, double* x, int_t incx)
  { dneg (n, x, incx); }
  static void neg (int_t n, int_t*  x, int_t incx)
  { ineg (n, x, incx); }
  static void neg (int_t n, float*  x, int_t incx)
  { sneg (n, x, incx); }


  static void vneg (int_t n,
		    const double* x, int_t incx, double* y, int_t incy)
  { dvneg (n, x, incx, y, incy); }
  static void vneg (int_t n,
		    const int_t*  x, int_t incx, int_t*  y, int_t incy)
  { ivneg (n, x, incx, y, incy); }
  static void vneg (int_t n,
		    const float*  x, int_t incx, float*  y, int_t incy)
  { svneg (n, x, incx, y, incy); }


  static void vsgn (int_t n,
		    const double* x, int_t incx, double* y, int_t incy)
  { dvsgn (n, x, incx, y, incy); }
  static void vsgn (int_t n,
		    const int_t*  x, int_t incx, int_t*  y, int_t incy)
  { ivsgn (n, x, incx, y, incy); }
  static void vsgn (int_t n,
		    const float*  x, int_t incx, float*  y, int_t incy)
  { svsgn (n, x, incx, y, incy); }

  
  static void sadd (int_t n, double alpha, const double* x, int_t incx, 
		                                 double* y, int_t incy)
  { dsadd (n, alpha, x, incx, y, incy); }
  static void sadd (int_t n, int_t  alpha, const int_t*  x, int_t incx,
		                                 int_t*  y, int_t incy)
  { isadd (n, alpha, x, incx, y, incy); }
  static void sadd (int_t n, float  alpha, const float*  x, int_t incx,
		                                 float*  y, int_t incy)
  { ssadd (n, alpha, x, incx, y, incy); }


  static void spow (const int_t n, const double alpha,
		    const double* x, int_t incx, double* y, int_t incy)
  { dspow (n, alpha, x, incx, y, incy); }
  static void spow (const int_t n, const float alpha,
		    const float *x, int_t incx, float *y, int_t incy)
  { sspow (n, alpha, x, incx, y, incy); }


  static void vadd (int_t n, const double* x, int_t incx,
		             const double* y, int_t incy,
                                   double* z, int_t incz)
  { dvadd (n, x, incx, y, incy, z, incz); }
  static void vadd (int_t n, const int_t*  x, int_t incx,
		             const int_t*  y, int_t incy,
                                   int_t*  z, int_t incz)
  { ivadd (n, x, incx, y, incy, z, incz); }
  static void vadd (int_t n, const float*  x, int_t incx,
		             const float*  y, int_t incy,
                                   float*  z, int_t incz)
  { svadd (n, x, incx, y, incy, z, incz); }

 
  static void ssub (int_t n, double alpha, const double* x, int_t incx, 
		                                 double* y, int_t incy)
  { dssub (n, alpha, x, incx, y, incy); }
  static void ssub (int_t n, int_t  alpha, const int_t*  x, int_t incx,
		                                 int_t*  y, int_t incy)
  { issub (n, alpha, x, incx, y, incy); }
  static void ssub (int_t n, float  alpha, const float*  x, int_t incx,
		                                 float*  y, int_t incy)
  { sssub (n, alpha, x, incx, y, incy); }

  
  static void vsub (int_t n, const double* x, int_t incx,
		             const double* y, int_t incy,
                                   double* z, int_t incz)
  { dvsub (n, x, incx, y, incy, z, incz); }
  static void vsub (int_t n, const int_t*  x, int_t incx,
		             const int_t*  y, int_t incy,
                                   int_t*  z, int_t incz)
  { ivsub (n, x, incx, y, incy, z, incz); }
  static void vsub (int_t n, const float*  x, int_t incx,
		             const float*  y, int_t incy,
                                   float*  z, int_t incz)
  { svsub (n, x, incx, y, incy, z, incz); }

  
  static void  smul  (int_t n, double alpha, const double* x, int_t incx,
		                                   double* y, int_t incy)
  { dsmul (n, alpha, x, incx, y, incy); }
  static void  smul  (int_t n, int_t  alpha, const int_t*  x, int_t incx,
		                                   int_t*  y, int_t incy)
  { ismul (n, alpha, x, incx, y, incy); }
  static void  smul  (int_t n, float  alpha,  const float* x, int_t incx,
		                                   float*  y, int_t incy)
  { ssmul (n, alpha, x, incx, y, incy); }

  
  static void vmul (int_t n, const double* x, int_t incx,
		             const double* y, int_t incy,
                                   double* z, int_t incz)
  { dvmul (n, x, incx, y, incy, z, incz); } 
  static void vmul (int_t n, const int_t*  x, int_t incx,
		             const int_t*  y, int_t incy,
                                   int_t*  z, int_t incz)
  { ivmul (n, x, incx, y, incy, z, incz); }
  static void vmul (int_t n, const float*  x, int_t incx,
		             const float*  y, int_t incy,
                                   float*  z, int_t incz)
  { svmul (n, x, incx, y, incy, z, incz); }


  static void sdiv (int_t n, double alpha, const double* x, int_t incx,
		                                 double* y, int_t incy)
  { dsdiv (n, alpha, x, incx, y, incy); }
  static void sdiv (int_t n, int_t  alpha, const int_t*  x, int_t incx,
		                                 int_t*  y, int_t incy)
  { isdiv (n, alpha, x, incx, y, incy); }
  static void  sdiv (int_t n, float alpha, const float*  x, int_t incx,
		                                 float*  y, int_t incy)
  { ssdiv (n, alpha, x, incx, y, incy); }

  
  static void vrecp (int_t n, const double* x, int_t incx,
                                    double* y, int_t incy)
  { dvrecp (n, x, incx, y, incy); }
  static void vrecp (int_t n, const float*  x, int_t incx,
                                    float*  y, int_t incy)
  { svrecp (n, x, incx, y, incy); }
  

  static void vdiv (int_t n, const double* x, int_t incx,
		             const double* y, int_t incy,
		                   double* z, int_t incz)
  { dvdiv (n, x, incx, y, incy, z, incz); }
  static void vdiv (int_t n, const float*  x, int_t incx,
		             const float*  y, int_t incy,
		                   float*  z, int_t incz)
  { svdiv (n, x, incx, y, incy, z, incz); }

  
  static void zero (int_t n, double* x, int_t incx)
  { dzero (n, x, incx); }
  static void zero (int_t n, int_t*  x, int_t incx)
  { izero (n, x, incx); }
  static void zero (int_t n, float*  x, int_t incx)
  { szero (n, x, incx); }


  // -- OTHER MATHEMATICAL FUNCTIONS:
  
  static void vabs (int_t n, const double*  x, int_t incx,
                                   double*  y, int_t incy)
  { dvabs (n, x, incx, y, incy); }
  static void vabs (int_t n, const int_t* x, int_t incx,
                                   int_t* y, int_t incy)
  { ivabs (n, x, incx, y, incy); }
  static void vabs (int_t n, const float*   x, int_t incx,
                                   float*   y, int_t incy)
  { svabs (n, x, incx, y, incy); }


  static void vamax (int_t n, const double*  x, int_t incx,
		              const double*  y, int_t incy,
		                    double*  z, int_t incz)
  { dvamax (n, x, incx, y, incy, z, incz); }
  static void vamax (int_t n, const int_t* x, int_t incx,
		              const int_t* y, int_t incy,
                                    int_t* z, int_t incz)
  { ivamax (n, x, incx, y, incy, z, incz); }
  static void vamax (int_t n, const float*   x, int_t incx,
		              const float*   y, int_t incy,
                                    float*   z, int_t incz)
  { svamax (n, x, incx, y, incy, z, incz); }

  
  static void vexp (int_t n, const double* x, int_t incx,
                                   double* y, int_t incy)
  { dvexp (n, x, incx, y, incy); }
  static void vexp (int_t n, const float*  x, int_t incx,
		                   float*  y, int_t incy)
  { svexp (n, x, incx, y, incy); }

  
  static void vlg10 (int_t n, const double* x, int_t incx,
		                    double* y, int_t incy)
  { dvlg10 (n, x, incx, y, incy); }
  static void vlg10 (int_t n, const float*  x, int_t incx,
		                    float*  y, int_t incy)
  { svlg10 (n, x, incx, y, incy); }


  static void vlog (int_t n, const double* x, int_t incx,
                                   double* y, int_t incy)
  { dvlog (n, x, incx, y, incy); }
  static void vlog (int_t n, const float*  x, int_t incx,
                                   float*  y, int_t incy)
  { svlog (n, x, incx, y, incy); }


  static void vatan (int_t n, const double* x, int_t incx,
	                            double* y, int_t incy)
  { dvatan (n, x, incx, y, incy); }
  static void vatan (int_t n, const float*  x, int_t incx,
                                    float*  y, int_t incy)
  { svatan (n, x, incx, y, incy); } 


  static void vatn2 (int_t n, const double* x, int_t incx,
		              const double* y, int_t incy,
                                    double* z, int_t incz)
  { dvatn2 (n, x, incx, y, incy, z, incz); }
  static void vatn2 (int_t n, const float*  x, int_t incx,
		              const float*  y, int_t incy,
                                    float*  z, int_t incz)
  { svatn2 (n, x, incx, y, incy, z, incz); }


  static void vcos (int_t n, const double* x, int_t incx,
		                   double* y, int_t incy)
  { dvcos (n, x, incx, y, incy); }
  static void vcos (int_t n, const float*  x, int_t incx,
		                   float*  y, int_t incy)
  { svcos (n, x, incx, y, incy); }


  static void vsin (int_t n, const double* x, int_t incx,
		                   double* y, int_t incy)
  { dvsin (n, x, incx, y, incy); }
  static void vsin (int_t n, const float*  x, int_t incx,
		                   float*  y, int_t incy)
  { svsin (n, x, incx, y, incy); }


  static void vsqrt (int_t n, const double* x, int_t incx,
		                    double* y, int_t incy)
  { dvsqrt (n, x, incx, y, incy); }
  static void vsqrt (int_t n, const float*  x, int_t incx,
		                    float*  y, int_t incy)
  { svsqrt (n, x, incx, y, incy); }


  static void vtanh (int_t n, const double* x, int_t incx,
		                    double* y, int_t incy)
  { dvtanh (n, x, incx, y, incy); }
  static void vtanh (int_t n, const float*  x, int_t incx,
		                    float*  y, int_t incy)
  { svtanh (n, x, incx, y, incy); }


  static void ranInit (int_t flag)
  { raninit (flag); }
  static double dranu ()         
  { return dranu (); }
  static float  sranu ()         
  { return sranu (); }
  static double drang ()         
  { return drang (); }
  static float  srang ()         
  { return srang (); }
  static double normal (double mean, double sdev)
  { return dnormal (mean, sdev); }
  static float  normal (float  mean, float  sdev)
  { return snormal (mean, sdev); }
  static void  vrandom (int_t n, double* x, int_t incx)
  { dvrandom (n, x, incx); }
  static void  vrandom (int_t n, float*  x, int_t incx)
  { svrandom (n, x, incx); }
  static void  vgauss  (int_t n, double* x, int_t incx)
  { dvgauss (n, x, incx); }
  static void  vgauss (int_t n, float*   x, int_t incx)
  { svgauss (n, x, incx); }
  static void  vnormal (int_t n, double mean, double sdev,
			double* x, int_t incx)
  { dvnormal (n, mean, sdev, x, incx); }
  static void  vnormal (int_t n, float  mean, float  sdev,
			float*  x, int_t incx)
  { svnormal (n, mean, sdev, x, incx); }


  static void vhypot (int_t n, const double* x, int_t incx,
		               const double* y, int_t incy,
                                     double* z, int_t incz)
  { dvhypot (n, x, incx, y, incy, z, incz); }
  static void vhypot (int_t n, const float*  x, int_t incx,
		               const float*  y, int_t incy,
                                     float*  z, int_t incz)
  { svhypot (n, x, incx, y, incy, z, incz); }
  static void vmag (int_t n, const double* w, int_t incw,
		    const double* x, int_t incx,
		    const double* y, int_t incy,
		          double* z, int_t incz)
  { dvmag (n, w, incw, x, incx, y, incy, z, incz); }
  static void vmag (int_t n, const float*  w, int_t incw,
		    const float*  x, int_t incx,
		    const float*  y, int_t incy,
		          float*  z, int_t incz)
  { svmag (n, w, incw, x, incx, y, incy, z, incz); }

  static void vpow (int_t n, const double* x, int_t incx,
		    const double* y, int_t incy, double* z, int_t incz)
  { dvpow (n, x, incx, y, incy, z, incz); } 
  static void vpow (int_t n, const float*  x, int_t incx,
		    const float*  y, int_t incy, float*  z, int_t incz)
  { svpow (n, x, incx, y, incy, z, incz); } 


  // -- TRIAD OPERATIONS:
  
  static void svmvt (int_t n, double alpha, const double*  x, int_t incx,
		                            const double*  y, int_t incy,
		                                  double*  z, int_t incz)
  { dsvmvt (n, alpha, x, incx, y, incy, z, incz); }
  static void svmvt (int_t n, float  alpha, const float*  x, int_t incx,
		                            const float*  y, int_t incy,
                                                  float*  z, int_t incz)
  { ssvmvt (n, alpha, x, incx, y, incy, z, incz); }

 
  static void svpvt (int_t n, double alpha, const double* x, int_t incx,
		                            const double* y, int_t incy,
                                              double* z, int_t incz)
  { dsvpvt (n, alpha, x, incx, y, incy, z, incz); }
  static void svpvt (int_t n, float  alpha, const float*  x, int_t incx,
		                            const float*  y, int_t incy,
                                                  float*  z, int_t incz)
  { ssvpvt (n, alpha, x, incx, y, incy, z, incz); }


  static void svtsp (int_t n, double alpha, double beta,
		     const double* x, int_t incx,
		           double* y, int_t incy)
  { dsvtsp (n, alpha, beta, x, incx, y, incy); }
  static void svtsp (int_t n, float  alpha, float  beta,
		     const float*  x, int_t incx,
		           float*  y, int_t incy)
  { ssvtsp (n, alpha, beta, x, incx, y, incy); }
 
 
  static void svtvm (int_t n, double alpha, const double* x, int_t incx,
		                            const double* y, int_t incy,
		                                  double* z, int_t incz)
  { dsvtvm (n, alpha, x, incx, y, incy, z, incz); }
  static void svtvm (int_t n, float  alpha, const float*  x, int_t incx,
		                            const float*  y, int_t incy,
                                                  float*  z, int_t incz)
  { ssvtvm (n, alpha, x, incx, y, incy, z, incz); }


  static void svtvp (int_t n, double alpha, const double* x, int_t incx,
		                            const double* y, int_t incy,
                                                  double* z, int_t incz)
  { dsvtvp (n, alpha, x, incx, y, incy, z, incz); }
  static void svtvp (int_t n, float  alpha, const float*  x, int_t incx,
		                            const float*  y, int_t incy,
                                                  float*  z, int_t incz)
  { ssvtvp (n, alpha, x, incx, y, incy, z, incz); }

  
  static void svvmt (int_t n, double alpha, const double* x, int_t incx,
		                            const double* y, int_t incy,
                                                  double* z, int_t incz)
  { dsvvmt (n, alpha, x, incx, y, incy, z, incz); }
  static void svvmt (int_t n, float  alpha, const float*  x, int_t incx,
		                            const float*  y, int_t incy,
		                                  float*  z, int_t incz)
  { ssvvmt (n, alpha, x, incx, y, incy, z, incz); }

 
  static void svvpt (int_t n, double alpha, const double* x, int_t incx,
		                            const double* y, int_t incy,
		                                  double* z, int_t incz)
  { dsvvpt (n, alpha, x, incx, y, incy, z, incz); }
  static void svvpt (int_t n, float  alpha, const float*  x, int_t incx,
		                            const float*  y, int_t incy,
                                                  float*  z, int_t incz)
  { ssvvpt (n, alpha, x, incx, y, incy, z, incz); }


  static void svvtm (int_t n, double alpha, const double* x, int_t incx,
		                            const double* y, int_t incy,
		                                  double* z, int_t incz)
  { dsvvtm (n, alpha, x, incx, y, incy, z, incz); }
  static void svvtm (int_t n, float  alpha, const float*  x, int_t incx,
		                            const float*  y, int_t incy,
		                                  float*  z, int_t incz)
  { ssvvtm (n, alpha, x, incx, y, incy, z, incz); }

 
  static void svvtp (int_t n, double alpha, const double* x, int_t incx,
		                            const double* y, int_t incy,
		                                  double* z, int_t incz)
  { dsvvtp (n, alpha, x, incx, y, incy, z, incz); }
  static void svvtp (int_t n, float  alpha, const float*  x, int_t incx,
		                            const float*  y, int_t incy,
		                                  float*  z, int_t incz)
  { ssvvtp (n, alpha, x, incx, y, incy, z, incz); }


  static void svvtt (int_t n, double alpha, const double* x, int_t incx,
		                            const double* y, int_t incy,
		                                  double* z, int_t incz)
  { dsvvtt (n, alpha, x, incx, y, incy, z, incz); }
  static void svvtt (int_t n, float  alpha, const float*  x, int_t incx,
		                            const float*  y, int_t incy,
		                                  float*  z, int_t incz)
  { ssvvtt (n, alpha, x, incx, y, incy, z, incz); }


  static void vvmvt (int_t n, const double* w, int_t incw,
		              const double* x, int_t incx,
		              const double* y, int_t incy,
		                    double* z, int_t incz)
  { dvvmvt (n, w, incw, x, incx, y, incy, z, incz); }
  static void vvmvt (int_t n, const float*  w, int_t incw,
		              const float*  x, int_t incx,
		              const float*  y, int_t incy,
		                    float*  z, int_t incz)
  { svvmvt (n, w, incw, x, incx, y, incy, z, incz); }

  
  static void vvpvt (int_t n, const double* w, int_t incw,
		              const double* x, int_t incx,
		              const double* y, int_t incy,
                                    double* z, int_t incz)
  { dvvpvt (n, w, incw, x, incx, y, incy, z, incz); }
  static void vvpvt (int_t n, const float*  w, int_t incw,
		              const float*  x, int_t incx,
		              const float*  y, int_t incy,
		                    float*  z, int_t incz)
  { svvpvt (n, w, incw, x, incx, y, incy, z, incz); }

 
  static void vvtvp (int_t n, const double* w, int_t incw,
		              const double* x, int_t incx,
		              const double* y, int_t incy,
		                    double* z, int_t incz)
  { dvvtvp (n, w, incw, x, incx, y, incy, z, incz); }
  static void vvtvp (int_t n, const float*  w, int_t incw,
		              const float*  x, int_t incx,
		              const float*  y, int_t incy,
		                    float*  z, int_t incz)
  { svvtvp (n, w, incw, x, incx, y, incy, z, incz); }

 
  static void vvtvm (int_t n, const double* w, int_t incw,
		              const double* x, int_t incx,
		              const double* y, int_t incy,
		                    double* z, int_t incz)
  { dvvtvm (n, w, incw, x, incx, y, incy, z, incz); }
  static void vvtvm (int_t n, const float*  w, int_t incw,
		              const float*  x, int_t incx,
		              const float*  y, int_t incy,
		                    float*  z, int_t incz)
  { svvtvm (n, w, incw, x, incx, y, incy, z, incz); }

  
  static void vvvtm (int_t n, const double* w, int_t incw,
		              const double* x, int_t incx,
		              const double* y, int_t incy,
		                    double* z, int_t incz)
  { dvvvtm (n, w, incw, x, incx, y, incy, z, incz); }
  static void vvvtm (int_t n, const float*  w, int_t incw,
		              const float*  x, int_t incx,
		              const float*  y, int_t incy,
		                    float*  z, int_t incz)
  { svvvtm (n, w, incw, x, incx, y, incy, z, incz); }
  

  static void vvvtt (int_t n, const double* w, int_t incw,
		              const double* x, int_t incx,
		              const double* y, int_t incy,
		                    double* z, int_t incz)
  { dvvvtt (n, w, incw, x, incx, y, incy, z, incz); }
  static void vvvtt (int_t n, const float*  w, int_t incw,
		              const float*  x, int_t incx,
		              const float*  y, int_t incy,
		                    float*  z, int_t incz)
  { svvvtt (n, w, incw, x, incx, y, incy, z, incz); }

 
  // -- RELATIONAL PRIMITIVE OPERATIONS:
  
  static void seq (int_t n, int_t alpha, const int_t* x, int_t incx,
                                               int_t* y, int_t incy)
  { iseq (n, alpha, x, incx, y, incy); }


  static void sge (int_t n, double alpha, const double* x, int_t incx,
		                                int_t*  y, int_t incy)
  { dsge (n, alpha, x, incx, y, incy); }
  static void sge (int_t n, int_t  alpha, const int_t*  x, int_t incx,
                                                int_t*  y, int_t incy)
  { isge (n, alpha, x, incx, y, incy); }
  static void sge (int_t n, float  alpha, const float*  x, int_t incx,
		                                int_t*  y, int_t incy)
  { ssge (n, alpha, x, incx, y, incy); }


  static void sle (int_t n, double alpha, const double* x, int_t incx,
		                                int_t*  y, int_t incy)
  { dsle (n, alpha, x, incx, y, incy); }
  static void sle (int_t n, int_t  alpha, const int_t*  x, int_t incx,
                                                int_t*  y, int_t incy)
  { isle (n, alpha, x, incx, y, incy); }
  static void sle (int_t n, float  alpha, const float*  x, int_t incx,
		                                int_t*  y, int_t incy)
  { ssle (n, alpha, x, incx, y, incy); }


  static void slt (int_t n, double alpha, const double* x, int_t incx,
		                                int_t*  y, int_t incy)
  { dslt (n, alpha, x, incx, y, incy); }
  static void slt (int_t n, int_t  alpha, const int_t*  x, int_t incx,
		                                int_t*  y, int_t incy)
  { islt (n, alpha, x, incx, y, incy); }
  static void slt (int_t n, float  alpha, const float*  x, int_t incx,
		                                int_t*  y, int_t incy)
  { sslt (n, alpha, x, incx, y, incy); }


  static void sne (int_t n, double alpha, const double* x, int_t incx,
		                                int_t*  y, int_t incy)
  { dsne (n, alpha, x, incx, y, incy); }
  static void sne (int_t n, int_t  alpha, const int_t*  x, int_t incx,
		                                int_t*  y, int_t incy)
  { isne (n, alpha, x, incx, y, incy); }
  static void sne (int_t n, float  alpha, const float*  x, int_t incx,
		                                int_t*  y, int_t incy)
  { ssne (n, alpha, x, incx, y, incy); } 
 

  // -- REDUCTION FUNCTIONS:
  
  static double sum (int_t n, const double* x, int_t incx)
  { return dsum (n, x, incx); }
  static int_t  sum (int_t n, const int_t*  x, int_t incx)
  { return isum (n, x, incx); }
  static float  sum (int_t n, const float*  x, int_t incx)
  { return ssum (n, x, incx); }


  static int_t imax (int_t n, const double* x, int_t incx)
  { return idmax (n, x, incx); }
  static int_t imax (int_t n, const int_t*  x, int_t incx)
  { return iimax (n, x, incx); }
  static int_t imax (int_t n, const float*  x, int_t incx)
  { return ismax (n, x, incx); }


  static int_t imin (int_t n, const double* x, int_t incx)
  { return idmin (n, x, incx); }
  static int_t imin (int_t n, const int_t*  x, int_t incx)
  { return iimin (n, x, incx); }
  static int_t imin (int_t n, const float*  x, int_t incx)
  { return ismin (n, x, incx); }


  static int_t count (int_t n, const int_t *x, int_t incx)
  { return icount (n, x, incx); }
  static int_t first (int_t n, const int_t *x, int_t incx)
  { return ifirst (n, x, incx); }
  static int_t any   (int_t n, const int_t *x, int_t incx)
  { return lany  (n, x, incx); }


  static int_t same (int_t n, const int_t* x, int_t incx,
                              const int_t* y, int_t incy)
  { return lisame (n, x, incx, y, incy); }
  static int_t same (int_t n, const double*  x, int_t incx,
                              const double*  y, int_t incy)
  { return ldsame (n, x, incx, y, incy); }
  static int_t same (int_t n, const float*   x, int_t incx,
                              const float*   y, int_t incy)
  { return lssame (n, x, incx, y, incy); }

    
  // -- CONVERSION PRIMITIVES:
  
  static void vdble (int_t n, const float*  x, int_t incx,
		                    double* y, int_t incy)
  { vdble (n, x, incx, y, incy); }
  static void vsngl (int_t n, const double* x, int_t incx,
		                    float*  y, int_t incy)
  { vsngl (n, x, incx, y, incy); }


  static void vfloa (int_t n, const int_t*  x, int_t incx,
                                    double* y, int_t incy)
  { dvfloa (n, x, incx, y, incy); }
  static void vfloa (int_t n, const int_t* x, int_t incx,
                                    float* y, int_t incy)
  { svfloa (n, x, incx, y, incy); }


  static int_t testFormat ()
  { return iformat (); }
  static void describeFormat (char* s)
  { format (s); }
  static void brev (int_t n, const double*  x, int_t incx,
		                   double*  y, int_t incy)
  { dbrev (n, x, incx, y, incy); }
  static void brev (int_t n, const int_t* x, int_t incx,
                                   int_t* y, int_t incy)
  { ibrev (n, x, incx, y, incy); }
  static void brev (int_t n, const float*   x, int_t incx,
                                   float*   y, int_t incy)
  { sbrev (n, x, incx, y, incy); }


   // -- MISCELLANEOUS FUNCTIONS:

  static double clock ()
  { return dclock () ; }
  

  static void scatr (int_t n, const double* x, const int_t* y, double* z)
  { dscatr (n, x, y, z); }
  static void scatr (int_t n, const int_t*  x, const int_t* y, int_t*  z)
  { iscatr (n, x, y, z); }
  static void scatr (int_t n, const float*  x, const int_t* y, float*  z)
  { sscatr (n, x, y, z); }


  static void gathr (int_t n, const double* x, const int_t* y, double* z)
  { dgathr (n, x, y, z); }
  static void gathr (int_t n, const int_t*  x, const int_t* y, int_t*  z)
  { igathr (n, x, y, z); }
  static void gathr (int_t n, const float*  x, const int_t* y, float*  z)
  { sgathr (n, x, y, z); }


  static void gathr_scatr (int_t n, const double* w, const int_t*  x,
			            const int_t*  y,       double* z)
  { dgathr_scatr (n, w, x, y, z); }
  static void gathr_scatr (int_t n, const int_t* w, const int_t* x,
			            const int_t* y,       int_t* z)
  { igathr_scatr (n, w, x, y, z); }
  static void gathr_scatr (int_t n, const float* w, const int_t* x,
 			            const int_t* y,       float* z)
  { sgathr_scatr (n, w, x, y, z); }

  
  static void scatr_sum (int_t n,
			 const double* x, const int_t *y, double* z)
  { dscatr_sum (n, x, y, z); }
  static void scatr_sum (int_t n,
			 const int_t*  x, const int_t *y, int_t*  z)
  { iscatr_sum (n, x, y, z); }
  static void scatr_sum (int_t n,
			 const float*  x, const int_t *y, float*  z)
  { sscatr_sum (n, x, y, z); }


  static void gathr_sum (int_t n,
			 const double* x, const int_t *y, double* z)
  { dgathr_sum (n, x, y, z); }
  static void gathr_sum (int_t n,
			 const int_t*  x, const int_t *y, int_t*  z)
  { igathr_sum (n, x, y, z); }
  static void gathr_sum (int_t n,
			 const float*  x, const int_t *y, float*  z)
  { sgathr_sum (n, x, y, z); }


  static void gathr_scatr_sum (int_t n, const double*  w, const int_t* x,
			                const int_t* y,       double*  z)
  { dgathr_scatr_sum (n, w, x, y, z); }
  static void gathr_scatr_sum (int_t n, const int_t* w, const int_t* x,
			                const int_t* y,       int_t* z)
  { igathr_scatr_sum (n, w, x, y, z); }
  static void gathr_scatr_sum (int_t n, const float*   w, const int_t* x,
			                const int_t* y,       float   *z)
  { sgathr_scatr_sum (n, w, x, y, z); }


  static void ramp (int_t n, double  alpha, double beta,
		    double*  x, int_t incx)
  { dramp (n, alpha, beta, x, incx); }
  static void ramp (int_t n, int_t   alpha, int_t  beta,
		    int_t* x, int_t incx)
  { iramp (n, alpha, beta, x, incx); }
  static void ramp (int_t n, float   alpha, float  beta,
		    float*   x, int_t incx)
  { sramp (n, alpha, beta, x, incx); }

  static void clip (int_t n, const double alpha, const double beta,
		    const double* x, int_t incx,
 		          double* y, int_t incy)
  { dclip (n, alpha, beta, x, incx, y, incy); }
  static void clip (int_t n, const  int_t alpha, const int_t  beta,
		    const int_t* x, int_t incx,
		          int_t* y, int_t incy)
  { iclip (n, alpha, beta, x, incx, y, incy); }
  static void clip (int_t n, const  float alpha, const float  beta,
		    const float* x,  int_t incx,
		          float* y, int_t incy)
  { sclip (n, alpha, beta, x, incx, y, incy); }

  static void clipup (int_t n, const double alpha,
		      const double* x, int_t incx,
		            double* y, int_t incy)
  { dclipup (n, alpha, x, incx, y, incy); }
  static void clipup (int_t n, const int_t alpha,
		      const int_t* x, int_t incx,
		            int_t* y, int_t incy)
  { iclipup (n, alpha, x, incx, y, incy); }
  static void clipup (int_t n, const float alpha,
		      const float* x, int_t incx,
		            float* y, int_t incy)
  { sclipup (n, alpha, x, incx, y, incy); }

  static void clipdn (int_t n, const double alpha,
		      const double* x, int_t incx,
		            double* y, int_t incy)
  { dclipdn (n, alpha, x, incx, y, incy); }
  static void clipdn (int_t n, const  int_t alpha,
		      const int_t* x, int_t incx,
		            int_t* y, int_t incy)
  { iclipdn (n, alpha, x, incx, y, incy); }
  static void clipdn (int_t n, const float alpha,
		      const float* x,  int_t incx,
		            float* y, int_t incy)
  { sclipdn (n, alpha, x, incx, y, incy); }

  static void iclip (int_t n, const double alpha,  const double beta,
		     const double* x, int_t incx,
		           double* y, int_t incy)
  { diclip (n, alpha, beta, x, incx, y, incy); }
  static void iclip (int_t n, const  int_t alpha, const int_t   beta,
		     const int_t* x, int_t incx,
		           int_t* y, int_t incy)
  { iiclip (n, alpha, beta, x, incx, y, incy); }
  static void iclip (int_t n, const float alpha,   const float  beta,
		     const float* x,  int_t incx,
		           float* y,  int_t incy)
  { siclip (n, alpha, beta, x, incx, y, incy); }

  static void cndst (int_t n, const double*  x, int_t incx,
		              const int_t* y, int_t incy,
		                    double*  z, int_t incz)
  { dcndst (n, x, incx, y, incy, z, incz); }
  static void cndst (int_t n, const int_t* x, int_t incx,
		              const int_t* y, int_t incy,
		                    int_t* z, int_t incz)
  { icndst (n, x, incx, y, incy, z, incz); }
  static void cndst (int_t n, const float*   x, int_t incx,
		              const int_t* y, int_t incy,
		                    float*   z, int_t incz)
  { scndst (n, x, incx, y, incy, z, incz); }
 

  static void mask (int_t n, const double*  w, int_t incw,
		             const double*  x, int_t incx,
		             const int_t* y, int_t incy,
		                   double*  z, int_t incz)
  { dmask (n, w, incw, x, incx, y, incy, z, incz); }
  static void mask (int_t n, const int_t* w, int_t incw,
		             const int_t* x, int_t incx,
		             const int_t* y, int_t incy,
		                   int_t *z, int_t incz)
  { imask (n, w, incw, x, incx, y, incy, z, incz); }
  static void mask (int_t n, const float*   w, int_t incw,
		             const float*   x, int_t incx,
		             const int_t* y, int_t incy, 
		                   float*   z, int_t incz)
  { smask (n, w, incw, x, incx, y, incy, z, incz); }

  
  static void vpoly (int_t n, const double* x, int_t incx, int_t m,
		              const double* c, int_t incc, 
		                    double* y, int_t incy)
  { dvpoly (n, x, incx, m, c, incc, y, incy); }
  static void vpoly (int_t n, const float*  x, int_t incx, int_t m,
		              const float*  c, int_t incc, 
		                    float*  y, int_t incy)
  { svpoly (n, x, incx, m, c, incc, y, incy); }  


  static double poly (int_t n, double x, const double* xp, const double* yp)
  { return dpoly (n, x, xp, yp); }
  static float  poly (int_t n, float  x, const float*  xp, const float*  yp)
  { return spoly (n, x, xp, yp); }


  static void polint (const double* xa, const double* ya, int_t n,
		            double  x,        double* y,  double* dy)
  { dpolint (xa, ya, n, x, y, dy); }
  static void polint (const float*  xa, const float*  ya, int_t n,
		            float   x,        float*  y, float*  dy)
  { spolint (xa, ya, n, x, y, dy); }
  

  static void spline (int_t n, double yp1, double ypn,
		      const double* x,  const double* y, double* y2)
  { dspline (n, yp1, ypn, x,  y, y2); }
  static void spline (int_t n, float yp1, float ypn,
		      const float*  x, const float*  y,  float*  y2)
  { sspline (n, yp1, ypn, x,  y, y2); }


  static double splint (int_t n, double x,
			const  double* xa, const double* ya, const double* y2a)
  { return dsplint (n, x, xa, ya, y2a); }
  static float  splint (int_t n, float  x,
			const  float*  xa, const float*  ya, const float*  y2a)
  { return ssplint (n, x, xa, ya, y2a); }


  static double splquad (const double* xa, const double* ya, const double* y2a,
			 const int_t   n,  const double xmn, const double  xmx)
  { return dsplquad (xa, ya, y2a, n, xmn, xmx); }
  static float  splquad (const float*  xa, const float*  ya, const float*  y2a,
			 const int_t   n,  const float  xmn, const float   xmx)
  { return ssplquad (xa, ya, y2a, n, xmn, xmx); }


  static void vvtvvtp (int_t n, const double* v, int_t incv,
          	                const double* w, int_t incw,
	                        const double* x, int_t incx,
	                        const double* y, int_t incy,
	                              double* z, int_t incz)
  { dvvtvvtp (n, v, incv, w, incw, x, incx, y, incy, z, incz); }
  static void vvtvvtp (int_t n, const float*  v, int_t incv,
	                        const float*  w, int_t incw,
	                        const float*  x, int_t incx,
	                        const float*  y, int_t incy,
		                      float*  z, int_t incz)
  { svvtvvtp (n, v, incv, w, incw, x, incx, y, incy, z, incz); }

  static void vvtvvtm (int_t n, const double* v, int_t incv,
          	                const double* w, int_t incw,
	                        const double* x, int_t incx,
	                        const double* y, int_t incy,
	                              double* z, int_t incz)
  { dvvtvvtm (n, v, incv, w, incw, x, incx, y, incy, z, incz); }
  static void vvtvvtm (int_t n, const float*  v, int_t incv,
	                        const float*  w, int_t incw,
	                        const float*  x, int_t incx,
	                        const float*  y, int_t incy,
		                      float*  z, int_t incz)
  { svvtvvtm (n, v, incv, w, incw, x, incx, y, incy, z, incz); }


  static void svvttvp (int_t n, const double  alpha,
          	                const double* w, int_t incw,
	                        const double* x, int_t incx,
	                        const double* y, int_t incy,
	                              double* z, int_t incz)
  { dsvvttvp (n, alpha, w, incw, x, incx, y, incy, z, incz); }
  static void svvttvp (int_t n, const float   alpha,
		                const float*  w, int_t incw,
		                const float*  x, int_t incx,
		                const float*  y, int_t incy,
		                      float*  z,       int_t incz)
  { ssvvttvp (n, alpha, w, incw, x, incx, y, incy, z, incz); }

};

#endif
