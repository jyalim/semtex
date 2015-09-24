#ifndef BLAS_H
#define BLAS_H
///////////////////////////////////////////////////////////////////////////////
// C++ header file to access BLAS routines.
//
// Notes:
// 1.  Class Blas simply serves to inline calls to appropriate FORTRAN
//     routines.  Overloading is used to resolve single or double precision.
// 2.  iamax, which aliases calls to BLAS1 idamax & isamax, reduces the
//     values returned by those routines by 1 to conform with C++ arrays.
// 3.  mxv computes a matrix-vector product for a contiguously-stored row-
//     major matrix and vectors with unity skips; it aliases calls to
//     BLAS2 dgemv & sgemv.
// 4.  mxm computes matrix-matrix product for contiguously-stored row-major
//     matrices; it aliases calls to BLAS3 dgemm & sgemm.  Conforming assumed.
//
// $Id: blas.h,v 8.1 2015/04/20 11:14:19 hmb Exp $
///////////////////////////////////////////////////////////////////////////////

#include <cfemdef.h>

extern "C" {
  double  F77NAME(ddot)  (const int_t& n,
			  const double* x, const int_t& incx,
			  const double* y, const int_t& incy);
  double  F77NAME(dasum) (const int_t& n,
			  const double* x, const int_t& incx);
  double  F77NAME(dnrm2) (const int_t& n,
			  const double* x, const int_t& incx);
  void    F77NAME(drotg) (const double& a, const double& b,
			  double& c, double& s);
  void    F77NAME(drot)  (const int_t& n, double* x, const int_t& incx,
			  double* y, const int_t& incy,
			  const double& c, const double& s);
  void    F77NAME(dswap) (const int_t& n, double* x, const int_t& incx,
			  double* y, const int_t& incy);
  void    F77NAME(dscal) (const int_t& n, const double& alpha, 
			  double* x, const int_t& incx);
  void    F77NAME(daxpy) (const int_t& n, const double& alpha,
			  const double* x, const int_t& incx, 
			  double* y, const int_t& incy);
  int_t F77NAME(idamax)(const int_t& n, const double* x,
			  const int_t& incx);
  float   F77NAME(sdot)  (const int_t& n,
			  const float* x, const int_t& incx,
			  const float* y, const int_t& incy);
  float   F77NAME(sasum) (const int_t& n,
			  const float* x, const int_t& incx);
  float   F77NAME(snrm2) (const int_t& n,
			  const float* x, const int_t& incx);
  void    F77NAME(srotg) (const float& a, const float& b, float& c, float& s);
  void    F77NAME(srot)  (const int_t& n, float* x, const int_t& incx,
			  float* y, const int_t& incy,
			  const float& c, const float& s);
  void    F77NAME(sswap) (const int_t& n, float* x, const int_t& incx,
			  float* y, const int_t& incy);
  void    F77NAME(sscal) (const int_t& n, const float& alpha,
			  float* x, const int_t& incx);
  void    F77NAME(saxpy) (const int_t& n, const float& alpha,
			  const float* x, const int_t& incx, 
			  float* y, const int_t& incy);
  int_t F77NAME(isamax)(const int_t& n,
			  const float* x, const int_t& incx);
  void    F77NAME(dgemv) (const char* trans, const int_t& m,
			  const int_t& n, const double& alpha,
			  const double* a, const int_t& lda,
			  const double* x, const int_t& incx,
			  const double& beta, double* y, const int_t& incy);
  void    F77NAME(sgemv) (const char* trans, const int_t& m,
			  const int_t& n, const float& alpha, 
			  const float* a, const int_t& lda,
			  const float* x, const int_t& incx,
			  const float& beta, float* y, const int_t& incy);
  void    F77NAME(dger)  (const int_t& m, const int_t& n,
			  const double& alpha,
			  const double* x, const int_t& incx,
			  const double* y, const int_t& incy,
			  double* a, const int_t& lda);
  void    F77NAME(sger)  (const int_t& m, const int_t& n,
			  const float& alpha,
			  const float* x, const int_t& incx,
			  const float* y, const int_t& incy,
			  float* a, const int_t& lda);
  void    F77NAME(dspmv) (const char* uplo, const int_t& n,
			  const double& alpha, const double* ap,
			  const double* x, const int_t& incx,
			  const double& beta, double* y, const int_t& incy);
  void    F77NAME(sspmv) (const char* uplo, const int_t& n,
			  const float& alpha, const float* ap,
			  const float* x, const int_t& incx,
			  const float& beta, float* y, const int_t& incy);
  void    F77NAME(dgemm) (const char* ta, const char* tb,
			  const int_t& m, const int_t& n,
			  const int_t& k, const double& alpha,
			  const double* a, const int_t& lda,
			  const double* b, const int_t& ldb,
			  const double& beta, double* c, const int_t& ldc);
  void   F77NAME(sgemm)  (const char* ta, const char* tb, 
			  const int_t& m, const int_t& n,
			  const int_t& k, const float& alpha,
			  const float* a, const int_t& lda,
			  const float* b, const int_t& ldb, 
			  const float& beta, float* c, const int_t& ldc);
  void   F77NAME(dmxm)   (const double* a, const int_t& nra,
			  const double* b, const int_t& nca,
			        double* c, const int_t& ncb);
  void   F77NAME(smxm)   (const float*  a, const int_t& nra,
			  const float*  b, const int_t& nca,
			        float*  c, const int_t& ncb);
  void   F77NAME(dmxma)  (const double* a, const int_t& nra,
			  const double* b, const int_t& nca,
			        double* c, const int_t& ncb);
  void   F77NAME(smxma)  (const float*  a, const int_t& nra,
			  const float*  b, const int_t& nca,
			        float*  c, const int_t& ncb);
  void   F77NAME(dmxms)  (const double* a, const int_t& nra,
			  const double* b, const int_t& nca,
			        double* c, const int_t& ncb);
  void   F77NAME(smxms)  (const float*  a, const int_t& nra,
			  const float*  b, const int_t& nca,
			        float*  c, const int_t& ncb);
  void   F77NAME(dmxmts) (const double* a, const int_t& nra,
			  const double* b, const int_t& nca,
			        double* c, const int_t& ncb);
  void   F77NAME(smxmts) (const float*  a, const int_t& nra,
			  const float*  b, const int_t& nca,
			        float*  c, const int_t& ncb);
  void   F77NAME(dmxv)   (const double* A, const int_t& nra,
			  const double* x, const int_t& nca, double* y);
  void   F77NAME(smxv)   (const float* A, const int_t& nra,
			  const float* x, const int_t& nca, float* y);
}


class Blas {
public:

  // -- Level 1:

  static double dot (const int_t& n, const double* x, const int_t& incx,
		     const double* y, const int_t& incy) {
    return F77NAME(ddot) (n,x,incx,y,incy);
  }
  static float  dot (const int_t& n, const float* x, const int_t& incx,
		     const float* y, const int_t& incy) {
    return F77NAME(sdot) (n,x,incx,y,incy);
  }


  static double asum (const int_t& n, const double* x, const int_t& incx) {
    return F77NAME(dasum) (n,x,incx);
  }
  static float  asum (const int_t& n, const float* x, const int_t& incx) {
    return F77NAME(sasum) (n,x,incx);
  }


  static double nrm2 (const int_t& n, const double* x, const int_t& incx) {
    return F77NAME(dnrm2) (n,x,incx);
  }
  static float  nrm2 (const int_t& n, const float* x, const int_t& incx) {
    return F77NAME(snrm2) (n,x,incx);
  }


  static void rotg (const double& a, const double& b,
		          double& c,       double& s) {
    F77NAME(drotg) (a,b,c,s);
  }
  static void rotg (const float&  a, const float&  b,
		          float&  c,       float&  s) {
    F77NAME(srotg) (a,b,c,s);
  }


  static void rot  (const int_t& n, double* x, const int_t& incx,
		    double* y, const int_t& incy,
		    const double& c, const double& s        ) {
    F77NAME(drot) (n,x,incx,y,incy,c,s);
  }
  static void rot  (const int_t& n, float* x, const int_t& incx,
		    float* y, const int_t& incy,
		    const float& c, const float& s        ) {
    F77NAME(srot) (n,x,incx,y,incy,c,s);
  }


  static void swap (const int_t& n, double* x, const int_t& incx,
		    double* y, const int_t& incy) {
    F77NAME(dswap) (n,x,incx,y,incy);
  }
  static void swap (const int_t& n, float* x, const int_t& incx,
		    float* y, const int_t& incy) {
    F77NAME(sswap) (n,x,incx,y,incy);
  }


  static void scal (const int_t& n, const double& alpha, 
		    double* x, const int_t& incx) {
    F77NAME(dscal) (n,alpha,x,incx);
  }
  static void scal (const int_t& n, const float& alpha, 
		    float* x, const int_t& incx) {
    F77NAME(sscal) (n,alpha,x,incx);
  }


  static void axpy (const int_t& n, const double& alpha,
		    const double* x, const int_t& incx, 
		    double* y, const int_t& incy) {
    F77NAME(daxpy) (n,alpha,x,incx,y,incy);
  }
  static void axpy (const int_t& n, const float& alpha,
		    const float* x, const int_t& incx, 
		    float* y, const int_t& incy) {
    F77NAME(saxpy) (n,alpha,x,incx,y,incy);
  }


  static int_t iamax (const int_t& n, const double* x, const int_t& incx)
  {
    return F77NAME(idamax) (n,x,incx) - 1;
  }
  static int_t iamax (const int_t& n, const float *x, const int_t& incx)
  {
    return F77NAME(isamax) (n,x,incx) - 1;
  }

  // -- BLAS level 2:

  static void gemv (const char* trans, const int_t& m, const int_t& n,
		    const double& alpha, const double* a, const int_t& lda,
		    const double* x, const int_t& incx, const double& beta,
		    double* y, const int_t& incy) {
    F77NAME(dgemv) (trans,m,n,alpha,a,lda,x,incx,beta,y,incy);
  }
  static void gemv (const char* trans, const int_t& m, const int_t& n,
		    const float&  alpha, const float* a, const int_t& lda,
		    const float* x, const int_t& incx, const float&  beta,
		    float* y, const int_t& incy) {
    F77NAME(sgemv) (trans,m,n,alpha,a,lda,x,incx,beta,y,incy);
  }


  static void ger (const int_t& m, const int_t& n, const double& alpha,
		   const double* x, const int_t& incx,
		   const double* y, const int_t& incy,
		   double* a, const int_t& lda) {
    F77NAME(dger) (m,n,alpha,x,incx,y,incy,a,lda);
  }
  static void ger (const int_t& m, const int_t& n, const float&  alpha,
		   const float* x, const int_t& incx,
		   const float* y, const int_t& incy,
		   float* a, const int_t& lda) {
    F77NAME(sger) (m,n,alpha,x,incx,y,incy,a,lda);
  }


  static void spmv (const char* uplo, const int_t& n, const double& alpha,
		    const double* ap, const double* x, const int_t& incx,
		    const double& beta, double* y, const int_t& incy) {
    F77NAME(dspmv) (uplo,n,alpha,ap,x,incx,beta,y,incy);
  }
  static void spmv (const char* uplo, const int_t& n, const float& alpha,
		    const float* ap, const float* x, const int_t& incx,
		    const float&  beta, float* y, const int_t& incy) {
    F77NAME(sspmv) (uplo,n,alpha,ap,x,incx,beta,y,incy);
  }

#if 0
  static void mxv (const double* A, const int_t& nra, const double* x,
		   const int_t& nca, double* y) {
    F77NAME(dmxv) (A,nra,x,nca,y);
  }
  static void mxv (const float* A, const int_t& nra, const float* x,
		   const int_t& nca, float* y) {
    F77NAME(smxv) (A,nra,x,nca,y);
  }
#else
  static void mxv (const double* A, const int_t& nra, const double* x,
		   const int_t& nca, double* y) {
    F77NAME(dgemv) ("T",nca,nra,1.0,A,nca,x,1,0.0,y,1);
  }
  static void mxv (const float* A, const int_t& nra, const float* x,
		   const int_t& nca, float* y) {
    F77NAME(sgemv) ("T",nca,nra,1.0,A,nca,x,1,0.0,y,1);
  }
#endif

  // -- BLAS level 3:

  static void gemm (const char* ta, const char* tb,
		    const int_t& m, const int_t& n, const int_t& k,
		    const double& alpha, const double* a, const int_t& lda,
		    const double* b, const int_t& ldb, const double& beta,
		    double* c, const int_t& ldc) {
    F77NAME(dgemm) (ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
  }
  static void gemm (const char* ta, const char* tb, 
		    const int_t& m, const int_t& n, const int_t& k,
		    const float&  alpha, const float* a, const int_t& lda,
		    const float* b, const int_t& ldb, const float&  beta,
		    float* c, const int_t& ldc) {
    F77NAME(sgemm) (ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
  }

#if defined (_SX)

  static void mxm (const double* A, const int_t& nra, const double* B,
		   const int_t& nca, double* C, const int_t& ncb) {
    F77NAME(dmxm) (A, nra, B, nca, C, ncb);
  }
  static void mxm (const float* A, const int_t& nra, const float* B,
		   const int_t& nca, float* C, const int_t& ncb) {
    F77NAME(smxm) (A, nra, B, nca, C, ncb);
  }

  static void mxma (const double* A, const int_t& nra, const double* B,
		    const int_t& nca, double* C, const int_t& ncb) {
    F77NAME(dmxma) (A, nra, B, nca, C, ncb);
  }
  static void mxma (const float* A, const int_t& nra, const float* B,
		    const int_t& nca, float* C, const int_t& ncb) {
    F77NAME(smxma) (A, nra, B, nca, C, ncb);
  }

  static void mxms (const double* A, const int_t& nra, const double* B,
		    const int_t& nca, double* C, const int_t& ncb) {
    F77NAME(dmxms) (A, nra, B, nca, C, ncb);
  }
  static void mxms (const float* A, const int_t& nra, const float* B,
		    const int_t& nca, float* C, const int_t& ncb) {
    F77NAME(smxms) (A, nra, B, nca, C, ncb);
  }

  static void mxmts(const double* A, const int_t& nra, const double* B,
		    const int_t& nca, double* C, const int_t& ncbt) {
    F77NAME(dmxmts)(A, nra, B, nca, C, ncbt);
  }
  static void mxmts(const float* A, const int_t& nra, const float* B,
		    const int_t& nca, float* C, const int_t& ncbt) {
    F77NAME(smxmts)(A, nra, B, nca, C, ncbt);
  }
#else
  static void mxm (const double* A, const int_t& nra, const double* B,
		   const int_t& nca, double* C, const int_t& ncb) {
    F77NAME(dgemm) ("N","N",ncb,nra,nca,1.0,B,ncb,A,nca,0.0,C,ncb);

  }
  static void mxm (const float* A, const int_t& nra, const float* B,
		   const int_t& nca, float* C, const int_t& ncb) {
    F77NAME(sgemm) ("N","N",ncb,nra,nca,1.0,B,ncb,A,nca,0.0,C,ncb);
  }

  static void mxma (const double* A, const int_t& nra, const double* B,
		    const int_t& nca, double* C, const int_t& ncb) {
    F77NAME(dgemm) ("N","N",ncb,nra,nca,1.0,B,ncb,A,nca,1.0,C,ncb);

  }
  static void mxma (const float* A, const int_t& nra, const float* B,
		    const int_t& nca, float* C, const int_t& ncb) {
    F77NAME(sgemm) ("N","N",ncb,nra,nca,1.0,B,ncb,A,nca,1.0,C,ncb);
  }

  static void mxms (const double* A, const int_t& nra, const double* B,
		    const int_t& nca, double* C, const int_t& ncb) {
    F77NAME(dgemm) ("N","N",ncb,nra,nca,1.0,B,ncb,A,nca,-1.0,C,ncb);

  }
  static void mxms (const float* A, const int_t& nra, const float* B,
		    const int_t& nca, float* C, const int_t& ncb) {
    F77NAME(sgemm) ("N","N",ncb,nra,nca,1.0,B,ncb,A,nca,-1.0,C,ncb);
  }

  static void mxmts(const double* A, const int_t& nra, const double* B,
		    const int_t& nca, double* C, const int_t& ncbt) {
    F77NAME(dgemm) ("T","N",ncbt,nra,nca,-1.0,B,nca,A,nca,1.0,C,ncbt);
  }
  static void mxmts(const float* A, const int_t& nra, const float* B,
		    const int_t& nca, float* C, const int_t& ncbt) {
    F77NAME(sgemm) ("T","N",ncbt,nra,nca,-1.0,B,nca,A,nca,1.0,C,ncbt);
  }
#endif
};

#endif
