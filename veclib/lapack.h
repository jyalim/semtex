#ifndef LAPACK_H
#define LAPACK_H
///////////////////////////////////////////////////////////////////////////////
// C++ header file to access LAPACK routines.
//
// Notes:
// 1.  Class Lapack simply serves to inline calls to appropriate FORTRAN
//     routines.  Overloading is used to resolve single or double precision.
// 2.  General matrices are assumed to have column-major (FORTRAN) ordering,
//     packed forms are described in LAPACK documentation.
//
// $Id: lapack.h,v 8.1 2015/04/20 11:14:19 hmb Exp $
///////////////////////////////////////////////////////////////////////////////

#include <cfemdef.h>


extern "C" {
  void F77NAME(dgetrf) (const int_t& m, const int_t& n, double* a,
			const int_t& lda, int_t* ipiv, int_t& info);
  void F77NAME(sgetrf) (const int_t& m, const int_t& n, float*  a,
			const int_t& lda, int_t* ipiv, int_t& info);
  void F77NAME(dgetrs) (const char* tr, const int_t& n, const int_t& nrhs,
			const double* a, const int_t& lda, int_t* ipiv,
			double* b, const int_t& ldb, int_t& info);
  void F77NAME(sgetrs) (const char* tr, const int_t& n, const int_t& nrhs,
			const float*  a, const int_t& lda, int_t* ipiv,
			float*  b, const int_t& ldb, int_t& info);
  void F77NAME(dgecon) (const char* norm, const int_t& n, const double* a,
			const double& anorm, double& rcond, double* rwrk,
			int_t* iwork, int_t& info);
  void F77NAME(sgecon) (const char* norm, const int_t& n, const float*  a,
			const float&  anorm, float&  rcond, float*  rwrk,
			int_t* iwork, int_t& info);
  void F77NAME(dgetri) (const int_t& n, double* a, const int_t& lda,
			int_t* ipiv,  
			double* work, const int_t& lwork, int_t& info);
  void F77NAME(sgetri) (const int_t& n, float*  a, const int_t& lda,
			int_t* ipiv,
			float*  work, const int_t& lwork, int_t& info);
  void F77NAME(dpbtrf) (const char* uplo, const int_t& n, const int_t& kd, 
			double* ab, const int_t& ldab, int_t& info);
  void F77NAME(spbtrf) (const char* uplo, const int_t& n, const int_t& kd,
			float*  ab, const int_t& ldab, int_t& info);
  void F77NAME(dpbcon) (const char* uplo, const int_t& n, const int_t& kd,
			const double* ab, const int_t& ldab,
			const double& anorm, double& rcond,
			double* rwrk, int_t* iwork, int_t& info);
  void F77NAME(spbcon) (const char* uplo, const int_t& n, const int_t& kd,
			const float*  ab, const int_t& ldab,
			const float&  anorm, float& rcond,
			float* rwrk, int_t* iwork, int_t& info);
  void F77NAME(dpbtrs) (const char* uplo, const int_t& n, 
			const int_t& kd, const int_t& nrhs,
			const double* ab, const int_t& ldab,
			double* b, const int_t& ldb, int_t& info);
  void F77NAME(spbtrs) (const char* uplo, const int_t& n,
			const int_t& kd, const int_t& nrhs,
			const float*  ab, const int_t& ldab,
			float*  b, const int_t& ldb, int_t& info);
  void F77NAME(dpptrf) (const char* uplo, const int_t& n,
			double* ap, int_t& info);
  void F77NAME(spptrf) (const char* uplo, const int_t& n, 
			float*  ap, int_t& info);
  void F77NAME(dppcon) (const char* uplo, const int_t& n, const double* ap,
			const double& anorm, double& rcond, double* rwrk,
			int_t* iwrk, int_t& info);
  void F77NAME(sppcon) (const char* uplo, const int_t& n, const float*  ap,
			const float&  anorm, float&  rcond, float*  rwrk,
			int_t* iwrk, int_t& info);
  void F77NAME(dpptrs) (const char* uplo, const int_t& n,
			const int_t& nrhs, const double* ap,
			double* b, const int_t& ldb, int_t& info);
  void F77NAME(spptrs) (const char* uplo, const int_t& n,
			const int_t& nrhs, const float*  ap, float*  b,
			const int_t& ldb, int_t& info);
  void F77NAME(dgesv)  (const int_t& n, const int_t& nrhs, double* a,
			const int_t& lda, int_t* ipiv, double* b,
			const int_t& ldb, int_t& info);
  void F77NAME(sgesv)  (const int_t& n, const int_t& nrhs, float* a,
			const int_t& lda, int_t* ipiv, float* b, 
			const int_t& ldb, int_t& info);
  void F77NAME(dgeev)  (const char* uplo, const char* lrev, const int_t& n,
			double* a, const int_t& lda, double* wr, double* wi,
			double* rev,  const int_t& ldr,
			double* lev,  const int_t& ldv,
			double* work, const int_t& lwork, int_t& info);
  void F77NAME(sgeev)  (const char* uplo, const char* lrev, const int_t& n,
			float*  a, const int_t& lda, float*  wr, float*  wi,
			float*  rev,  const int_t& ldr,
			float*  lev,  const int_t& ldv,
			float*  work, const int_t& lwork, int_t& info);
  void F77NAME(dspev)  (const char* jobz, const char* uplo, const int_t& n,
			double* ap, double* w, double* z, const int_t& ldz,
			double* work, int_t& info);
  void F77NAME(sspev)  (const char* jobz, const char* uplo, const int_t& n,
			float* ap, float* w, float* z, const int_t& ldz,
			float* work, int_t& info);
}


class Lapack {
public:

  // -- See LAPACK Users' Guide, \S\S 5.3.2 & 5.3.3 for 'U' storage schemes:

  static inline int_t band_addr (int_t i, int_t j, int_t bw)
  { return (j+1)*bw-j+i-1; }
  static inline int_t pack_addr (int_t i, int_t j)
  { return  (((j+1)*j)>>1)+i; }

  // -- L-U factor a general real matrix.

  static void getrf (const int_t& m, const int_t& n, double *a,
		     const int_t& lda, int_t *ipiv, int_t& info) {
    F77NAME(dgetrf) (m,n,a,lda,ipiv,info);
  }
  static void getrf (const int_t& m, const int_t& n, float  *a,
		     const int_t& lda, int_t *ipiv, int_t& info) {
    F77NAME(sgetrf) (m,n,a,lda,ipiv,info);
  }

  // -- Solve a system of equations from L-U factorization.

  static void getrs (const char* trans, const int_t& n, const int_t& nrhs,
		     const double* a, const int_t& lda, int_t* ipiv,
		     double* b, const int_t& ldb, int_t& info) {
    F77NAME(dgetrs) (trans,n,nrhs,a,lda,ipiv,b,ldb,info);
  }
  static void getrs (const char *trans, const int_t& n, const int_t& nrhs,
		     const float*  a, const int_t& lda, int_t* ipiv,
		     float*  b, const int_t& ldb, int_t& info) {
    F77NAME(sgetrs) (trans,n,nrhs,a,lda,ipiv,b,ldb,info);
  }
  // -- Estimate condition number of general matrix from factorization.

  static void gecon (const char* norm, const int_t& n, const double* a,
		     const double& anorm, double& rcond, double* rwrk,
		     int_t* iwork, int_t& info) {
    F77NAME(dgecon) (norm,n,a,anorm,rcond,rwrk,iwork,info);
  }
  static void gecon (const char* norm, const int_t& n, const float*  a,
		     const float&  anorm, float&  rcond, float*  rwrk,
		     int_t* iwork, int_t& info) {
    F77NAME(sgecon) (norm,n,a,anorm,rcond,rwrk,iwork,info);
  }

  // -- Invert a general real matrix from L-U factors.

  static void getri (const int_t& n, double *a, const int_t& lda,
		     int_t *ipiv, double *work, const int_t& lwork,
		     int_t& info) {
    F77NAME(dgetri) (n,a,lda,ipiv,work,lwork,info);
  }
  static void getri (const int_t& n, float  *a, const int_t& lda,
		     int_t *ipiv, float  *work, const int_t& lwork, 
		     int_t& info) {
    F77NAME(sgetri) (n,a,lda,ipiv,work,lwork,info);
  }

  // -- Cholesky factorize a real positive-definite band-symmetric matrix.

  static void pbtrf (const char *uplo, const int_t& n, const int_t& kd, 
		     double *ab, const int_t& ldab, int_t& info) {
    F77NAME(dpbtrf) (uplo,n,kd,ab,ldab,info);
  }
  static void pbtrf (const char *uplo, const int_t& n, const int_t& kd,
		     float  *ab, const int_t& ldab, int_t& info) {
    F77NAME(spbtrf) (uplo,n,kd,ab,ldab,info);
  }

  // -- Estimate condition number of a real P-D B-S matrix from factorization.

  static void pbcon (const char* uplo, const int_t& n, const int_t& kd,
		     const double* ab, const int_t& ldab,
		     const double& anorm, double& rcond,
		     double* rwrk, int_t* iwork, int_t& info) {
    F77NAME(dpbcon) (uplo,n,kd,ab,ldab,anorm,rcond,rwrk,iwork,info);
  }
  static void pbcon (const char* uplo, const int_t& n, const int_t& kd,
		     const float*  ab, const int_t& ldab,
		     const float&  anorm, float&  rcond,
		     float*  rwrk, int_t* iwork, int_t& info) {
    F77NAME(spbcon) (uplo,n,kd,ab,ldab,anorm,rcond,rwrk,iwork,info);
  }

  // -- Solve a real, P-D B-S matrix problem using Cholesky factorization.

  static void pbtrs (const char *uplo, const int_t& n,
		     const int_t& kd, const int_t& nrhs,
		     const double *ab, const int_t& ldab,
		     double *b, const int_t& ldb, int_t& info) {
    F77NAME(dpbtrs) (uplo,n,kd,nrhs,ab,ldab,b,ldb,info);
  }
  static void pbtrs (const char *uplo, const int_t& n,
		     const int_t& kd, const int_t& nrhs,
		     const float  *ab, const int_t& ldab,
		     float  *b, const int_t& ldb, int_t& info) {
    F77NAME(spbtrs) (uplo,n,kd,nrhs,ab,ldab,b,ldb,info);
  }

  // -- Cholesky factor a real P-D packed-symmetric matrix.

  static void pptrf (const char *uplo, const int_t& n,
		     double *ap, int_t& info) {
    F77NAME(dpptrf) (uplo,n,ap,info);
  }
  static void pptrf (const char *uplo, const int_t& n,
		     float  *ap, int_t& info) {
    F77NAME(spptrf) (uplo,n,ap,info);
  }  

  // -- Estimate condition number of a P-D P-S matrix from factorization.

  static void ppcon (const char* uplo, const int_t& n, const double* ap,
		     const double& anorm, double& rcond, double* rwrk,
		     int_t* iwrk, int_t& info) {
    F77NAME(dppcon) (uplo,n,ap,anorm,rcond,rwrk,iwrk,info);
  }
  static void ppcon (const char* uplo, const int_t& n, const float*  ap,
		     const float&  anorm, float&  rcond, float*  rwrk,
		     int_t* iwrk, int_t& info) {
    F77NAME(sppcon) (uplo,n,ap,anorm,rcond,rwrk,iwrk,info);
  }

  // -- Solve a real P-D P-S matrix problem using Cholesky factorization.

  static void pptrs (const char *uplo, const int_t& n, const int_t& nrhs,
		     const double *ap, double *b, const int_t& ldb,
		     int_t& info) {
    F77NAME(dpptrs) (uplo,n,nrhs,ap,b,ldb,info);
  }
  static void pptrs (const char *uplo, const int_t& n, const int_t& nrhs,
		     const float  *ap, float  *b, const int_t& ldb,
		     int_t& info) {
    F77NAME(spptrs) (uplo,n,nrhs,ap,b,ldb,info);
  }

  // -- Solve a general system of real equations.

  static void gesv (const int_t& n, const int_t& nrhs,
		    double* a, const int_t& lda, int_t* ipiv,
		    double* b, const int_t& ldb, int_t& info) {
    F77NAME(dgesv) (n,nrhs,a,lda,ipiv,b,ldb,info);
  }
  static void gesv (const int_t& n, const int_t& nrhs,
		    float*  a, const int_t& lda, int_t* ipiv,
		    float*  b, const int_t& ldb, int_t& info) {
    F77NAME(sgesv) (n,nrhs,a,lda,ipiv,b,ldb,info);
  }

  // -- Solve general real matrix eigenproblem.

  static void geev (const char* uplo, const char* lrev, const int_t& n,
		    double* a, const int_t& lda, double* wr, double* wi,
		    double* rev,  const int_t& ldr,
		    double* lev,  const int_t& ldv,
		    double* work, const int_t& lwork, int_t& info) {
    F77NAME(dgeev) (uplo, lrev, n, a, lda, wr, wi, rev,
			 ldr, lev, ldv, work, lwork, info);
  }
  static void geev (const char* uplo, const char* lrev, const int_t& n,
		    float*  a, const int_t& lda, float*  wr, float*  wi,
		    float*  rev,  const int_t& ldr,
		    float*  lev,  const int_t& ldv,
		    float*  work, const int_t& lwork, int_t& info) {
    F77NAME(sgeev) (uplo, lrev, n, a, lda, wr, wi, rev,
			 ldr, lev, ldv, work, lwork, info);
  }

  // -- Solve packed-symmetric real matrix eigenproblem.

  static void spev (const char* jobz, const char* uplo, const int_t& n,
		    double* ap, double* w, double* z, const int_t& ldz,
		    double* work, int_t& info) {
    F77NAME(dspev) (jobz, uplo, n, ap, w, z, ldz, work, info);
  }
  static void spev (const char* jobz, const char* uplo, const int_t& n,
		    float* ap, float* w, float* z, const int_t& ldz,
		    float* work, int_t& info) {
    F77NAME(sspev) (jobz, uplo, n, ap, w, z, ldz, work, info);
  }

};

#endif
