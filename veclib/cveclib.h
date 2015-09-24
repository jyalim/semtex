#ifndef VECLIB_H
#define VECLIB_H
/*****************************************************************************
 *                           V E C L I B . H
 *
 * $Id: cveclib.h,v 8.1 2015/04/20 11:14:19 hmb Exp $
 *****************************************************************************/

#include <stdarg.h>

#include <cfemdef.h>

#ifndef BLAS			/* Assume vendor-supplied FORTRAN BLAS       */
#define BLAS 3			/* library will be available.  If not,       */
#endif                          /* define BLAS to be 0 before this point.    */

#ifndef LAPACK			/* Assume vendor-supplied FORTRAN LAPACK     */
#define LAPACK 1	        /* library will be available.  If not,       */
#endif                          /* define LAPACK to be 0 before this point.  */

#define NVREG  8		/* Registers reserved for FORTRAN linkage.   */

extern int_t _vecIreg[NVREG];
extern char    _vecCreg[NVREG];
extern float   _vecSreg[NVREG];
extern double  _vecDreg[NVREG];

#define SQR(a)     ((a) * (a))
#define SGN(a)     ((a) < 0.0 ?    -1.0 : 1.0)
#define MIN(a, b)  ((a) < (b) ?     (a) : (b))
#define MAX(a, b)  ((a) > (b) ?     (a) : (b))
#define SIGN(a, b) ((b) > 0.0 ? fabs(a) : -fabs(a))

#ifndef STR_MAX
#define STR_MAX 2048
#endif

#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif

#ifndef TWOPI
#define TWOPI  6.28318530717958647692
#endif

#ifndef PI_180
#define PI_180 0.01745329251994329576
#endif

#define EPSm3   1.0e-3
#define EPSm4   1.0e-4
#define EPSm5   1.0e-5
#define EPSm6   1.0e-6
#define EPSSP   6.0e-7
#define EPSm7   1.0e-7
#define EPSm8   1.0e-8
#define EPSm12  1.0e-12
#define EPSDP   6.0e-14
#define EPSm14  1.0e-14
#define EPSm20  1.0e-20
#define EPSm30  1.0e-30
 
typedef struct  {float  Re, Im;}         complex;
typedef struct  {double Re, Im;}         zomplex;
enum    err_lev {WARNING, ERROR, REMARK};


/* ------------------------------------------------------------------------- *
 * UTILITIES:
 * ------------------------------------------------------------------------- */

char buf[STR_MAX];

void    message (const char *routine, const char *txt, int level);
FILE*   efopen  (const char *file,    const char *mode);
double  dclock  (void);
float   sclock  (void);

void    printDvector (FILE *fp,  int_t width, int_t prec,
		      int_t ntot,  int_t inc,   int_t nfield, ...);
void    printIvector (FILE *fp,  int_t width,           
		      int_t ntot,  int_t inc,   int_t nfield, ...);
void    printSvector (FILE *fp,  int_t width, int_t prec,
		       int_t ntot,  int_t inc,   int_t nfield, ...);

/* ------------------------------------------------------------------------- *
 * MEMORY MANAGEMENT:
 * ------------------------------------------------------------------------- */

#define tempVector(v, n) double *v = dvector (0, n-1)
#define freeVector(v)    free (v)

complex   *cvector  (int_t nl, int_t nh);
complex  **cmatrix  (int_t rl, int_t rh, int_t cl, int_t ch);
complex ***c3matrix (int_t rl, int_t rh, int_t cl, int_t ch,
		     int_t dl, int_t dh);

double    *dvector  (int_t nl, int_t nh);
double   **dmatrix  (int_t rl, int_t rh, int_t cl, int_t ch);
double  ***d3matrix (int_t rl, int_t rh, int_t cl, int_t ch,
		     int_t dl, int_t dh);

float     *svector  (int_t nl, int_t nh);
float    **smatrix  (int_t rl, int_t rh, int_t cl, int_t ch);
float   ***s3matrix (int_t rl, int_t rh, int_t cl, int_t ch,
		     int_t dl, int_t dh);

int_t   *ivector  (int_t nl, int_t nh);
int_t  **imatrix  (int_t rl, int_t rh, int_t cl, int_t ch);
int_t ***i3matrix (int_t rl, int_t rh, int_t cl, int_t ch,
		     int_t dl, int_t dh);

zomplex   *zvector  (int_t nl, int_t nh);
zomplex  **zmatrix  (int_t rl, int_t rh, int_t cl, int_t ch);
zomplex ***z3matrix (int_t rl, int_t rh, int_t cl, int_t ch,
		     int_t dl, int_t dh);

void freeCvector  (complex   *v, int_t nl);
void freeCmatrix  (complex  **m, int_t nrl, int_t ncl);
void freeC3matrix (complex ***t, int_t nrl, int_t ncl, int_t ndl);

void freeDvector  (double    *v, int_t nl);
void freeDmatrix  (double   **m, int_t nrl, int_t ncl);
void freeD3matrix (double  ***t, int_t nrl, int_t ncl, int_t ndl);

void freeSvector  (float     *v, int_t nl);
void freeSmatrix  (float    **m, int_t nrl, int_t ncl);
void freeS3matrix (float   ***t, int_t nrl, int_t ncl, int_t ndl);

void freeIvector  (int_t   *v, int_t nl);
void freeImatrix  (int_t  **m, int_t nrl, int_t ncl);
void freeI3matrix (int_t ***t, int_t nrl, int_t ncl, int_t ndl);

void freeZvector  (zomplex   *v, int_t nl);
void freeZmatrix  (zomplex  **m, int_t nrl, int_t ncl);
void freeZ3matrix (zomplex ***t, int_t nrl, int_t ncl, int_t ndl);


/* ------------------------------------------------------------------------- *
 * MATHEMATICAL PRIMITIVES:
 * ------------------------------------------------------------------------- */

void dcopy (int_t n,
	    const double*  x, int_t incx, double*  y, int_t incy);
void icopy (int_t n, 
	    const int_t* x, int_t incx, int_t* y, int_t incy);
void scopy (int_t n,
	    const float*   x, int_t incx, float*   y, int_t incy);

void dfill (int_t n, double  alpha, double*  x, int_t incx);
void ifill (int_t n, int_t alpha, int_t* x, int_t incx);
void sfill (int_t n, float   alpha, float*   x, int_t incx);

void dznan (int_t n, double*  x, int_t incx);
void sznan (int_t n, float*   x, int_t incx);

void dneg (int_t n, double*  x, int_t incx);
void ineg (int_t n, int_t* x, int_t incx);
void sneg (int_t n, float*   x, int_t incx);

void dvneg (int_t n,
	    const double*  x, int_t incx, double*  y, int_t incy);
void ivneg (int_t n,
	    const int_t* x, int_t incx, int_t* y, int_t incy);
void svneg (int_t n,
	    const float*   x, int_t incx, float*   y, int_t incy);

void dvsgn (int_t n,
	    const double*  x, int_t incx, double*  y, int_t incy);
void ivsgn (int_t n,
	    const int_t* x, int_t incx, int_t* y, int_t incy);
void svsgn (int_t n,
	    const float*   x, int_t incx, float*   y, int_t incy);

void dsadd (int_t n, double  alpha, const double*  x, int_t incx, 
	                                    double*  y, int_t incy);
void isadd (int_t n, int_t alpha, const int_t* x, int_t incx,
	                                    int_t* y, int_t incy);
void ssadd (int_t n, float   alpha, const float*   x, int_t incx,
	                                    float*   y, int_t incy);

void dspow (const int_t n, const double alpha,
	    const double* x, int_t incx, double* y, int_t incy);
void sspow (const int_t n, const float  alpha,
	    const float*  x, int_t incx, float*  y, int_t incy);

void dvadd (int_t n, const double*  x, int_t incx,
	    const double*  y, int_t incy, double*  z, int_t incz);
void ivadd (int_t n, const int_t* x, int_t incx,
	    const int_t* y, int_t incy, int_t* z, int_t incz);
void svadd (int_t n, const float*   x, int_t incx,
	    const float*   y, int_t incy, float*   z, int_t incz);

void dssub (int_t n, double  alpha, const double*  x, int_t incx, 
	                                    double*  y, int_t incy);
void issub (int_t n, int_t alpha, const int_t* x, int_t incx,
	                                    int_t* y, int_t incy);
void sssub (int_t n, float   alpha, const float*   x, int_t incx,
	                                    float*   y, int_t incy);

void dvsub (int_t n, const double*  x, int_t incx,
	    const double*  y, int_t incy, double*  z, int_t incz);
void ivsub (int_t n, const int_t* x, int_t incx,
	    const int_t* y, int_t incy, int_t* z, int_t incz);
void svsub (int_t n, const float*   x, int_t incx,
	    const float*   y, int_t incy, float*   z, int_t incz);

void dsmul (int_t n, double  alpha, const double*  x, int_t incx,
	                                    double*  y, int_t incy);
void ismul (int_t n, int_t alpha, const int_t* x, int_t incx,
	                                    int_t* y, int_t incy);
void ssmul (int_t n, float   alpha, const float*   x, int_t incx,
	                                    float*   y, int_t incy);

void dvmul (int_t n, const double*  x, int_t incx,
	    const double*  y, int_t incy, double*  z, int_t incz);
void ivmul (int_t n, const int_t* x, int_t incx,
	    const int_t* y, int_t incy, int_t* z, int_t incz);
void svmul (int_t n, const float*   x, int_t incx,
	    const float*   y, int_t incy, float*   z, int_t incz);

void dsdiv (int_t n, double  alpha, const double*  x, int_t incx,
	                                    double*  y, int_t incy);
void isdiv (int_t n, int_t alpha, const int_t* x, int_t incx,
	                                    int_t* y, int_t incy);
void ssdiv (int_t n, float   alpha, const float*   x, int_t incx,
	                                    float*   y, int_t incy);

void dvrecp (int_t n,
	     const double* x, int_t incx, double* y, int_t incy);
void svrecp (int_t n,
	     const float*  x, int_t incx, float*  y, int_t incy);

void dvdiv (int_t n, const double*  x, int_t incx,
	    const double*  y, int_t incy, double*  z, int_t incz);
void svdiv (int_t n, const float*   x, int_t incx,
	    const float*   y, int_t incy, float*   z, int_t incz);

void dzero (int_t n, double*  x, int_t incx);
void izero (int_t n, int_t* x, int_t incx);
void szero (int_t n, float*   x, int_t incx);


/* ------------------------------------------------------------------------- *
 * OTHER MATHEMATICAL FUNCTIONS:
 * ------------------------------------------------------------------------- */

void dvabs (int_t n, const double*  x, int_t incx,
                             double*  y, int_t incy);
void ivabs (int_t n, const int_t* x, int_t incx,
	                     int_t* y, int_t incy);
void svabs (int_t n, const float*   x, int_t incx,
	                     float*   y, int_t incy);

void dvamax (int_t n, const double*  x, int_t incx,
	     const double*  y, int_t incy, double*  z, int_t incz);
void ivamax (int_t n, const int_t* x, int_t incx,
	     const int_t* y, int_t incy, int_t* z, int_t incz);
void svamax (int_t n, const float*   x, int_t incx,
	     const float*   y, int_t incy, float*   z, int_t incz);

void dvexp (int_t n, const double* x, int_t incx, double* y, int_t incy);
void svexp (int_t n, const float*  x, int_t incx, float*  y, int_t incy);

void dvlg10 (int_t n, const double* x, int_t incx,
	                      double* y, int_t incy);
void svlg10 (int_t n, const float*  x, int_t incx,
	                      float*  y, int_t incy);

void dvlog (int_t n, const double* x, int_t incx, double* y, int_t incy);
void svlog (int_t n, const float*  x, int_t incx, float*  y, int_t incy);

void dvatan (int_t n, const double* x, int_t incx,
	                      double* y, int_t incy);
void svatan (int_t n, const float*  x, int_t incx,
	                      float*  y, int_t incy);

void dvatn2 (int_t n, const double* x, int_t incx,
	     const double* y, int_t incy, double* z, int_t incz);
void svatn2 (int_t n, const float*  x, int_t incx,
	     const float*  y, int_t incy, float*  z, int_t incz);

void dvcos (int_t n, const double* x, int_t incx, double* y, int_t incy);
void svcos (int_t n, const float*  x, int_t incx, float*  y, int_t incy);

void dvsin (int_t n, const double* x, int_t incx, double* y, int_t incy);
void svsin (int_t n, const float*  x, int_t incx, float*  y, int_t incy);

void dvsqrt (int_t n, const double* x, int_t incx,
	                    double* y, int_t incy);
void svsqrt (int_t n, const float*  x, int_t incx,
	                    float*  y, int_t incy);

void dvtanh (int_t n, const double* x, int_t incx,
	                    double* y, int_t incy);
void svtanh (int_t n, const float*  x, int_t incx,
	                    float*  y, int_t incy);

void   raninit  (int_t flag);
double dranu    (void);
float  sranu    (void);
double drang    (void);
float  srang    (void);
double dnormal  (double mean, double sdev);
float  snormal  (float  mean, float  sdev);
void   dvrandom (int_t n, double* x, int_t incx);
void   svrandom (int_t n, float*  x, int_t incx);
void   dvgauss  (int_t n, double* x, int_t incx);
void   dsgauss  (int_t n, float*  x, int_t incx);
void   dvnormal (int_t n, double mean, double sdev, double* x, int_t incx);
void   svnormal (int_t n, float  mean, float  sdev, float*  x, int_t incx);

int_t ispow2   (int_t k);
int_t irpow2   (int_t k);
void zpreft (int_t k,                          zomplex *wtab, int_t sign);
void cpreft (int_t k,                          complex *wtab, int_t sign);
void zfft   (int_t n, zomplex *x,
	     int_t l, const zomplex *wtab, int_t dir);
void cfft   (int_t n, complex *x,
	     int_t l, const complex *wtab, int_t dir);
void dzfft  (int_t n, zomplex *x,
	     int_t l, const zomplex *wtab, int_t dir);
void scfft  (int_t n, complex *x,
	     int_t l, const complex *wtab, int_t dir);
void zpfft  (int_t n, zomplex *x,
	     int_t l, const zomplex *wtab, int_t dir);
void spfft  (int_t n, complex *x,
	     int_t l, const complex *wtab, int_t dir);

void dvhypot (int_t n, const double* x, int_t incx,
	      const double* y, int_t incy, double* z, int_t incz);
void svhypot (int_t n, const float*  x, int_t incx,
	      const float*  y, int_t incy, float*  z, int_t incz);
void dvmag   (int_t n, const double* w, int_t incw, 
	      const double* x, int_t incx, const double* y,
	      int_t incy, double* z, int_t incz);
void svmag   (int_t n, const float* w, int_t incw, 
	      const float*  x, int_t incx, const float*  y,
	      int_t incy, float*  z, int_t incz);

void dvpow (int_t n, const double* x, int_t incx,
	    const double* y, int_t incy, double* z, int_t incz);
void svpow (int_t n, const float*  x, int_t incx,
	    const float*  y, int_t incy, float*  z, int_t incz);


/* ------------------------------------------------------------------------- *
 * TRIAD OPERATIONS:
 * ------------------------------------------------------------------------- */

void dsvmvt (int_t n, double alpha, const double* x, int_t incx,
                                    const double* y, int_t incy,
	                                  double* z, int_t incz);
void ssvmvt (int_t n, float  alpha, const float*  x, int_t incx,
                                    const float*  y, int_t incy,
	                                  float*  z, int_t incz);

void dsvpvt (int_t n, double alpha, const double* x, int_t incx,
                                    const double* y, int_t incy,
	                                  double* z, int_t incz);
void ssvpvt (int_t n, float  alpha, const float*  x, int_t incx,
                                    const float*  y, int_t incy,
	                                  float*  z, int_t incz);

void dsvtsp (int_t n, double alpha, double beta,
	     const double* x, int_t incx, double* y, int_t incy);
void ssvtsp (int_t n, float  alpha, float  beta,
	     const float*  x, int_t incx, float*  y, int_t incy);

void dsvtvm (int_t n, double alpha, const double* x, int_t incx,
                                    const double* y, int_t incy,
	                                  double* z, int_t incz);
void ssvtvm (int_t n, float  alpha, const float*  x, int_t incx,
                                    const float*  y, int_t incy,
	                                  float*  z, int_t incz);

void dsvtvp (int_t n, double alpha, const double* x, int_t incx,
                                    const double* y, int_t incy,
	                                  double* z, int_t incz);
void ssvtvp (int_t n, float  alpha, const float*  x, int_t incx,
                                    const float*  y, int_t incy,
	                                  float*  z, int_t incz);

void dsvvmt (int_t n, double alpha, const double* x, int_t incx,
                                    const double* y, int_t incy,
	                                  double* z, int_t incz);
void ssvvmt (int_t n, float  alpha, const float*  x, int_t incx,
                                    const float*  y, int_t incy,
	                                  float*  z, int_t incz);

void dsvvpt (int_t n, double alpha, const double* x, int_t incx,
                                    const double* y, int_t incy,
	                                  double* z, int_t incz);
void ssvvpt (int_t n, float  alpha, const float*  x, int_t incx,
                                    const float*  y, int_t incy,
	                                  float*  z, int_t incz);

void dsvvtm (int_t n, double alpha, const double* x, int_t incx,
                                    const double* y, int_t incy,
	                                  double* z, int_t incz);
void ssvvtm (int_t n, float  alpha, const float*  x, int_t incx,
                                    const float*  y, int_t incy,
	                                  float*  z, int_t incz);

void dsvvtp (int_t n, double alpha, const double* x, int_t incx,
                                    const double* y, int_t incy,
	                                  double* z, int_t incz);
void ssvvtp (int_t n, float  alpha, const float*  x, int_t incx,
                                    const float*  y, int_t incy,
	                                  float*  z, int_t incz);

void dsvvtt (int_t n, double alpha, const double* x, int_t incx,
                                    const double* y, int_t incy,
	                                  double* z, int_t incz);
void ssvvtt (int_t n, float  alpha, const float*  x, int_t incx,
                                    const float*  y, int_t incy,
	                                  float*  z, int_t incz);

void dvvmvt (int_t n, const double* w, int_t incw, 
	              const double* x, int_t incx,
	              const double* y, int_t incy,
	                    double* z, int_t incz);
void svvmvt (int_t n, const float*  w, int_t incw,
                      const float*  x, int_t incx,
                      const float*  y, int_t incy,
	                    float*  z, int_t incz);

void dvvpvt (int_t n, const double* w, int_t incw,
	              const double* x, int_t incx,
                      const double* y, int_t incy,
	                    double* z, int_t incz);
void svvpvt (int_t n, const float*  w, int_t incw,
	              const float*  x, int_t incx,
                      const float*  y, int_t incy,
	                    float*  z, int_t incz);

void dvvpvt (int_t n, const double* w, int_t incw,
	              const double* x, int_t incx,
                      const double* y, int_t incy,
	                      double* z, int_t incz);
void svvpvt (int_t n, const float*  w, int_t incw,
	              const float*  x, int_t incx,
                      const float*  y, int_t incy,
	                    float*  z, int_t incz);

void dvvtvp (int_t n, const double* w, int_t incw,
	              const double* x, int_t incx,
                      const double* y, int_t incy,
	                    double* z, int_t incz);
void svvtvp (int_t n, const float*  w, int_t incw,
	              const float*  x, int_t incx,
                      const float*  y, int_t incy,
	                    float*  z, int_t incz);

void dvvtvm (int_t n, const double* w, int_t incw,
	              const double* x, int_t incx,
                      const double* y, int_t incy,
	                    double* z, int_t incz);
void svvtvm (int_t n, const float*  w, int_t incw,
	              const float*  x, int_t incx,
                      const float*  y, int_t incy,
	                    float*  z, int_t incz);

void dvvvtm (int_t n, const double* w, int_t incw,
	              const double* x, int_t incx,
                      const double* y, int_t incy,
	                    double* z, int_t incz);
void svvvtm (int_t n, const float*  w, int_t incw,
	              const float*  x, int_t incx,
                      const float*  y, int_t incy,
	                    float*  z, int_t incz);

void dvvvtt (int_t n, const double* w, int_t incw,
	              const double* x, int_t incx,
                      const double* y, int_t incy,
	                    double* z, int_t incz);
void svvvtt (int_t n, const float*  w, int_t incw,
	              const float*  x, int_t incx,
                      const float*  y, int_t incy,
	                    float*  z, int_t incz);


/* ------------------------------------------------------------------------- *
 * RELATIONAL PRIMITIVE OPERATIONS:
 * ------------------------------------------------------------------------- */

void iseq (int_t n, int_t alpha,
	   const int_t* x, int_t incx, int_t *y, int_t incy);
void dsge (int_t n, double  alpha,
	   const double*  x, int_t incx, int_t *y, int_t incy);
void isge (int_t n, int_t alpha,
	   const int_t* x, int_t incx, int_t *y, int_t incy);
void ssge (int_t n, float   alpha,
	   const float*   x, int_t incx, int_t *y, int_t incy);
void dsle (int_t n, double  alpha,
	   const double*  x, int_t incx, int_t *y, int_t incy);
void isle (int_t n, int_t alpha,
	   const int_t* x, int_t incx, int_t *y, int_t incy);
void ssle (int_t n, float   alpha,
	   const float*   x, int_t incx, int_t *y, int_t incy);
void dslt (int_t n, double  alpha,
	   const double*  x, int_t incx, int_t *y, int_t incy);
void islt (int_t n, int_t alpha,
	   const int_t* x, int_t incx, int_t *y, int_t incy);
void sslt (int_t n, float   alpha,
	   const float*   x, int_t incx, int_t *y, int_t incy);
void dsne (int_t n, double  alpha,
	   const double*  x, int_t incx, int_t *y, int_t incy);
void isne (int_t n, int_t alpha,
	   const int_t* x, int_t incx, int_t *y, int_t incy);
void ssne (int_t n, float   alpha,
	   const float*   x, int_t incx, int_t *y, int_t incy);


/* ------------------------------------------------------------------------- *
 * REDUCTION FUNCTIONS:
 * ------------------------------------------------------------------------- */

double dsum   (int_t n, const double*  x, int_t incx);
int_t  isum   (int_t n, const int_t* x, int_t incx);
float  ssum   (int_t n, const float*   x, int_t incx);
int_t  idmax  (int_t n, const double*  x, int_t incx);
int_t  iimax  (int_t n, const int_t* x, int_t incx);
int_t  ismax  (int_t n, const float*   x, int_t incx);
int_t  idmin  (int_t n, const double*  x, int_t incx);
int_t  iimin  (int_t n, const int_t* x, int_t incx);
int_t  ismin  (int_t n, const float*   x, int_t incx);
int_t  icount (int_t n, const int_t* x, int_t incx);
int_t  ifirst (int_t n, const int_t* x, int_t incx);
int_t  lany   (int_t n, const int_t* x, int_t incx);
int_t  lisame (int_t n, const int_t* x, int_t incx,
  	                const int_t* y, int_t incy);
int_t  ldsame (int_t n, const double*  x, int_t incx,
		        const double*  y, int_t incy);
int_t lssame (int_t n,  const float*   x, int_t incx,
		        const float*   y, int_t incy);

/* ------------------------------------------------------------------------- *
 * CONVERSION PRIMITIVES:
 * ------------------------------------------------------------------------- */

void vdble  (int_t n, const float*   x, int_t incx,
	                      double*  y, int_t incy);
void vsngl  (int_t n, const double*  x, int_t incx,
	                      float*   y, int_t incy);

void dvfloa (int_t n, const int_t* x, int_t incx,
	                      double*  y, int_t incy);
void svfloa (int_t n, const int_t* x, int_t incx,
	                      float*   y, int_t incy);

int_t iformat(void);
void format (char*);
void dbrev  (int_t n, const double*  x, int_t incx,
	                      double*  y, int_t incy);
void ibrev  (int_t n, const int_t* x, int_t incx,
	                      int_t* y, int_t incy);
void sbrev  (int_t n, const float*   x, int_t incx,
	                      float*   y, int_t incy);


/* ------------------------------------------------------------------------- *
 * MISCELLANEOUS FUNCTIONS:
 *
 * NB: xmxm & xmxv operations are replaced by macros that call BLAS xgemm,
 * xgemv routines (generally faster).
 * ------------------------------------------------------------------------- */

void dscatr (int_t n, double*  x, int_t *y, double*  z);
void iscatr (int_t n, int_t* x, int_t *y, int_t* z);
void sscatr (int_t n, float*   x, int_t *y, float*   z);

void dgathr (int_t n, double*  x, int_t *y, double*  z);
void igathr (int_t n, int_t* x, int_t *y, int_t* z);
void sgathr (int_t n, float*   x, int_t *y, float*   z);

void dramp (int_t n, double  alpha, double  beta, double*  x, int_t incx);
void iramp (int_t n, int_t alpha, int_t beta, int_t* x, int_t incx);
void sramp (int_t n, float   alpha, float   beta, float*   x, int_t incx);

void dclip (int_t n, const  double alpha,  const double   beta,
	    const double* x,  int_t incx,  double*  y, int_t incy);
void iclip (int_t n, const  int_t alpha, const int_t  beta,
	    const int_t* x, int_t incx,  int_t* y, int_t incy);
void sclip (int_t n, const  float alpha,   const float    beta,
	    const float* x,  int_t incx,   float*   y, int_t incy);

void dclipup (int_t n, const  double alpha,
	      const double* x,  int_t incx,  double*  y, int_t incy);
void iclipup (int_t n, const  int_t alpha,
	      const int_t* x, int_t incx,  int_t* y, int_t incy);
void sclipup (int_t n, const  float alpha,
	      const float* x,  int_t incx,   float*   y, int_t incy);

void dclipdn (int_t n, const  double alpha,
	      const double* x,  int_t incx,  double*  y, int_t incy);
void iclipdn (int_t n, const  int_t alpha,
	      const int_t* x, int_t incx,  int_t* y, int_t incy);
void sclipdn (int_t n, const  float alpha,
	      const float* x,  int_t incx,   float*   y, int_t incy);

void diclip (int_t n, const  double alpha,  const double   beta,
	     const double* x,  int_t incx,  double*  y, int_t incy);
void iiclip (int_t n, const  int_t alpha, const int_t  beta,
	     const int_t* x, int_t incx,  int_t* y, int_t incy);
void siclip (int_t n, const  float alpha,   const float    beta,
	     const float* x,  int_t incx,   float*   y, int_t incy);

void dcndst (int_t n, const  double*  x, int_t incx,
	     const int_t* y, int_t incy, double*  z, int_t incz);
void icndst (int_t n, const  int_t* x, int_t incx,
	     const int_t* y, int_t incy, int_t* z, int_t incz);
void scndst (int_t n, const  float*   x, int_t incx,
	     const int_t* y, int_t incy, float*   z, int_t incz);

void dmask (int_t n, const double*  w, int_t incw,
	               const double*  x, int_t incx,
	               const int_t* y, int_t incy,
	                     double*  z, int_t incz);
void imask (int_t n, const int_t* w, int_t incw,
	               const int_t* x, int_t incx,
	               const int_t* y, int_t incy,
	                     int_t* z, int_t incz);
void smask (int_t n, const float*   w, int_t incw,
	               const float*   x, int_t incx,
	               const int_t* y, int_t incy,
	                     float*   z, int_t incz);

void dvpoly (int_t n, const double* x, int_t incx, int_t m,
	                const double* c, int_t incc, 
		              double* y, int_t incy);
void svpoly (int_t n, const float*  x, int_t incx, int_t m,
	                const float*  c, int_t incc, 
		              float*  y, int_t incy);

double dpoly   (int_t n, double x, const double* xp, const double* yp);
float  spoly   (int_t n, float  x, const float*  xp, const float*  yp);
void   dpolint (const double* xa, const double* ya, int_t n,
	              double x,        double* y,  double* dy);
void   spolint (const float*  xa, const float*  ya, int_t n,
	              float   x,        float*  y, float*  dy);

void   dspline  (int_t n, double yp1, double ypn,
	         const double* x, const double* y,  double* y2);
void   sspline  (int_t n, float yp1, float ypn,
	         const float*  x, const float*  y,  float*  y2);
double dsplint  (int_t n, double x,
	         const  double* xa, const double* ya, const double* y2a);
float  ssplint  (int_t n, float  x,
	         const  float*  xa, const float*  ya, const float*  y2a);
double dsplquad (const double* xa, const double* ya,   const double* y2a,
		 const int_t n,      const double xmin,  const double xmax);
float  ssplquad (const float*  xa, const float*  ya,   const float*  y2a,
		 const int_t n,      const float  xmin,  const float  xmax);

void dmxm (double* A, int_t nra, double* B, int_t nca,
	   double* C, int_t ncb);
void smxm (float*  A, int_t nra, float*  B, int_t nca,
	   float*  C, int_t ncb);

void dmxv (double* A, int_t nra, double* B, int_t nca, double* C);
void smxv (float*  A, int_t nra, float*  B, int_t nca, float*  C);

void dmxva (double* A, int_t iac, int_t iar, double* B, int_t ib,
	    double* C, int_t ic,  int_t nra, int_t nca);
void smxva (float*  A, int_t iac, int_t iar, float*  B, int_t ib,
	    float*  C, int_t ic,  int_t nra, int_t nca);

void dvvtvvtp (int_t n, const double* v, int_t incv,
          	          const double* w, int_t incw,
	                  const double* x, int_t incx,
	                  const double* y, int_t incy,
	                        double* z, int_t incz);
void svvtvvtp (int_t n, const float*  v, int_t incv,
	                  const float*  w, int_t incw,
	                  const float*  x, int_t incx,
	                  const float*  y, int_t incy,
	                        float*  z, int_t incz);
void dvvtvvtm (int_t n, const double* v, int_t incv,
          	          const double* w, int_t incw,
	                  const double* x, int_t incx,
	                  const double* y, int_t incy,
	                        double* z, int_t incz);
void svvtvvtm (int_t n, const float*  v, int_t incv,
	                  const float*  w, int_t incw,
	                  const float*  x, int_t incx,
	                  const float*  y, int_t incy,
	                        float*  z, int_t incz);

void dsvvttvp (int_t n, const double  alpha,
		          const double* w, int_t incw,
		          const double* x, int_t incx,
		          const double* y, int_t incy,
		                double* z, int_t incz);
void ssvvttvp (int_t n, const float   alpha,
	                  const float*  w, int_t incw,
		          const float*  x, int_t incx,
		          const float*  y, int_t incy,
		                float*  z, int_t incz);  

/*****************************************************************************
 * Prototypes and macros for vendor-supplied FORTRAN libraries follow.
 *****************************************************************************/


/* ------------------------------------------------------------------------- *
 * BLAS level 1 selection.
 *
 * Note that the regular (FORTRAN) BLAS routines idamax & isamax have their
 * return values reduced by 1 to conform with usual C practice.
 * ------------------------------------------------------------------------- */

#if BLAS >= 1

double ddot_  (int_t* n, double* x, int_t* incx, double* y, int_t* incy);
double dasum_ (int_t* n, double* x, int_t* incx);
double dnrm2_ (int_t* n, double* x, int_t* incx);

void drotg_  (double* a, double* b, double* c, double* s);
void drot_   (int_t* n, double* x, int_t* incx, double* y, int_t* incy,
	      double* c, double* s);
void dswap_  (int_t* n, double* x, int_t* incx, double* y, int_t* incy);
void dscal_  (int_t* n, double* alpha, double* x, int_t* incx);
void daxpy_  (int_t* n, double* alpha, double* x, int_t* incx, 
		                         double* y, int_t* incy);
int_t idamax_ (int_t* n, double* x, int_t* incx);

float sdot_  (int_t* n, float*  x, int_t* incx, float*  y, int_t* incy);
float sasum_ (int_t* n, float*  x, int_t* incx);
float snrm2_ (int_t* n, float*  x, int_t* incx);

void srotg_  (float*  a, float*  b, float*  c, float*  s);
void srot_   (int_t* n, float *x, int_t* incx, float *y, int_t* incy,
	      float *c, float *s);
void sswap_  (int_t* n, float *x, int_t* incx, float *y, int_t* incy);
void sscal_  (int_t* n, float *alpha, float *x, int_t* incx);
void saxpy_  (int_t* n, float *alpha, float *x, int_t* incx, 
		                        float *y, int_t* incy);
int_t isamax_ (int_t* n, float *x, int_t* incx);


#define drotg(a, b, c, s) ( drotg_(a, b, c, s) )
#define dswap(n, x, incx, y, incy) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecIreg[2]=incy, \
   dswap_(_vecIreg, x ,_vecIreg+1, y, _vecIreg+2) )
#define dscal(n, alpha, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecDreg[0]=alpha, \
   dscal_(_vecIreg, _vecDreg, x, _vecIreg+1) )
#define daxpy(n, alpha, x, incx, y, incy) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecIreg[2]=incy, _vecDreg[0]=alpha, \
   daxpy_(_vecIreg, _vecDreg, x, _vecIreg+1, y, _vecIreg+2) )
#define ddot(n, x, incx, y, incy) \
  (_vecIreg[0]=n,_vecIreg[1]=incx,_vecIreg[2]=incy,\
   ddot_ (_vecIreg, x, _vecIreg+1, y, _vecIreg+2) )
#define dasum(n, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, dasum_(_vecIreg, x, _vecIreg+1) )
#define dnrm2(n, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, dnrm2_(_vecIreg, x, _vecIreg+1) )
#define drot(n, x, incx, y, incy, c, s) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecIreg[2]=incy, \
   _vecDreg[0]=c, _vecDreg[1]=s, \
   drot_(_vecIreg, x, _vecIreg+1, y, _vecIreg+2, _vecDreg, _vecDreg+1) )
#define idamax(n, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, idamax_(_vecIreg, x, _vecIreg+1)-1 )

#define srotg(a, b, c, s) ( srotg_(a, b, c, s) )
#define sswap(n, x, incx, y, incy) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecIreg[2]=incy, \
   sswap_(_vecIreg, x ,_vecIreg+1, y, _vecIreg+2) )
#define sscal(n, alpha, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecSreg[0]=alpha, \
   sscal_(_vecIreg, _vecSreg, x, _vecIreg+1) )
#define saxpy(n, alpha, x, incx, y, incy) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecIreg[2]=incy, _vecSreg[0]=alpha, \
   saxpy_(_vecIreg, _vecSreg, x, _vecIreg+1, y, _vecIreg+2) )
#define sdot(n, x, incx, y, incy) \
  (_vecIreg[0]=n,_vecIreg[1]=incx,_vecIreg[2]=incy,\
   sdot_ (_vecIreg, x, _vecIreg+1, y, _vecIreg+2) )
#define sasum(n, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, sasum_(_vecIreg, x, _vecIreg+1) )
#define snrm2(n, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, snrm2_(_vecIreg, x, _vecIreg+1) )
#define srot(n, x, incx, y, incy, c, s) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, _vecIreg[2]=incy, \
   _vecSreg[0]=c, _vecSreg[1]=s,                      \
   srot_(_vecIreg, x, _vecIreg+1, y, _vecIreg+2, _vecSreg, _vecSreg+1) )
#define isamax(n, x, incx) \
  (_vecIreg[0]=n, _vecIreg[1]=incx, isamax_(_vecIreg, x, _vecIreg+1)-1 )

#endif

/* ------------------------------------------------------------------------- *
 * BLAS level 2 selection.
 * ------------------------------------------------------------------------- */

#if BLAS >= 2

void dgemv_ (char *trans, int_t* m, int_t* n, double* alpha,
	     double* a, int_t* lda,
	     double* x, int_t* incx, double* beta, double* y, int_t* incy);
void sgemv_ (char *trans, int_t* m, int_t* n, float*  alpha,
	     float*  a, int_t* lda,
	     float*  x, int_t* incx, float*  beta, float*  y, int_t* incy);

void dger_ (int_t* m, int_t* n, double* alpha, double* x, int_t* incx,
	    double* y, int_t* incy, double* a, int_t* lda);
void sger_ (int_t* m, int_t* n, float*  alpha, float*  x, int_t* incx,
	    float*  y, int_t* incy, float*  a, int_t* lda);

void dspmv_(char *uplo, int_t* n, double* alpha, double* ap, double* x,
	    int_t* incx, double* beta, double* y, int_t* incy);
void sspmv_(char *uplo, int_t* n, float*  alpha, float*  ap, float*  x,
	    int_t* incx, float*  beta, float*  y, int_t* incy);

#define dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)           \
  (_vecCreg[0]=trans,_vecIreg[0]=m,_vecIreg[1]=n,_vecIreg[2]=lda, \
   _vecIreg[3]=incx,_vecIreg[4]=incy,_vecDreg[0]=alpha,           \
   _vecDreg[1]=beta,                                              \
   dgemv_(_vecCreg,_vecIreg,_vecIreg+1,_vecDreg,a,_vecIreg+2,x,   \
	  _vecIreg+3,_vecDreg+1,y,_vecIreg+4))
#define sgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)           \
  (_vecCreg[0]=trans,_vecIreg[0]=m,_vecIreg[1]=n,_vecIreg[2]=lda, \
   _vecIreg[3]=incx,_vecIreg[4]=incy,_vecSreg[0]=alpha,           \
   _vecSreg[1]=beta,                                              \
   sgemv_(_vecCreg,_vecIreg,_vecIreg+1,_vecSreg,a,_vecIreg+2,x,   \
	  _vecIreg+3,_vecSreg+1,y,_vecIreg+4))

#define dspmv(uplo,n,alpha,a,x,incx,beta,y,incy)        \
  (_vecCreg[0]=uplo,_vecIreg[0]=n,_vecIreg[1]=incx,     \
   _vecIreg[2]=incy,_vecDreg[0]=alpha,_vecDreg[1]=beta, \
   dspmv_(_vecCreg,_vecIreg,_vecDreg,a,x,               \
	  _vecIreg+1,_vecDreg+1,y,_vecIreg+2))
#define sspmv(uplo,n,alpha,a,x,incx,beta,y,incy)        \
  (_vecCreg[0]=uplo,_vecIreg[0]=n,_vecIreg[1]=incx,     \
   _vecIreg[2]=incy,_vecSreg[0]=alpha,_vecSreg[1]=beta, \
   sspmv_(_vecCreg,_vecIreg,_vecSreg,a,x,               \
	  _vecIreg+1,_vecSreg+1,y,_vecIreg+2))

#define dger(m,n,alpha,x,incx,y,incy,a,lda)           \
  (_vecIreg[0]=m,_vecIreg[1]=n,_vecDreg[0]=alpha,     \
   _vecIreg[2]=incx,_vecIreg[3]=incy,_vecIreg[4]=lda, \
   dger_(_vecIreg,_vecIreg+1,_vecDreg,x,_vecIreg+2,y,_vecIreg+3,a,_vecIreg+4))
#define sger(m,n,alpha,x,incx,y,incy,a,lda)           \
  (_vecIreg[0]=m,_vecIreg[1]=n,_vecSreg[0]=alpha,     \
   _vecIreg[2]=incx,_vecIreg[3]=incy,_vecIreg[4]=lda, \
   dger_(_vecIreg,_vecIreg+1,_vecSreg,x,_vecIreg+2,y,_vecIreg+3,a,_vecIreg+4))

#define dmxv(A,nra,x,nca,y)                                       \
  (_vecCreg[0]='T',_vecIreg[0]=nra,_vecIreg[1]=nca,_vecIreg[2]=1, \
   _vecDreg[0]=1.0,_vecDreg[1]=0.0,                               \
   dgemv_(_vecCreg,_vecIreg+1,_vecIreg,_vecDreg,A,_vecIreg+1,     \
	  x,_vecIreg+2,_vecDreg+1,y,_vecIreg+2))
#define smxv(A,nra,x,nca,y)\
  (_vecCreg[0]='T',_vecIreg[0]=nra,_vecIreg[1]=nca,_vecIreg[2]=1, \
   _vecSreg[0]=1.0,_vecSreg[1]=0.0,                               \
   sgemv_(_vecCreg,_vecIreg+1,_vecIreg,_vecSreg,A,_vecIreg+1,     \
	  x,_vecIreg+2,_vecSreg+1,y,_vecIreg+2))

#endif

/* ------------------------------------------------------------------------- *
 * BLAS level 3 selection.
 * ------------------------------------------------------------------------- */

#if BLAS == 3

void dgemm_ (char *ta, char *tb, int_t* m, int_t* n, int_t* k,
	     double* alpha,
	     double* a, int_t* lda,
	     double* b, int_t* ldb, double* beta,
	     double* c, int_t* ldc);
void sgemm_ (char *ta, char *tb, int_t* m, int_t* n, int_t* k,
	     float *alpha,
	     float *a, int_t* lda, float *b,
	     int_t* ldb, float *beta,
	     float *c, int_t* ldc);
 
#define dgemm(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)           \
  (_vecCreg[0]=ta,_vecCreg[1]=tb,_vecIreg[0]=m,_vecIreg[1]=n,     \
   _vecIreg[2]=k,_vecIreg[3]=lda,_vecIreg[4]=ldb,_vecIreg[5]=ldc, \
   _vecDreg[0]=alpha,_vecDreg[1]=beta,                            \
   dgemm_(_vecCreg,_vecCreg+1,_vecIreg,_vecIreg+1,_vecIreg+2,     \
	  _vecDreg,a,_vecIreg+3,b,_vecIreg+4,_vecDreg+1,c,        \
	  _vecIreg+5))
#define sgemm(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)           \
  (_vecCreg[0]=ta,_vecCreg[1]=tb,_vecIreg[0]=m,_vecIreg[1]=n,     \
   _vecIreg[2]=k,_vecIreg[3]=lda,_vecIreg[4]=ldb,_vecIreg[5]=ldc, \
   _vecSreg[0]=alpha,_vecSreg[1]=beta,                            \
   sgemm_(_vecCreg,_vecCreg+1,_vecIreg,_vecIreg+1,_vecIreg+2,     \
	  _vecSreg,a,_vecIreg+3,b,_vecIreg+4,_vecSreg+1,c,        \
	  _vecIreg+5))

#define dmxm(A,nra,B,nca,C,ncb)                                     \
  (_vecCreg[0]='N',_vecCreg[1]='N',_vecIreg[0]=nra,_vecIreg[1]=nca, \
   _vecIreg[2]=ncb,_vecDreg[0]=1.0,_vecDreg[1]=0.0,                 \
   dgemm_(_vecCreg,_vecCreg+1,_vecIreg+2,_vecIreg,_vecIreg+1,       \
	  _vecDreg,B,_vecIreg+2,A,_vecIreg+1,_vecDreg+1,C,_vecIreg+2))
#define smxm(A,nra,B,nca,C,ncb)                                     \
  (_vecCreg[0]='N',_vecCreg[1]='N',_vecIreg[0]=nra,_vecIreg[1]=nca, \
   _vecIreg[2]=ncb,_vecSreg[0]=1.0,_vecSreg[1]=0.0,                 \
   dgemm_(_vecCreg,_vecCreg+1,_vecIreg+2,_vecIreg,_vecIreg+1,       \
	  _vecSreg,B,_vecIreg+2,A,_vecIreg+1,_vecSreg+1,C,_vecIreg+2))


#endif

/* ------------------------------------------------------------------------- *
 * LAPACK selection.
 * ------------------------------------------------------------------------- */

#if LAPACK == 1

/* Factor a general matrix. */

void dgetrf_ (int_t* m, int_t* n, double* a, int_t* lda,
	      int_t* ipiv, int_t* info);
void sgetrf_ (int_t* m, int_t* n, float*  a, int_t* lda,
	      int_t* ipiv, int_t* info);

#define dgetrf(m,n,a,lda,ipiv,info)               \
  (_vecIreg[0]=m, _vecIreg[1]=n, _vecIreg[2]=lda, \
   dgetrf_(_vecIreg, _vecIreg+1, a, _vecIreg+2, ipiv, info))
#define sgetrf(m,n,a,lda,ipiv,info)               \
  (_vecIreg[0]=m, _vecIreg[1]=n, _vecIreg[2]=lda, \
   sgetrf_(_vecIreg, _vecIreg+1, a, _vecIreg+2, ipiv, info))

/* Invert a general matrix from factors. */

void dgetri_ (int_t* n,    double* a,     int_t* lda,   int_t* ipiv,  
	      double* work, int_t* lwork, int_t* info);
void sgetri_ (int_t* n,    float*  a,     int_t* lda,   int_t* ipiv,
	      float*  work, int_t* lwork, int_t* info);

#define dgetri(n,a,lda,ipiv,work,lwork,info)          \
  (_vecIreg[0]=n, _vecIreg[1]=lda, _vecIreg[2]=lwork, \
   dgetri_(_vecIreg, a, _vecIreg+1, ipiv, work, _vecIreg+2, info))
#define sgetri(n,a,lda,ipiv,work,lwork,info)          \
  (_vecIreg[0]=n, _vecIreg[1]=lda, _vecIreg[2]=lwork, \
   sgetri_(_vecIreg, a, _vecIreg+1, ipiv, work, _vecIreg+2, info))

/* Factor and solve band-symmetric matrix problem. */

void dpbtrf_ (char *uplo, int_t* n, int_t* kd,
	      double* ab, int_t* ldab, int_t* info);
void spbtrf_ (char *uplo, int_t* n, int_t* kd,
	      float*  ab, int_t* ldab, int_t* info);

#define dpbtrf(uplo,n,kd,ab,ldab,info)                                \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, _vecIreg[1]=kd, _vecIreg[2]=ldab, \
   dpbtrf_(_vecCreg, _vecIreg, _vecIreg+1, ab, _vecIreg+2, info))
#define spbtrf(uplo,n,kd,ab,ldab,info)                                \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, _vecIreg[1]=kd, _vecIreg[2]=ldab, \
   spbtrf_(_vecCreg, _vecIreg, _vecIreg+1, ab, _vecIreg+2, info))

void dpbtrs_ (char *uplo, int_t* n, int_t* kd, int_t* nrhs, double* ab,
	      int_t* ldab, double* b, int_t* ldb, int_t* info);
void spbtrs_ (char *uplo, int_t* n, int_t* kd, int_t* nrhs, float*  ab,
	      int_t* ldab, float*  b, int_t* ldb, int_t* info);

#define dpbtrs(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)                      \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, _vecIreg[1]=kd, _vecIreg[2]=nrhs,  \
   _vecIreg[3]=ldab, _vecIreg[4]=ldb,                                  \
   dpbtrs_(_vecCreg, _vecIreg, _vecIreg+1, _vecIreg+2, ab, _vecIreg+3, \
	   b, _vecIreg+4, info))
#define spbtrs(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)                      \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, _vecIreg[1]=kd, _vecIreg[2]=nrhs,  \
   _vecIreg[3]=ldab, _vecIreg[4]=ldb,                                  \
   spbtrs_(_vecCreg, _vecIreg, _vecIreg+1, _vecIreg+2, ab, _vecIreg+3, \
	   b, _vecIreg+4, info))

/* Factor and solve packed-symmetric matrix problem. */

void dpptrf_(char *uplo, int_t* n, double* ap, int_t* info);
void spptrf_(char *uplo, int_t* n, double* ap, int_t* info);

#define dpptrf(uplo,n,ap,info) \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, dpptrf_(_vecCreg, _vecIreg, ap, info))
#define spptrf(uplo,n,ap,info) \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, spptrf_(_vecCreg, _vecIreg, ap, info))

void dpptrs_(char *uplo, int_t* n, int_t* nrhs, double* ap,
	     double* b, int_t* ldb, int_t* info);
void spptrs_(char *uplo, int_t* n, int_t* nrhs, float*  ap,
	     float*  b, int_t* ldb, int_t* info);

#define dpptrs(uplo,n,nrhs,ap,b,ldb,info)                              \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, _vecIreg[1]=nrhs, _vecIreg[2]=ldb, \
   dpptrs_(_vecCreg, _vecIreg, _vecIreg+1, ap, b, _vecIreg+2, info))
#define spptrs(uplo,n,nrhs,ap,b,ldb,info)                              \
  (_vecCreg[0]=uplo, _vecIreg[0]=n, _vecIreg[1]=nrhs, _vecIreg[2]=ldb, \
   spptrs_(_vecCreg, _vecIreg, _vecIreg+1, ap, b, _vecIreg+2, info))

#endif

#endif

