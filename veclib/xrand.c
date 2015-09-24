/*****************************************************************************
 *                         RANDOM NUMBER GENERATION
 *
 * The following set of routines provides several types of random number
 * generation.  The only routines provided as part of VECLIB generate a
 * set of numbers distributed uniformly on (0,1), but here an extension to
 * normally-distributed RVs has been implemented.
 *
 * Source: Numerical Recipes 2nd edn.
 *
 * $Id: xrand.c,v 8.1 2015/04/20 11:14:20 hmb Exp $
 *****************************************************************************/

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

static double  UD     (double, double);
static double  GD     (double, double);
static double  ran1   (long *);	/* These are provided as alternatives */
static double  ran2   (long *);	/* but rand2 is the default. Change   */
static double  ran3   (long *);	/* _GENERATOR_ below to use another.  */
static double  gasdev (long *);
static long    iseed  = 0;

#define _GENERATOR_ ran2


void raninit (int_t flag)
/* ------------------------------------------------------------------------- *
 * Initialise random number generator.  Non-positive numbers
 * initialise the generator directly (with supplied value); if flag is
 * positive, the seed is generated from the wall-clock time (the value
 * of flag is irrelevant in this case).
 * ------------------------------------------------------------------------- */
{
  iseed = (flag > 0) ? -time (NULL) : flag;

  (void) _GENERATOR_ (&iseed);
}

 
double dranu (void)
/* ------------------------------------------------------------------------- *
 * Provide a single random number UD on (0, 1). 
 * ------------------------------------------------------------------------- */
{
  return UD (0.0, 1.0);
}


float sranu (void)
/* ------------------------------------------------------------------------- *
 * Provide a single random number UD on (0, 1). 
 * ------------------------------------------------------------------------- */
{
  return (float) UD (0.0, 1.0);
}

 
double drang (void)
/* ------------------------------------------------------------------------- *
 * Provide a single random number, Normal(0, 1). 
 * ------------------------------------------------------------------------- */
{
  return gasdev (&iseed);
}


float srang (void)
/* ------------------------------------------------------------------------- *
 * Provide a single random number, Normal(0, 1). 
 * ------------------------------------------------------------------------- */
{
  return (float) gasdev (&iseed);
}


double dnormal (double mean, double sdev)
/* ------------------------------------------------------------------------- *
 * Provide a single random number, Normal(mean, sdev). 
 * ------------------------------------------------------------------------- */
{
  return GD(mean, sdev);
}


float snormal (float mean, float sdev)
/* ------------------------------------------------------------------------- *
 * Provide a single random number, Normal(mean, sdev). 
 * ------------------------------------------------------------------------- */
{
  return (float) GD((double) mean, (double) sdev);
}


void dvrandom (int_t n, double* x, int_t incx)
/* ------------------------------------------------------------------------- *
 * Randomize vector x, UD on (0, 1).
 * ------------------------------------------------------------------------- */
{
  register int_t i;
  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = UD (0.0, 1.0);
}


void svrandom (int_t n, float* x, int_t incx)
/* ------------------------------------------------------------------------- *
 * Randomize vector x, UD on (0, 1).
 * ------------------------------------------------------------------------- */
{
  register int_t i;
  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = (float) UD (0.0, 1.0);
}


void dvgauss (int_t n, double* x, int_t incx)
/* ------------------------------------------------------------------------- *
 * Randomize vector x, Normal (0, 1).
 * ------------------------------------------------------------------------- */
{
  register int_t i;
  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = gasdev (&iseed);
}


void svaguss (int_t n, float* x, int_t incx)
/* ------------------------------------------------------------------------- *
 * Randomize vector x, Normal (0, 1).
 * ------------------------------------------------------------------------- */
{
  register int_t i;
  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = (float) gasdev (&iseed);
}


void dvnormal (int_t n, double mean, double sdev, double* x, int_t incx)
/* ------------------------------------------------------------------------- *
 * Randomize vector x, Normal(mean, sdev).
 * ------------------------------------------------------------------------- */
{
  register int_t i;
  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) x[i*incx] = GD (mean, sdev);
























}


void svnormal (int_t n, float mean, float sdev, float* x, int_t incx)
/* ------------------------------------------------------------------------- *
 * Randomize vector x, Normal(mean, sdev).
 * ------------------------------------------------------------------------- */
{
  register int_t i;
  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) x[i*incx] = (float) GD ((double) mean, (double) sdev);
}


/*****************************************************************************
 * Remaining routines only have file scope.
 *****************************************************************************/


static double UD (double low, double high)
/* ------------------------------------------------------------------------- *
 * Return a random number UD on (low, high).
 * ------------------------------------------------------------------------- */
{
  return _GENERATOR_ (&iseed) * (high - low) + low;
}


static double GD (double mean, double sdev)
/* ------------------------------------------------------------------------- *
 * Return normally distributed random deviate with specified mean and
 * standard deviation.
 * ------------------------------------------------------------------------- */
{
  return (mean + sdev * (gasdev (&iseed)));
}


static double gasdev (long* idum)
/* ------------------------------------------------------------------------- *
 * Returns a normally distributed deviate with zero mean and unit variance,
 * using _GENERATOR_(idum) as the source of uniform deviates. NR 2e.
 * ------------------------------------------------------------------------- */
{
  static int_t iset = 0;
  static double  gset;
  double         fac, rsq, v1, v2;

  if (iset == 0) {
    do {
      v1 = 2.0 * _GENERATOR_ (idum) - 1.0;
      v2 = 2.0 * _GENERATOR_ (idum) - 1.0;
      rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0 * log (rsq) / rsq);
    gset = v1 * fac;
    iset = 1;
    return v2 * fac;
  } else {
    iset = 0;
    return gset;
  }
}


#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

static double ran1 (long* idum)
/* ------------------------------------------------------------------------- *
 * Ran1 from NR.
 * Set idum to any negative value to initialize or re-initialize sequence.
 * ------------------------------------------------------------------------- */
{
  int    j;
  long   k;
  static long iy = 0;
  static long iv[NTAB];
  float  temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum = 1;
    else *idum = -(*idum);
    for (j = NTAB+7; j >= 0; j--) {
      k = (*idum)/IQ;
      *idum = IA*(*idum - k*IQ) - IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = (*idum)/IQ;
  *idum = IA*(*idum - k*IQ) - IR*k;
  if (*idum < 0) *idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if ((temp = AM*iy) > RNMX) return RNMX;
  else return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


#define IM1   2147483563
#define IM2   2147483399
#define AM    (1.0/IM1)
#define IMM1  (IM1-1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
#define NTAB  32
#define NDIV  (1+IMM1/NTAB)
#define EPS   1.2e-7
#define RNMX  (1.0-EPS)

static double ran2 (long* idum)
/* ------------------------------------------------------------------------- *
 * Ran2 from NR 2e.  Returns a uniform random deviate between 0.0 &
 * 1.0 (exclusive of endpoints).  Call with idum a negative int_t to
 * initialize; thereafter, do not alter idum between successive
 * deviates in a sequence.  RNMX should approximate the largest
 * floating value that is less than 1.
 * ------------------------------------------------------------------------- */
{
  int         j;
  long        k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float       temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum = 1;
    else              *idum = -(*idum);
    idum2 = (*idum);
    for (j=NTAB+7; j>=0; j--) {
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k*IQ1) - k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }

  k = (*idum) / IQ1;
  *idum = IA1*(*idum - k*IQ1) - k*IR1;
  if (*idum < 0) *idum += IM1;

  k = idum2 / IQ2;
  idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
  if (idum2 < 0) idum2 += IM2;

  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;

  if ((temp=AM*iy) > RNMX) return RNMX;
  else                     return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3 (long* idum)
/* ------------------------------------------------------------------------- *
 * Knuth's random number generator.
 * Set idum to any negative value to initialize or re-initialize sequence.
 * ------------------------------------------------------------------------- */
{
  static int  inext, inextp;
  static long ma[56];
  static int  iff = 0;
  long        mj, mk;
  int         i, ii, k;

  if (*idum < 0 || iff == 0) {
    iff = 1;
    mj = MSEED - (*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    for (i = 1; i <= 54; i++) {
      ii = (21*i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk < MZ) mk += MBIG;
      mj = ma[ii];
    }
    for (k = 1; k <= 4; k++) 
      for (i = 1; i <= 55; i++) {
	ma[i] -= ma[1 + (i + 30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext  = 0;
    inextp = 31;
    *idum  = 1;
  }
  if (++inext  == 56) inext  = 1;
  if (++inextp == 56) inextp = 1;
  mj = ma[inext] - ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext] = mj;
  return mj * FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
