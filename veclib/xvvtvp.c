/*****************************************************************************
 * xvvtvp:  z[i] = (w[i] * x[i]) + y[i]. 
 *
 * $Id: xvvtvp.c,v 8.1 2015/04/20 11:14:21 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

void dvvtvp (int_t n,
	     const double* w, int_t incw,
	     const double* x, int_t incx,
	     const double* y, int_t incy,
	           double* z, int_t incz)
{
#if defined(__DECC)
  /* -- DEC OSF 1, with Loop unrolling from KAPC. */

    register int i;
    register unsigned int _Kii1;
    register unsigned int _Kii2;
    register unsigned int _Kii3;
    register unsigned int _Kii4;
    int    _Kii5;
    double _Kdd1;
    double _Kdd2;
    double _Kdd3;
    double _Kdd4;
    
    if (incw == 1 && incx == 1 && incy == 1 && incz == 1) {
      if (!(((void  *)(z + n - 1) < (void  *)(w + (int )(0)) ||
	     (void  *)(z + (int )(0)) > (void  *)(w + n - 1)) &&
	    ((void  *)(z + n - 1) < (void  *)(x + (int )(0)) ||
	     (void  *)(z + (int )(0)) > (void  *)(x + n - 1)) && 
	    ((void  *)(z + n - 1) < (void  *)(y + (int )(0)) ||
	     (void  *)(z + (int )(0)) > (void  *)(y + n - 1)))) {
	for ( _Kii5 = (int )(0); _Kii5<=n - 4; _Kii5+=4 ) {
	  z[_Kii5] = y[_Kii5] + w[_Kii5] * x[_Kii5];
	  z[_Kii5+1] = y[_Kii5+1] + w[_Kii5+1] * x[_Kii5+1];
	  z[_Kii5+2] = y[_Kii5+2] + w[_Kii5+2] * x[_Kii5+2];
	  z[_Kii5+3] = y[_Kii5+3] + w[_Kii5+3] * x[_Kii5+3];
	}
	for ( ; _Kii5<n; _Kii5++ ) {
	  z[_Kii5] = y[_Kii5] + w[_Kii5] * x[_Kii5];
	}
      } else {
	for ( _Kii5 = (int )(0); _Kii5<=n - 4; _Kii5+=4 ) {
	  _Kdd2 = w[_Kii5] * x[_Kii5];
	  _Kdd3 = w[_Kii5+1] * x[_Kii5+1];
	  _Kdd4 = w[_Kii5+2] * x[_Kii5+2];
	  _Kdd1 = w[_Kii5+3] * x[_Kii5+3];
	  z[_Kii5] = y[_Kii5] + _Kdd2;
	  z[_Kii5+1] = y[_Kii5+1] + _Kdd3;
	  z[_Kii5+2] = y[_Kii5+2] + _Kdd4;
	  z[_Kii5+3] = y[_Kii5+3] + _Kdd1;
	}
	for ( ; _Kii5<n; _Kii5++ ) {
	  _Kdd1 = w[_Kii5] * x[_Kii5];
	  z[_Kii5] = y[_Kii5] + _Kdd1;
	}
      } 
    } else {
      if (incw < 0) {
	_Kii1 = (unsigned int )((1 - n) * incw);
      } else {
	_Kii1 = 0;
      } 
      w +=  _Kii1;
      if (incx < 0) {
	_Kii2 = (unsigned int )((1 - n) * incx);
      } else {
	_Kii2 = 0;
      } 
      x +=  _Kii2;
      if (incy < 0) {
	_Kii3 = (unsigned int )((1 - n) * incy);
      } else {
	_Kii3 = 0;
      } 
      y +=  _Kii3;
      if (incz < 0) {
	_Kii4 = (unsigned int )((1 - n) * incz);
      } else {
	_Kii4 = 0;
      } 
      z +=  _Kii4;
      for ( _Kii5 = (int )(0); _Kii5<=n - 4; _Kii5+=4 ) {
	z[_Kii5*incz]     =
	  y[_Kii5*incy]     + w[_Kii5*incw] * x[_Kii5*incx];
	z[(_Kii5+1)*incz] =
	  y[(_Kii5+1)*incy] + w[(_Kii5+1)*incw] * x[(_Kii5+1)*incx];
	z[(_Kii5+2)*incz] =
	  y[(_Kii5+2)*incy] + w[(_Kii5+2)*incw] * x[(_Kii5+2)*incx];
	z[(_Kii5+3)*incz] =
	  y[(_Kii5+3)*incy] + w[(_Kii5+3)*incw] * x[(_Kii5+3)*incx];
      }
      for ( ; _Kii5<n; _Kii5++ ) {
	z[_Kii5*incz]     =
	  y[_Kii5*incy]      + w[_Kii5*incw] * x[_Kii5*incx];
      }
    } 

#else
  register int_t i;

  if (incw == 1 && incx == 1 && incy == 1 && incz == 1) 
   for (i = 0; i < n; i++) z[i] = (w[i] * x[i]) + y[i]; 

  else {

    w += (incw<0) ? (-n+1)*incw : 0;
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    z += (incz<0) ? (-n+1)*incz : 0;

    for (i = 0; i < n; i++) z[i*incz] = (w[i*incw] * x[i*incx]) + y[i*incy];
  }

#endif
}


void svvtvp (int_t n,
	     const float* w, int_t incw,
	     const float* x, int_t incx,
	     const float* y, int_t incy,
	           float* z, int_t incz)
{
#if defined(__DECC)
    register int i;
    register unsigned int _Kii1;
    register unsigned int _Kii2;
    register unsigned int _Kii3;
    register unsigned int _Kii4;
    int _Kii5;
    float _Krr1;
    float _Krr2;
    float _Krr3;
    float _Krr4;
    
    if (incw == 1 && incx == 1 && incy == 1 && incz == 1) {
        if (!(((void  *)(z + n - 1) < (void  *)(w + (int )(0)) ||
	       (void  *)(z + (int )(0)) > (void  *)(w + n - 1)) &&
	      ((void  *)(z + n - 1) < (void  *)(x + (int )(0)) ||
	       (void  *)(z + (int )(0)) > (void  *)(x + n - 1)) &&
	      ((void  *)(z + n - 1) < (void  *)(y + (int )(0)) ||
	       (void  *)(z + (int )(0)) > (void  *)(y + n - 1)))) {
            for ( _Kii5 = (int )(0); _Kii5<=n - 4; _Kii5+=4 ) {
                z[_Kii5] = y[_Kii5] + w[_Kii5] * x[_Kii5];
                z[_Kii5+1] = y[_Kii5+1] + w[_Kii5+1] * x[_Kii5+1];
                z[_Kii5+2] = y[_Kii5+2] + w[_Kii5+2] * x[_Kii5+2];
                z[_Kii5+3] = y[_Kii5+3] + w[_Kii5+3] * x[_Kii5+3];
            }
            for ( ; _Kii5<n; _Kii5++ ) {
                z[_Kii5] = y[_Kii5] + w[_Kii5] * x[_Kii5];
            }
        } else {
            for ( _Kii5 = (int )(0); _Kii5<=n - 4; _Kii5+=4 ) {
                _Krr2 = w[_Kii5] * x[_Kii5];
                _Krr3 = w[_Kii5+1] * x[_Kii5+1];
                _Krr4 = w[_Kii5+2] * x[_Kii5+2];
                _Krr1 = w[_Kii5+3] * x[_Kii5+3];
                z[_Kii5] = y[_Kii5] + _Krr2;
                z[_Kii5+1] = y[_Kii5+1] + _Krr3;
                z[_Kii5+2] = y[_Kii5+2] + _Krr4;
                z[_Kii5+3] = y[_Kii5+3] + _Krr1;
            }
            for ( ; _Kii5<n; _Kii5++ ) {
                _Krr1 = w[_Kii5] * x[_Kii5];
                z[_Kii5] = y[_Kii5] + _Krr1;
            }
        } 
    } else {
        if (incw < 0) {
            _Kii1 = (unsigned int )((1 - n) * incw);
        } else {
            _Kii1 = 0;
        } 
        w +=  _Kii1;
        if (incx < 0) {
            _Kii2 = (unsigned int )((1 - n) * incx);
        } else {
            _Kii2 = 0;
        } 
        x +=  _Kii2;
        if (incy < 0) {
            _Kii3 = (unsigned int )((1 - n) * incy);
        } else {
            _Kii3 = 0;
        } 
        y +=  _Kii3;
        if (incz < 0) {
            _Kii4 = (unsigned int )((1 - n) * incz);
        } else {
            _Kii4 = 0;
        } 
        z +=  _Kii4;
        for ( _Kii5 = (int )(0); _Kii5<=n - 4; _Kii5+=4 ) {
            z[_Kii5*incz]     =
	      y[_Kii5*incy]     + w[_Kii5*incw] * x[_Kii5*incx];
            z[(_Kii5+1)*incz] =
	      y[(_Kii5+1)*incy] + w[(_Kii5+1)*incw] * x[(_Kii5+1)*incx];
            z[(_Kii5+2)*incz] =
	      y[(_Kii5+2)*incy] + w[(_Kii5+2)*incw] * x[(_Kii5+2)*incx];
            z[(_Kii5+3)*incz] =
	      y[(_Kii5+3)*incy] + w[(_Kii5+3)*incw] * x[(_Kii5+3)*incx];
        }
        for ( ; _Kii5<n; _Kii5++ ) {
            z[_Kii5*incz]     =
	      y[_Kii5*incy]      + w[_Kii5*incw] * x[_Kii5*incx];
        }
    }
 
#else
  register int_t i;

  if (incw == 1 && incx == 1 && incy == 1 && incz == 1) 
   for (i = 0; i < n; i++) z[i] = w[i] * x[i] + y[i]; 
  
  else {
    w += (incw<0) ? (-n+1)*incw : 0;
    x += (incx<0) ? (-n+1)*incx : 0;
    y += (incy<0) ? (-n+1)*incy : 0;
    z += (incz<0) ? (-n+1)*incz : 0;

    for (i = 0; i < n; i++) z[i*incz] = w[i*incw] * x[i*incx] + y[i*incy];
  }

#endif
}
