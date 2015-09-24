/*****************************************************************************
 * fourier.c
 *
 * Copyright (c) 1999<-->$Date: 2015/04/20 11:14:14 $, Hugh Blackburn
 *
 * 1D Fourier transform routines for real data fields based on FFTPACK
 * or Temperton FFT routines, or vendor-supplied alternatives.
 * NB: different restrictions may apply to input args depending on
 * selected routine. 
 *
 * --
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
 *****************************************************************************/

/* $Id: fourier.c,v 8.1 2015/04/20 11:14:14 hmb Exp $ */

#include <stdio.h>
#include <stdlib.h>

#include <cfemdef.h>
#include <cveclib.h>
#include <cfemlib.h>

void dDFTr (double*     data,
	    const int_t tlen,
	    const int_t ntrn,
	    const int_t sign)
/* ------------------------------------------------------------------------- *
 * Carry out multiple 1D single--complex Fourier transforms of data.
 * Data is to be Fourier transformed in the direction normal to the most
 * rapid traverse through memory, with sucessive points in the transform
 * separated by ntrn.  Data has a total of tlen * ntrn real points.
 *
 * Input parameters:
 * ----------------
 * tlen: number of real data in each transform, product of prime numbers.
 * ntrn: number of transforms to perform (also skip in data).
 * sign: transform direction: +1 ==> r-->c, -1 ==> c-->r.
 *
 * Notes:
 * -----
 * (1) Data are scaled/normalized with 1/tlen when sign is +1, so that
 *     the zeroth Fourier mode contains the spatial average value.
 * (2) After forward (r-->c) transform, data are ordered so that within
 *     each transform, the zeroth mode datum comes first.  The zeroth
 *     mode is followed by the real datum from the maximum frequency mode,
 *     after which the real and imaginary parts for each mode alternate.
 * ------------------------------------------------------------------------- */
{
  const char      routine[] = "dDFTr";
  char            err[STR_MAX];
  const int_t     ntot = tlen * ntrn;
  register int_t  i;
  int_t           dum, ip, iq, ir, ipqr2, *ifax;
  register double *work, *Wtab, *ptr;

  if (tlen < 2 || !ntrn) return;

#if defined(_SX)  /* -- Use NEC FFT routines. */

  ifax = ivector (0, 63);
  work = dvector (0, ntot + tlen - 1);
  Wtab = work + ntot;

  rftfax (tlen, ifax, Wtab);

  if (ifax[0] == -99)
    message (routine, "tlen needs prime factors 2, 3, 5", ERROR);

  if (sign == FORWARD) {
    rfft  (data, work, Wtab, ifax, tlen, ntrn, 1.0 / (double) tlen);
    dcopy ((tlen - 2) * ntrn, data +              ntrn, 1, work,            1);
    dcopy (             ntrn, data + (tlen - 1) * ntrn, 1, data + ntrn,     1);
    dcopy ((tlen - 2) * ntrn, work,                     1, data + 2 * ntrn, 1);
  } else {
    dcopy ((tlen - 2) * ntrn, data + 2 * ntrn, 1, work,                     1);
    dcopy (             ntrn, data + ntrn,     1, data + (tlen - 1) * ntrn, 1);
    dcopy ((tlen - 2) * ntrn, work,            1, data +              ntrn, 1);
    rfft  (data, work, Wtab, ifax, tlen, ntrn, -1.0);
  }

  freeIvector (ifax, 0);
  freeDvector (work, 0);

#elif defined(DEBUG_FFT) /* -- Unvectorized FFTPACK routines. */

  work = dvector (0, 3 * tlen - 1);
  ifax = ivector (0, 14);
  Wtab = work + 2 * tlen;
  ptr  = data;

  drffti (tlen, Wtab, ifax);

  switch (sign) {

  case FORWARD:
    for (i = 0; i < ntrn; i++, ptr++) {
      dcopy  (tlen, ptr, ntrn, work, 1);
      drfftf (tlen, work, work + tlen, Wtab, ifax);
      dcopy  (tlen - 2, work + 1, 1, ptr + 2 * ntrn, ntrn);
      ptr[0]    = work[0];
      ptr[ntrn] = work[tlen - 1];
    }
    dscal (ntot, 1.0 / tlen, data, 1);
    break;

  case INVERSE:
    for (i = 0; i < ntrn; i++, ptr++) {
      work[tlen - 1] = ptr[ntrn];
      work[0]        = ptr[0];
      dcopy  (tlen - 2, ptr + 2 * ntrn, ntrn, work + 1, 1);
      drfftb (tlen, work, work + tlen, Wtab, ifax);
      dcopy  (tlen, work, 1, ptr, ntrn);
    }
    break;

  default:
    message (routine, "illegal direction flag", ERROR);
    break;
  }
  
  freeDvector (work, 0);
  freeIvector (ifax, 0);

#else /* -- Temperton FFT routine is default. */

  dum = tlen;
  prf235 (&dum, &ip, &iq, &ir, &ipqr2);
  
  if (!dum) {
    sprintf (err, "transform length (%1d) needs prime factors 2, 3, 5", tlen);
    message (routine, err, ERROR);
  }
  if (ntrn & 1) {
    sprintf (err, "number of transforms (%1d) must be even", ntrn);
    message (routine, err, ERROR);
  }

  work = dvector (0, ntot + ipqr2 - 1);
  Wtab = work + ntot;

  dsetpf (Wtab, tlen, ip, iq, ir);
  dmpfft (data, work, ntrn, tlen, ip, iq, ir, Wtab, sign);
  if (sign == FORWARD) dscal (ntot, 1.0 / tlen, data, 1);

  freeDvector (work, 0);

#endif
}
