/*****************************************************************************
 * filter.c
 *
 * Copyright (c) 1999<-->$Date: 2015/04/20 11:14:14 $, Hugh Blackburn
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
 * Routines for computing spectral filters.
 *
 * $Id: filter.c,v 8.1 2015/04/20 11:14:14 hmb Exp $
 *****************************************************************************/

#include <math.h>
#include <stdio.h>

#include <cfemdef.h>
#include <cveclib.h>


void bvdFilter (const int_t  N     ,
		const real_t order ,
		const real_t lag   ,
		const real_t attn  ,
		real_t*      filter)
/* ------------------------------------------------------------------------- *
 * Load filter with the Boyd--Vandeven (i.e. erfc) filter [0, N] of given
 * order (p) and lag (s).
 *
 * NB: N should be one less than the number of coefficients to
 * which the filter will be applied.
 * 
 * Input parameter attn gives the attenuation at high
 * wavenumbers. 0<=attn<=1, with attn = 1 giving complete attenuation.
 * 
 * Input parameter lag is the proportion (0--1) of the length 0--N at
 * which filter rolloff starts.
 *
 * Reference:
 * J.G. Levin & M. Iskandarani & D.B. Haidvogel (1997), A spectral
 * filtering procedure for eddy-resolving simulations with a spectral
 * element ocean model, JCP V137, 130--154.
 * ------------------------------------------------------------------------- */
{
  int_t        i;
  real_t       arg, theta, chi, omega;
  const real_t shift = N * lag, EPS = EPSSP;

  for (i = 0; i <= N; i++) {
    if (i < shift) 
      filter[i] = 1.0;

    else {
      theta = (i - shift) / (N - shift);

      if      (fabs (theta - 0.0) < EPS)
	filter[i] = 1.0;

      else if (fabs (theta - 1.0) < EPS)
	filter[i] = 1.0 - attn;

      else if (fabs (theta - 0.5) < EPS)
	filter[i] = 1.0 - attn * 0.5;

      else {
	omega     = fabs (theta) - 0.5;
	arg       = 1.0 - 4.0 * SQR (omega);
	chi       = sqrt (-log (arg) / (4.0 * SQR (omega)));
	filter[i] = (1.0 - attn) + attn* 0.5*erfc (2.0*sqrt(order)*chi*omega);
      }
    }
  }
}
