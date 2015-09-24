///////////////////////////////////////////////////////////////////////////////
// integration.C:  supply coefficients for discrete time integration schemes.
//
// Copyright (c) 1994 <--> $Date: 2015/04/20 11:14:18 $, Hugh Blackburn
//
// Maximum time order supported is 3 (4 for implicit Adams--Moulton methods).
// Coefficients for all schemes can be found in Gear's book, "Numerical
// Initial Value Problems in Ordinary Differential Equations", 1971.
//
// --
// This file is part of Semtex.
// 
// Semtex is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// Semtex is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Semtex (see the file COPYING); if not, write to the Free
// Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
// 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: integration.cpp,v 8.1 2015/04/20 11:14:18 hmb Exp $";

#include <sem.h>


const int_t Integration::OrderMax = 4;


void Integration::AdamsBashforth  (const int_t n    ,
				   real_t*     coeff)
// ---------------------------------------------------------------------------
// Adams--Bashforth (predictor) coefficients of order n.  Gear, Table 7.3.
// ---------------------------------------------------------------------------
{
  char routine[] = "Integration::AdamsBashforth";

  switch (n) {
  case 1:
    coeff[0] =  1.0;
    break;
  case 2:
    coeff[0] =  1.5;
    coeff[1] = -0.5;
    break;
  case 3:
    coeff[0] =  23.0 / 12.0;
    coeff[1] = -16.0 / 12.0;
    coeff[2] =   5.0 / 12.0;
    break;
  default:
    message (routine, "requested order out of range", ERROR);
    break;
  }
}


void Integration::AdamsMoulton (const int_t n    ,
				real_t*     coeff)
// ---------------------------------------------------------------------------
// Adams--Moulton (corrector) coefficients of order n.  Gear, Table 7.5.
// ---------------------------------------------------------------------------
{
  char routine[] = "Integration::AdamsMoulton";

  switch (n) {
  case 1:
    coeff[0] =  1.0;
    break;
  case 2:
    coeff[0] =  0.5;
    coeff[1] =  0.5;
    break;
  case 3:
    coeff[0] =  5.0 / 12.0;
    coeff[1] =  8.0 / 12.0;
    coeff[2] = -1.0 / 12.0;
    break;
  case 4:
    coeff[0] =  9.0 / 24.0;
    coeff[1] = 19.0 / 24.0;
    coeff[2] = -5.0 / 24.0;
    coeff[3] =  1.0 / 24.0;
    break;
  default:
    message (routine, "requested order out of range", ERROR);
    break;
  }
}


void Integration::StifflyStable (const int_t n    ,
				 real_t*     coeff)
// ---------------------------------------------------------------------------
// "Stiffly-stable" backwards differentiation coefficients of order n.
// NB: vector coeff must be of length n + 1.  First coefficient in each
// case applies to the new time level.  Gear, Table 11.1.
// 
// NB: Karniadakis, Israeli & Orszag JCP 97 (1991), KIO91, also use
// these coefficients but their nomenclature differs. Their gamma_0 is
// the same as coeff[0], and their alpha vector holds the *negatives*
// of the remaining coefficients.
// ---------------------------------------------------------------------------
{
  char routine[] = "Integration::StifflyStable";

  switch (n) {
  case 1:
    coeff[0] =  1.0;
    coeff[1] = -1.0;
    break;
  case 2:
    coeff[0] =  1.5;
    coeff[1] = -2.0;
    coeff[2] =  0.5;
    break;
  case 3:
    coeff[0] =  11.0 / 6.0;
    coeff[1] = -3.0;
    coeff[2] =  1.5;
    coeff[3] = -1.0 / 3.0;
    break;
  default:
    message (routine, "requested order out of range", ERROR);
    break;
  }
}


void Integration::Extrapolation (const int_t n    ,
				 real_t*     coeff)
// ---------------------------------------------------------------------------
// Coefficients of order n for explicit extrapolation to end of timestep.
//
// These are the coefficients beta in KIO91.
// ---------------------------------------------------------------------------
{
  char routine[] = "Integration::Extrapolation";

  switch (n) {
  case 1:
    coeff[0] =  1.0;
    break;
  case 2:
    coeff[0] =  2.0;
    coeff[1] = -1.0;
    break;
  case 3:
    coeff[0] =  3.0;
    coeff[1] = -3.0;
    coeff[2] =  1.0;
    break;
  default:
    message (routine, "requested order out of range", ERROR);
    break;
  }
}
