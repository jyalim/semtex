//////////////////////////////////////////////////////////////////////////////
// helmholtz.C:  routines to solve elliptic problems in one variable.
//
// Copyright (c) 1994<-->$Date: 2015/04/20 11:14:14 $, Hugh Blackburn
//
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
// 02110-1301 USA
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: helmholtz.cpp,v 8.1 2015/04/20 11:14:14 hmb Exp $";

#include <sem.h>


void Helmholtz (Domain*   D,
		AuxField* F)
// ---------------------------------------------------------------------------
// Solve Helmholtz's equation
//                                  2
//               div grad u - LAMBDA  u = f(x,y)
//
// subject to BCs.
// ---------------------------------------------------------------------------
{
  const real_t lambda2 = Femlib::value ("LAMBDA2");
  const real_t beta    = Femlib::value ("BETA");
  const int_t  nmodes  = Geometry::nModeProc();
  const int_t  base    = Geometry::baseMode();
  const int_t  nz      = Geometry::nZProc();
  SolverKind   method  = (Femlib::ivalue("ITERATIVE")) ? JACPCG : DIRECT;

  ModalMatrixSys* M = new ModalMatrixSys
    (lambda2, beta, base, nmodes, D -> elmt, D -> b[0], method);

  D -> u[0] -> solve (F, M);

  D -> step++;
}
