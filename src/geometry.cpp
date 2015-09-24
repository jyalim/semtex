///////////////////////////////////////////////////////////////////////////////
// geometry.C: define geometrical properties for 2D quad X Fourier spaces.
//
// Copyright (c) 1994 <--> $Date: 2015/04/20 11:14:18 $, Hugh Blackburn
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
// --
//
// Most routines are inlined in header file geometry.h
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: geometry.cpp,v 8.1 2015/04/20 11:14:18 hmb Exp $";

#include <cstdio>
#include <iostream>

#include <cfemdef.h>
#include <utility.h>
#include <geometry.h>
#include <femlib.h>

int_t Geometry::_pid   = 0;
int_t Geometry::_nproc = 0;
int_t Geometry::_ndim  = 0;
int_t Geometry::_np    = 0;
int_t Geometry::_nz    = 0;
int_t Geometry::_nzp   = 0;
int_t Geometry::_nel   = 0;
int_t Geometry::_psize = 0;
Geometry::CoordSys Geometry::_csys  = Geometry::Cartesian;


void Geometry::set (const int_t    NP,
		    const int_t    NZ,
		    const int_t    NE,
		    const CoordSys CS)
// ---------------------------------------------------------------------------
// Load values of static internal variables.
//
// The number of processors is restricted: it must either be 1 or an
// even number, and it must be less than or equal to the number of
// planes / 2.  Furthermore, the number of planes on a processor must
// be even, unless NZ == 1.
//
// NB: the value of psize is the value of nPlane, but rounded up if
// necessary to be an even number and also an int_t multiple of the
// number of processors.  The even number restriction is to simplify
// the handling of Fourier transforms, which can be based on a
// real--complex transform on some platforms.  The restriction to be
// an int_t multiple of the number of processors is to simplify the
// structure of memory exchanges required for Fourier transforms.
// ---------------------------------------------------------------------------
{
  static char routine[] = "Geometry::set", err[StrMax];

  _pid   = Femlib::ivalue ("I_PROC");
  _nproc = Femlib::ivalue ("N_PROC");

  _np   = NP; _nz = NZ; _nel = NE; _csys = CS;
  _nzp  = _nz / _nproc;
  _ndim = (_nz > 1) ? 3 : 2;

  if (_nz > 1 && _nz & 1) {	// -- 3D problems must have NZ even.
    sprintf (err, "N_Z must be even (%1d)", _nz);
    message (routine, err, ERROR);
  }

  if (_nproc > 1) {		// -- Concurrent execution restrictions.
    if (_nproc & 1) {
      sprintf (err, "No. of processors must be even (%1d)",
	       _nproc);
      message (routine, err, ERROR);
    }

    if (_nproc << 1 > _nz) {
      sprintf (err, "No. of processors (%1d) can at most be half N_Z (%1d)",
	       _nproc, _nz);
      message (routine, err, ERROR);
    }

    if (_nz % (2 * _nproc)) {
      sprintf (err, "No. of planes (%1d) per processor (%1d) must be even",
	       _nz, _nproc);
      message (routine, err, ERROR);
    }

    _psize  = nPlane();
    _psize += 2 * _nproc - nPlane() % (2 * _nproc);

  } else {

    _psize = nPlane() + (nPlane() % 2);
  }
}
