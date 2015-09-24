///////////////////////////////////////////////////////////////////////////////
// bsys.C: BoundarySys class functions.
//
// Copyright (c) 1999 <--> $Date: 2015/04/20 11:14:17 $, Hugh Blackburn
//
// The information to be returned by class functions are the global
// numbering scheme and vector of Boundary*'s for a given Field and
// Fourier mode number.  There is one BoundarySys for each Field, but
// a possible modal dependence for the appropriate BCs (in fact, only
// for 3D cylindrical coordinate systems in which the axis appears).
//
// Use of cylindrical coordinates is flagged by the Geometry class
// variable.  In the case where the number of space dimensions is also
// 3, the number of boundary frames and numbering systems is set to 3,
// for the 0th, 1st and 2nd (and higher) modes, irrespective of the
// number of Fourier modes actually used.  See bcmgr.C.
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

static char RCS[] = "$Id: bsys.cpp,v 8.1 2015/04/20 11:14:17 hmb Exp $";

#include <sem.h>


BoundarySys::BoundarySys (BCmgr*                  bcmgr,
			  const vector<Element*>& elmt ,
			  const char              name ) :
// ---------------------------------------------------------------------------
// Construct internal storage for boundary systems for all modes.
// ---------------------------------------------------------------------------
  _field_name (name),
  _nbound     (bcmgr -> nBCedges()),
  _mixed      (false)
{
  const int_t                 verbose = Femlib::ivalue ("VERBOSE");
  const int_t                 np = Geometry::nP();
  vector<BCtriple*>::iterator edge;
  BCtriple*                   BCT;
  const Condition*            C;
  const char*                 S;
  char                        buf[StrMax], group;
  int_t                       i, j, k, offset;

  _number = new NumberSys* [3];
  for (i = 0; i < 3; i++) _number[i] = bcmgr -> getNumberSys (_field_name, i);

  _boundary = new vector<Boundary*> [3];

  if (!_nbound) { for (i = 0; i < 3; i++) _boundary[i].resize (0); return; }

  // -- Construct vectors of Boundary pointers using bcmgr.

  for (i = 0; i < 3; i++) _boundary[i].resize (_nbound);
 
  // -- Mode 0 boundaries, and default settings for other modes.

  edge = bcmgr -> getBCedges().begin();  
  for (offset = 0, i = 0; i < _nbound; i++, edge++, offset += np) {
    BCT   = *edge;
    group = BCT -> group;
    j     = BCT -> elmt;
    k     = BCT -> side;
    S     = bcmgr -> groupInfo    (group);
    C     = bcmgr -> getCondition (group, _field_name, 0);
    
    C -> describe (buf);
    if (strstr (buf, "mixed") || strstr (buf, "convective")) _mixed = true;
    
    _boundary[0][i] =
    _boundary[1][i] =
    _boundary[2][i] = new Boundary (i, S, C, elmt[j], k);
  }

  if (!(Geometry::system() == Geometry::Cylindrical && Geometry::nDim() == 3))
    return;

  // -- Mode 1 boundaries, adjusted on axis.

  edge = bcmgr -> getBCedges().begin();  
  for (offset = 0, i = 0; i < _nbound; i++, edge++, offset += np) {
    BCT   = *edge;
    group = BCT -> group;
    j     = BCT -> elmt;
    k     = BCT -> side;
    S     = bcmgr -> groupInfo    (group);
    C     = bcmgr -> getCondition (group, _field_name, 1);

    if (strstr (S, "axis"))
      _boundary[1][i] = new Boundary (i, "axis", C, elmt[j], k);
  }

  // -- Mode 2 boundaries,adjusted on axis.
  
  edge = bcmgr -> getBCedges().begin();
  for (offset = 0, i = 0; i < _nbound; i++, edge++, offset += np) {
    BCT   = *edge;
    group = BCT -> group;
    j     = BCT -> elmt;
    k     = BCT -> side;
    S     = bcmgr -> groupInfo    (group);
    C     = bcmgr -> getCondition (group, _field_name, 2);

    if (strstr (bcmgr -> groupInfo (group), "axis"))
      _boundary[2][i] = new Boundary (i, "axis", C, elmt[j], k);
  }
}


const vector<Boundary*>& BoundarySys::BCs (const int_t mode) const
// ---------------------------------------------------------------------------
// Return appropriate vector of Boundary*'s for field name, according
// to Fourier mode.  Mode number is the actual number, counting from
// zero, not modulo number of modes on this process.
// ---------------------------------------------------------------------------
{
  return
    _boundary [clamp (mode,static_cast<int_t>(0),static_cast<int_t>(2))];
}


const NumberSys* BoundarySys::Nsys (const int_t mode) const
// ---------------------------------------------------------------------------
// Return appropriate NumberSystem* for field name, according to
// Fourier mode.
// ---------------------------------------------------------------------------
{
  return
    _number [clamp (mode,static_cast<int_t>(0),static_cast<int_t>(2))];
}


const real_t* BoundarySys::Imass (const int_t mode) const
// ---------------------------------------------------------------------------
// Return appropriate inverse mass matrix for field name, according to
// Fourier mode.
// ---------------------------------------------------------------------------
{
  return 
    _number [clamp(mode,static_cast<int_t>(0),static_cast<int_t>(2))]->imass();
}
