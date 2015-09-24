///////////////////////////////////////////////////////////////////////////////
// history.C
//
// Copyright (c) 1994 <--> $Date: 2015/04/20 11:14:18 $, Hugh Blackburn
//
// Routines to provide history point information at x, y, z locations.
//
// History points are nominated in the session file, e.g.
// <HISTORY NUMBER=2>
// #    tag   x     y    z
//      1     0.0   0.0  0.0
//      2     0.3  -0.1 -2
// </HISTORY>
//
// The z locational information is ignored for 2D.  For 3D the periodic
// assumption is used in the z-direction in case the z location does not
// fall within the first spatial period.
//
// Output to a file called session.his.
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

static char RCS[] = "$Id: history.cpp,v 8.1 2015/04/20 11:14:18 hmb Exp $";

#include <sem.h>


const Element* HistoryPoint::locate (const real_t      x   ,
				     const real_t      y   ,
				     vector<Element*>& Esys,
				     real_t&           r   ,
				     real_t&           s   )
// ---------------------------------------------------------------------------
// Static class member function which tries to locate x, y, point within
// an element E, and return its location in r, s coordinates within E.
// ---------------------------------------------------------------------------
{
  register int_t i;
  const int_t    NEL   = Esys.size();
  const bool     guess = true;
  const Element* E;
  vector<real_t> work (max (2*Geometry::nTotElmt(), 5*Geometry::nP() + 6));

  for (E = 0, i = 0; E == 0 && i < NEL; i++)
    if (Esys[i] -> locate (x, y, r=0.0, s=0.0, &work[0], guess)) E = Esys[i];

  return E;
}


void HistoryPoint::extract (vector<AuxField*>& u  ,
			    real_t*            tgt) const
// ---------------------------------------------------------------------------
// Load tgt with information extracted from each AuxField in u, tgt is
// assumed to have sufficient storage to suit.
// ---------------------------------------------------------------------------
{
  register int_t i;
  const int_t    N = u.size();

  for (i = 0; i < N; i++) tgt[i] = u[i] -> probe (_E, _r, _s, _z);
}

