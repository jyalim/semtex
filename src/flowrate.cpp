//////////////////////////////////////////////////////////////////////////////
// flowrate.C: routines to estimate and set flow rate using Stokes
// Green's functions.
//
// Copyright (c) 2003 <--> $Date: 2015/04/20 11:14:18 $, Hugh Blackburn
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
// The flow rate per unit width is measured across a cut through the
// domain, defined in session's "CUT" section, which contains a set of
// element sides over which the flow rate per unit width is measured.
//
// SYNTAX:
//
// <CUT NUMBER=4>
//   1   4  2
//   2   8  2
//   3  12  2
//   4  16  2
// </CUT>
//
// On each line, the first number is a reference tag, and the second
// and third numbers give the element and side number of element
// boundaries along the cut. Note that the orientation of the unit
// outward normals should be continuous (not reversed) from one
// element to the next, and the cut should span the domain, roughly
// normal to the dominant flow direction but certainly not parallel to
// it.
//
// Also used is the file session.grn, which contains field data for
// the Stokes Green's function for the flow.
//
// The TOKENS section should define the FLOWRATE, the desired flowrate
// per unit width.
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: flowrate.cpp,v 8.1 2015/04/20 11:14:18 hmb Exp $";

#include <sem.h>


Flowrate::Flowrate (Domain* D,
		    FEML*   F) :
// ---------------------------------------------------------------------------
// This sets up the data extraction surface and stores the Stokes
// Green's function from file.
// ---------------------------------------------------------------------------
  _src  (D),
  _refQ (Femlib::value ("FLOWRATE"))
{
  // -- We shouldn't be sent here if there is no CUT section; check anyway.

  if (!(F -> seek ("CUT"))) return;

  const char routine[] = "Flowrate::Flowrate";
  char       buf[StrMax];
  int_t      i, id, elmt, side, Nedge;

  // -- Get the information for data extraction surface.

  _curve.resize (Nedge = F -> attribute ("CUT", "NUMBER"));

  for (i = 0; i < Nedge; i++) {
    while (F -> stream().peek() == '#') // -- Skip comments.
      F -> stream().ignore (StrMax, '\n');
    F -> stream() >> id >> elmt >> side;
    
    if (elmt > Geometry::nElmt() || elmt < 1) {
      sprintf (buf, "egde %1d: element no. %1d out of range [1, %1d]",
	       id, elmt, Geometry::nElmt());
      message (routine, buf, ERROR);
    }
    
    if (side > 4 || side < 1) {
      sprintf (buf, "egde %1d: side no. %1d out of range [1, 4]",
	       id, side);
      message (routine, buf, ERROR);
    }
    
    _curve[i] = new Edge ("flowrate", _src -> elmt[--elmt], --side);
  }

  // -- Now pull in Stokes Green's function from session.grn

  char     file[StrMax];
  ifstream grnfunc (strcat (strcpy (file, _src -> name), ".grn"));
  real*    alloc = new real [2 * Geometry::planeSize()];
  Header   header;

  if (!grnfunc) {
    sprintf (buf, "can't open Stokes Green's function file %s", file);
    message (routine, buf, ERROR);
  }

  // -- Check the header indicates it conforms with session file.

  grnfunc >> header;

  if (header.nr != Geometry::nP() || header.nel != Geometry::nElmt())
    message (routine, "session and Stokes forcing flow do not conform", ERROR);
  if (header.nz != 1)
    message (routine, "Stokes forcing flow must be 2D (N_Z=1)", ERROR);
  if (strcmp (header.flds, "uvp"))
    message (routine, "Stokes forcing flow must have only u, v, p", ERROR);

  _green.resize (2);
  _green[0] = new AuxField (alloc,                      1, _src -> elmt, 'u');
  _green[1] = new AuxField (alloc + Geometry::nPlane(), 1, _src -> elmt, 'v');

  grnfunc >> *_green[0];
  grnfunc >> *_green[1];

  if (header.swab()) {
    _green[0] -> reverse();
    _green[1] -> reverse();
  }

  Femlib::synchronize();	// -- Wait until root process has done this.
}


real Flowrate::getQ () const
// ---------------------------------------------------------------------------
// Estimate the volume flow rate per unit width across the cut.  This
// is done for the zeroth Fourier mode, so should be called only on
// the root process in parallel operation.
// ---------------------------------------------------------------------------
{
  vector<Edge*>::const_iterator c;
  vector<real> work (2 * Geometry::nP());
  real Q = 0.0;
  
  for (c = _curve.begin(); c != _curve.end(); c++)
    Q += (*c)->vectorFlux ("flowrate", _src->udat[0], _src->udat[1], &work[0]);
  
  return Q;
}


