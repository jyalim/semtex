#ifndef BSYS_H
#define BSYS_H

class BoundarySys
// ===========================================================================
// This class automates the retrieval of the boundary condition
// applicators (Boundary objects), global numbering schemes
// (NumberSys) and inverse mass matrix for a given Field and Fourier
// mode.
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
// ===========================================================================
{
public:
  BoundarySys  (BCmgr*, const vector<Element*>&, const char);
  ~BoundarySys () { };

  char                     field () const { return _field_name; }
  int_t                    nSurf () const { return _nbound; }
  bool                     mixBC () const { return _mixed; }
  const vector<Boundary*>& BCs   (const int_t) const;
  const NumberSys*         Nsys  (const int_t) const;
  const real_t*            Imass (const int_t) const;

private:
  char               _field_name;
  int_t              _nbound    ;  // Number of element edges with BCs.
  bool               _mixed     ;  // Flags presence of mixed BC type.
  vector<Boundary*>* _boundary  ;  // Boundary*'s  for modes 0, 1, 2.
  NumberSys**        _number    ;  // NumberSys*'s for modes 0, 1, 2.
};

#endif
