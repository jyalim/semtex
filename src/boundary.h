#ifndef BOUNDARY_H
#define BOUNDARY_H

class Boundary : public Edge
// ===========================================================================
// Physical field element-wise boundary class.
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
  Boundary (const int_t id, const char* group, const Condition* bcond,
	    const Element* elmt, const int_t side):
  Edge (group, elmt, side), _id (id), _bcond (bcond) { }

  int_t ID        () const { return _id; }
  void  print     () const;

  void  evaluate  (const Field*,const int_t,const int_t,const bool,
		   real_t*)                                              const;

  // -- Impose essential BCs:
  void  set       (const real_t*,const int_t*,real_t*)                   const;

  // -- Apply natural BCs:
  void  sum       (const real_t*,const int_t*,real_t*,real_t*)           const;

  // -- Apply mixed BCs:
  void  augmentSC (const int_t,const int_t,const int_t*,real_t*,real_t*) const;
  void  augmentOp (const int_t*,const real_t*,real_t*)                   const;
  void  augmentDg (const int_t*,real_t*)                                 const;

private:
  int_t            _id   ;	// Ident number.
  const Condition* _bcond;	// Boundary condition.
};

#endif
