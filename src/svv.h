#ifndef SVV_H
#define SVV_H

#include <cfemdef.h>

// ===========================================================================
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

namespace SVV {

  const real_t* coeffs    (const int_t);
  const real_t* coeffs_z  (const int_t);

  void          operators (const int_t, const real_t**, const real_t**);
}

#endif
