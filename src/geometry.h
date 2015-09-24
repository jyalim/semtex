#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "cfemdef.h"

class Geometry
// ===========================================================================
// Details of the logical/storage (as opposed to geometric)
// representation used for scalar fields.  Static functions make
// information globally accessible.
//
// Copyright (c) 1994 <--> $Date: 2015/04/20 11:14:18 $, Hugh Blackburn
//
// In all cases, 2D quad elements are employed, with a possible
// extension by Fourier expansions in the third dimension.  While the
// representation implied by this class is not necessarily conforming,
// the same order of interpolation is used in each element.
// Equal-order interpolation is used in each direction on faces of
// quads.
//
// With the introduction of concurrent execution, the concept of
// geometry has been extended to include the processor ID, number of
// processors, number of data planes per processor, etc.
//
// $Id: geometry.h,v 8.1 2015/04/20 11:14:18 hmb Exp $
// ===========================================================================
{
public:
  enum CoordSys { Cartesian, Cylindrical };

  static void set (const int_t, const int_t, const int_t, const CoordSys);

  static CoordSys system    () { return _csys;                 }  
  static bool     cylindrical () { return _csys == Geometry::Cylindrical; }

  static int_t  nP        () { return _np;                   }
  static int_t  nZ        () { return _nz;                   }
  static int_t  nElmt     () { return _nel;                  }
  static int_t  nTotElmt  () { return _np * _np;             }
  static int_t  nExtElmt  () { return 4 * (_np - 1);         }
  static int_t  nIntElmt  () { return (_np - 2) * (_np - 2); }
  static int_t  nMode     () { return (_nz + 1) >> 1;        }
  static int_t  nDim      () { return _ndim;                 }
  static int_t  nPlane    () { return _nel * nTotElmt();     }
  static int_t  nBnode    () { return _nel * nExtElmt();     }
  static int_t  nInode    () { return _nel * nIntElmt();     }
  static int_t  nTot      () { return _nz  * nPlane();       }
  static int_t  planeSize () { return _psize;                }
  static int_t  nTotal    () { return _nz * _psize;          }

  static int_t  nProc     () { return _nproc;                }
  static int_t  procID    () { return _pid;                  }
  static int_t  nZProc    () { return _nzp;                  }
  static int_t  nZ32      () { return (_nproc > 1) ? _nzp : (3 * _nz) >> 1; }
  static int_t  nTotProc  () { return _nzp * _psize;         }
  static int_t  nModeProc () { return nMode() / _nproc;      }
  static int_t  baseMode  () { return _pid * nModeProc();    }
  static int_t  basePlane () { return _pid * _nzp;           }
  static int_t  nBlock    () { return _psize / _nproc;       }

private:
  static int_t    _nproc ;	// Number of processors.
  static int_t    _pid   ;	// ID for this processor, starting at 0.
  static int_t    _ndim  ;	// Number of space dimensions.
  static int_t    _np    ;	// Number of points along element edge.
  static int_t    _nz    ;	// Number of planes (total).
  static int_t    _nzp   ;	// Number of planes per processor.
  static int_t    _nel   ;	// Number of elements.
  static int_t    _psize ;	// nPlane rounded up to suit restrictions.
  static CoordSys _csys  ;	// Coordinate system (Cartesian/cylindrical).

};
#endif
