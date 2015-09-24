///////////////////////////////////////////////////////////////////////////////
// matrix.C: routines that generate solvers for Helmholtz problems.
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
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: matrix.cpp,v 8.1 2015/04/20 11:14:18 hmb Exp $";

#include <sem.h>


static vector<MatrixSys*> MS;


ModalMatrixSys::ModalMatrixSys (const real_t            lambda2 ,
				const real_t            beta    ,
				const int_t             baseMode,
				const int_t             numModes,
				const vector<Element*>& Elmt    ,
				const BoundarySys*      Bsys    ,
				const SolverKind        method  )
// ---------------------------------------------------------------------------
// Generate or retrieve from internal database MS the vector of
// MatrixSys's which will be used to solve all the Fourier-mode
// discrete Helmholtz problems for the associated scalar Fields
// (called out by names).
//
// Input variables:
//   lambda2 : Helmholtz constant for the problem,	
//   beta    : Fourier length scale = TWOPI / Lz,
//   nmodes  : number of Fourier modes which will be solved,
//   Elmt    : vector of Element*'s used to make local Helmholtz matrices,
//   Bsys    : boundary system for this field.
//   method  : specify the kind of solver we want (Cholesky, PCG ...).
// ---------------------------------------------------------------------------
{
  const char name = Bsys -> field();
  int_t      mode;
  bool       found;

  MatrixSys* M;
  vector<MatrixSys*>::iterator m;

  _fields = new char [strlen (Bsys -> Nsys (0) -> fields()) + 1];
  strcpy (_fields, Bsys -> Nsys (0) -> fields());
  _Msys.resize (numModes);

  if (method == DIRECT) {
    ROOTONLY cout << "-- Installing matrices for field '" << name << "' [";
    cout.flush();
    Femlib::synchronize();
  }

  for (mode = baseMode; mode < baseMode + numModes; mode++) {
    const int_t      modeIndex = mode * Femlib::ivalue ("BETA");
    const NumberSys* N         = Bsys -> Nsys(modeIndex);
    const real_t     betak2    = sqr (Field::modeConstant(name,mode,beta));
    const int_t      localMode = mode - baseMode;

    // Multiply Helmholtz constant with SVV-specific weight:
    //   betak2_svv = betak2 * (1 + eps_N/nu * Q) for modes k > SVV_MZ
    //   and lambda2 > 0 (i.e. only for the velocity components)
 
    const real_t* S = SVV::coeffs_z (numModes);
    const real_t  betak2_svv = (lambda2>EPSDP)?(betak2*S[localMode]) : betak2; 

    for (found = false, m = MS.begin(); !found && m != MS.end(); m++) {
      M     = *m;
      found = M -> match (lambda2, betak2_svv, N, method);
    }
    if (found) {
      _Msys[localMode] = M;
      if (method == DIRECT) { cout << '.'; cout.flush(); }
    } else {
      _Msys[localMode] =
	new MatrixSys (lambda2, betak2_svv, modeIndex, Elmt, Bsys, method);
      MS.insert (MS.end(), _Msys[localMode]);
      if (method == DIRECT) { cout << '*'; cout.flush(); }
    }
  }

  if (method == DIRECT) {
    Femlib::synchronize();
    ROOTONLY cout << "]" << endl;
    cout.flush();
  }
}


ModalMatrixSys::~ModalMatrixSys ()
// ---------------------------------------------------------------------------
// Destructor hands off calls to MatrixSys::~MatrixSys.  Note there
// can be side effects here, owing to the multiple use of MatrixSys*'s
// in different ModalMatrixSys's.  If multiple related
// ModalMatrixSys's exist, it is advisable to delete and recreate all
// of them before attempting reuse.
// ---------------------------------------------------------------------------
{
  int_t N = _Msys.size();
  vector<MatrixSys*>::iterator p;
  while (N--) {
    for (p = MS.begin(); p != MS.end(); p++)
      if (*p == _Msys[N]) { delete (_Msys[N]); MS.erase(p); break; }
  }
}


// ===========================================================================
// The routines above deal with vectors of MatrixSys class elements.
// The routines below deal with these elements themselves.
// ===========================================================================


MatrixSys::MatrixSys (const real_t            lambda2,
		      const real_t            betak2 ,
		      const int_t             mode   ,
		      const vector<Element*>& elmt   ,
		      const BoundarySys*      bsys   ,
		      const SolverKind        method ) :
// ---------------------------------------------------------------------------
// Initialize and factorise matrices in this system.
//
// For method == DIRECT:
//   Matrices are assembled using LAPACK-compatible ordering systems.
//   Global Helmholtz matrix uses symmetric-banded format; elemental
//   Helmholtz matrices (hii & hbi) use column-major formats.
// For method == JACPCG:
//   Build and invert diagonal preconditioner matrix.
//
// The Fourier-modal dependence of BCs and numbering is only really
// relevant for cylindrical systems (and at the axis); for Cartesian
// systems there is only a single (though replicated) set.
// --------------------------------------------------------------------------
// NB: these get evaluated in the order they appear in the class
// definition!:
  _HelmholtzConstant (lambda2),
  _FourierConstant   (betak2 ),
  _BC                (bsys -> BCs  (mode)),
  _NS                (bsys -> Nsys (mode)),
  _nel               (Geometry::nElmt()),
  _nglobal           (_NS -> nGlobal()),
  _singular          ((_HelmholtzConstant + _FourierConstant) < EPSSP &&
		      !_NS -> fmask() && !bsys -> mixBC()),
  _nsolve            ((_singular) ? _NS -> nSolve() - 1 : _NS -> nSolve()),
  _method            (method),
  _nband             (_NS -> nBand()),
  _npack             (_nband * _nsolve),
  _H                 (0),
  _hbi               (0),
  _hii               (0),
  _bipack            (0),
  _iipack            (0),
  _npts              (_nglobal + Geometry::nInode()),
  _PC                (0)
{
  const char     routine[] = "MatrixSys::MatrixSys";
  const int_t    verbose   = Femlib::ivalue ("VERBOSE");
  const int_t    np        = Geometry::nP();
  const int_t    next      = Geometry::nExtElmt();
  const int_t    nint      = Geometry::nIntElmt();
  const int_t    npnp      = Geometry::nTotElmt();
  const int_t*   bmap;
  register int_t i, j, k, m, n;

  switch (_method) {

  case DIRECT: {
    vector<real_t> work     (sqr (next) + sqr (np) + sqr (npnp));
    vector<int_t>  pivotmap (nint);
    real_t*        hbb  = &work[0];
    real_t*        rmat = hbb  + sqr (next);
    real_t*        rwrk = rmat + sqr (np);
    int_t*         ipiv = &pivotmap[0];
    int_t          info;

    _hbi    = new real_t*[static_cast<size_t>(_nel)];
    _hii    = new real_t*[static_cast<size_t>(_nel)];
    _bipack = new int_t  [static_cast<size_t>(_nel)];
    _iipack = new int_t  [static_cast<size_t>(_nel)];

    if (_nsolve) {
      _H = new real_t [static_cast<size_t>(_npack)];
      Veclib::zero (_npack, _H, 1);

      if (verbose > 1)
	cout << endl
	     << "Helmholtz constant (lambda2): " << setw(10) << lambda2
	     << ", Fourier constant (betak2): "  << setw(10) << betak2;
      if (verbose)
	cout << endl << "System matrix: " << _nsolve << "x" << _nband
	     << "\t(" << _npack << " words)";
    }

    // -- Loop over elements, creating & posting elemental Helmholtz matrices.

    for (bmap = _NS -> btog(), j = 0; j < _nel; j++, bmap += next) {
      _bipack[j] = next * nint;
      _iipack[j] = nint * nint;

      if (nint) {
	_hbi[j] = new real_t [static_cast<size_t>(_bipack[j])];
	_hii[j] = new real_t [static_cast<size_t>(_iipack[j])];
	Veclib::zero (_bipack[j], _hbi[j], 1);
	Veclib::zero (_iipack[j], _hii[j], 1);
      } else
	_hbi[j] = _hii[j] = 0;
      

      elmt[j]->HelmholtzSC(lambda2,betak2,hbb,_hbi[j],_hii[j],rmat,rwrk,ipiv);

      for (i = 0; i < next; i++)
	if ((m = bmap[i]) < _nsolve)
	  for (k = 0; k < next; k++)
	    if ((n = bmap[k]) < _nsolve && n >= m)
	      _H[Lapack::band_addr (m, n, _nband)] +=
		hbb[Veclib::row_major (i, k, next)];

      Family::adopt (_bipack[j], _hbi + j);
      Family::adopt (_iipack[j], _hii + j);

    }
    if (_nsolve) {
      // -- Loop over BCs and add diagonal contribution from mixed BCs.

      if (bsys -> mixBC()) {
	const int_t  nbound = bsys -> nSurf();
	const int_t* bmap   = _NS  -> btog();
	for (i = 0; i < nbound; i++)
	  _BC[i] -> augmentSC (_nband, _nsolve, bmap, rwrk, _H);
      }

      // -- Cholesky factor global banded-symmetric Helmholtz matrix.
    
      Lapack::pbtrf ("U",_nsolve,_nband-1,_H,_nband,info);

      if (info) message (routine, "failed to factor Helmholtz matrix", ERROR);

      Family::adopt (_npack, &_H);

      if (verbose) {
	real_t cond;
	pivotmap.resize (_nsolve);  ipiv = &pivotmap[0];
	work.resize (3 * _nsolve);  rwrk = &work[0];

	Lapack::pbcon ("U",_nsolve,_nband-1,_H,_nband,1.0,cond,rwrk,ipiv,info);
	cout << ", (inverse) condition number: " << cond << endl;
      }
    }
  } break;

  case JACPCG: {
    const int_t    nbound = _BC.size();   
    real_t*        PCi;
    vector<real_t> work (2 * npnp + np);
    real_t         *ed = &work[0], *ewrk = &work[0] + npnp;

    if (verbose > 1)
      cout << "PCG "
	   << "Helmholtz const (lambda2): " << setw(10) << lambda2
	   << ", Fourier const (betak2): "  << setw(10) << betak2 << endl;


    _PC = new real_t [static_cast<size_t>(_npts)];

    Veclib::zero (_npts, _PC, 1);

    PCi  = _PC + _NS -> nGlobal();
    bmap = _NS -> btog();

    // -- Mixed BC contributions.

    if (bsys -> mixBC())
      for (i = 0; i < nbound; i++)
	_BC[i] -> augmentDg (bmap, _PC);

    // -- Element contributions.

    for (i = 0; i < _nel; i++, bmap += next, PCi += nint) {
      elmt[i] -> HelmholtzDiag (lambda2, betak2, ed, ewrk);
      Veclib::scatr_sum (next, ed,  bmap,    _PC);
      Veclib::copy      (nint, ed + next, 1, PCi, 1);
    }

#if 1
    Veclib::vrecp (_npts, _PC, 1, _PC, 1);
#else  // -- Turn off preconditioner for testing.
    Veclib::fill  (_npts, 1.0, _PC, 1);
#endif

    Family::adopt (_npts, &_PC);

  } break;

  default:
    message (routine, "no solver of type requested -- never happen", ERROR);
    break;
  }
}


bool MatrixSys::match (const real_t     lambda2,
		       const real_t     betak2 ,
		       const NumberSys* nScheme,
		       const SolverKind method ) const
// ---------------------------------------------------------------------------
// The unique identifiers of a MatrixSys are presumed to be given
// by the constants and the numbering system used.  Other things that
// could be checked but aren't (yet) include geometric systems and
// quadrature schemes.
// ---------------------------------------------------------------------------
{
  if (fabs (_HelmholtzConstant - lambda2) < EPSDP                       &&
      fabs (_FourierConstant   - betak2 ) < EPSDP                       &&
      _NS -> nGlobal() == nScheme -> nGlobal()                          &&
      _NS -> nSolve()  == nScheme -> nSolve()                           &&
      Veclib::same (_NS->nGlobal(), _NS->btog(), 1, nScheme->btog(), 1) &&
      _method == method                                                  )
    return true;

  else
    return false;
}


MatrixSys::~MatrixSys()
// ---------------------------------------------------------------------------
// Destructor.  Because there may be aliases to the internal vector
// storage we use the family class routines.
// ---------------------------------------------------------------------------
{
  switch (_method) {
  case JACPCG:
    Family::abandon (&_PC);
    break;
  case DIRECT: {
    int_t i;
    for (i = 0; i < _nel; i++) {
      Family::abandon (_hbi + i);
      Family::abandon (_hii + i);
    }
    Family::abandon (&_H);
    delete[] _bipack;
    delete[] _iipack;
  } break;
  default:
    break;
  }
}
