///////////////////////////////////////////////////////////////////////////////
// particle.C: functions for integrating positions of particles.
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
// See also analysis.C
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: particle.cpp,v 8.1 2015/04/20 11:14:18 hmb Exp $";

#include <sem.h>

Domain* FluidParticle::_Dom     = 0;
int_t   FluidParticle::_NCOM    = 0;
int_t   FluidParticle::_NEL     = 0;
int_t   FluidParticle::_NZ      = 0;
int_t   FluidParticle::_TORD    = 0;
int_t   FluidParticle::_ID_MAX  = 0;
real_t* FluidParticle::_P_coeff = 0;
real_t* FluidParticle::_C_coeff = 0;
real_t* FluidParticle::_Work    = 0;
real_t  FluidParticle::_DT      = 0.0;
real_t  FluidParticle::_Lz      = 0.0;


FluidParticle::FluidParticle (Domain*     d,
			      const int_t i,
			      Point&      p) :
// ---------------------------------------------------------------------------
// Initially particle is located at p.  Find it in the 2D mesh.  Trim to
// periodic length in 3D if required.
// ---------------------------------------------------------------------------
  _id    (i),
  _ctime (d -> time),
  _step  (0),
  _p     (p)
{
  int_t      k;
  const bool guess = true;

  if (!_Dom) {			// -- Set up first time through.
    _Dom     = d;
    _NCOM    = (d -> nField() > 3) ? 3 : 2;
    _NEL     = Geometry::nElmt();
    _NZ      = Geometry::nZ();
    _DT      = Femlib::value  ("D_T");
    _Lz      = Femlib::value  ("TWOPI / BETA");
    _TORD    = Femlib::ivalue ("N_TIME");
    _P_coeff = new real_t [static_cast<size_t>(_TORD + _TORD)];
    _C_coeff = new real_t [static_cast<size_t>(_TORD*(_TORD + 1))];
    _Work    = new real_t [static_cast<size_t>(max 
			(2*Geometry::nTotElmt(), 5 * Geometry::nP() + 6))];

    Veclib::zero (_TORD*_TORD,     _P_coeff, 1);
    Veclib::zero (_TORD*(_TORD+1), _C_coeff, 1);
    
    for (k = 0; k < _TORD; k++) {
      Integration::AdamsBashforth (k+1, _P_coeff+k* _TORD   );
      Integration::AdamsMoulton   (k+2, _C_coeff+k*(_TORD+1));
    }

    Blas::scal (_TORD*_TORD,     _DT, _P_coeff, 1);
    Blas::scal (_TORD*(_TORD+1), _DT, _C_coeff, 1);
  }

  // -- Try to locate particle, stop if can't.

  _E = 0;
  for (k = 0; k < _NEL; k++) {
    _r = _s = 0.0;
    if (_Dom -> elmt[k] -> locate (_p.x, _p.y, _r, _s, _Work, guess)) {
      _E = _Dom -> elmt[k];
      break;
    }
  }

  if (!_E) return;
  if (_id > _ID_MAX) _ID_MAX = _id;

  if (_NCOM == 2) {
    _p.z = 0.0;

    _u = new real_t [static_cast<size_t>(_TORD + _TORD)];
    _v = _u + _TORD;

  } else {
    if   (_p.z < 0.0) _p.z = _Lz - fmod (fabs (_p.z), _Lz);
    else              _p.z = fmod (_p.z, _Lz);
    
    _u = new real_t [static_cast<size_t>(_TORD + _TORD + _TORD)];
    _v = _u + _TORD;
    _w = _v + _TORD;
  }
}


void FluidParticle::integrate ()
// ---------------------------------------------------------------------------
// Integrate massless particle's position using predictor--corrector
// scheme.  If particles leave 2D mesh they are marked by setting E =
// 0.  For 3D, they get put back into the fundamental period of the
// solution if they leave in the z-direction.
//
// NB: The domain velocity fields must be in Fourier space prior to call (3D).
// ---------------------------------------------------------------------------
{
  if (!_E) return;
#if defined (DEBUG)
  const char routine[] = "FluidParticle::integrate";
#endif
  register int_t i;
  const int_t    N     = min (++_step, _TORD);
  const int_t    NP    = N + 1;
  const int_t    NM    = N - 1;
  const bool     guess = true;
  real_t         xp, yp, zp, up, vp, wp;
  real_t         *predictor, *corrector;

  predictor = _P_coeff + NM *  _TORD;
  corrector = _C_coeff + NM * (_TORD+1);

  if (_NCOM == 2) {		// -- 2D integration.
    
    // -- Predictor.

    _u[0] = _Dom -> u[0] -> probe (_E, _r, _s, static_cast<int_t>(0));
    _v[0] = _Dom -> u[1] -> probe (_E, _r, _s, static_cast<int_t>(0));

    xp = _p.x;
    yp = _p.y;
    for (i = 0; i < N; i++) {
      xp += predictor[i] * _u[i];
      yp += predictor[i] * _v[i];
    }

    if (!_E -> locate (xp, yp, _r, _s, _Work)) {
      _E = 0;
      for (i = 0; i < _NEL; i++) {
	_r = _s = 0.0;
	if (_Dom -> elmt[i] -> locate (xp, yp, _r, _s, _Work, guess)) {
	  _E = _Dom -> elmt[i];
	  break;
	}
      }
      if (!_E) {
#if defined (DEBUG)
	if (Femlib::ivalue ("VERBOSE") > 3) {
	  char     str[StrMax];
	  sprintf (str, "Particle %1d at (%f, %f, %f) left mesh",
		   _id, _p.x, _p.y, _p.z);
	  message (routine, str, WARNING);
	}
#endif
	return;
      }
    }

    // -- Corrector.

    up = _Dom -> u[0] -> probe (_E, _r, _s, static_cast<int_t>(0));
    vp = _Dom -> u[1] -> probe (_E, _r, _s, static_cast<int_t>(0));

    _p.x += corrector[0] * up;
    _p.y += corrector[0] * vp;
    for (i = 1; i < NP; i++) {
      _p.x += corrector[i] * _u[i - 1];
      _p.y += corrector[i] * _v[i - 1];
    }

    if (!_E -> locate (_p.x, _p.y, _r, _s, _Work)) {
      _E = 0;
      for (i = 0; i < _NEL; i++) {
	_r = _s = 0.0;
	if (_Dom -> elmt[i] -> locate (_p.x, _p.y, _r, _s, _Work, guess)) {
	  _E = _Dom -> elmt[i];
	  break;
	}
      }
      if (!_E) {
#if defined (DEBUG)
	if (Femlib::ivalue ("VERBOSE") > 3) {
	  char     str[StrMax];
	  sprintf (str, "Particle %1d at (%f, %f, %f) left mesh",
		   _id, _p.x, _p.y, _p.z);
	  message (routine, str, WARNING);
	}
#endif
	return;
      }
    }

    // -- Maintain multilevel storage.

    rollv (_u, _TORD);
    rollv (_v, _TORD);

  } else {			// -- 3D integration.
    
    // -- Predictor.

    _u[0] = _Dom -> u[0] -> probe (_E, _r, _s, _p.z);
    _v[0] = _Dom -> u[1] -> probe (_E, _r, _s, _p.z);
    _w[0] = _Dom -> u[2] -> probe (_E, _r, _s, _p.z);

    xp = _p.x;
    yp = _p.y;
    zp = _p.z;
    for (i = 0; i < N; i++) {
      xp += predictor[i] * _u[i];
      yp += predictor[i] * _v[i];
      zp += predictor[i] * _w[i];
    }

    if (!_E -> locate (xp, yp, _r, _s, _Work)) {
      _E = 0;
      for (i = 0; i < _NEL; i++) {
	_r = _s = 0.0;
	if (_Dom -> elmt[i] -> locate (xp, yp, _r, _s, _Work, guess)) {
	  _E = _Dom -> elmt[i];
	  break;
	}
      }
      if (!_E) return;
    }

    // -- Corrector.

    up = _Dom -> u[0] -> probe (_E, _r, _s, zp);
    vp = _Dom -> u[1] -> probe (_E, _r, _s, zp);
    wp = _Dom -> u[2] -> probe (_E, _r, _s, zp);

    _p.x += corrector[0] * up;
    _p.y += corrector[0] * vp;
    _p.z += corrector[0] * wp;
    for (i = 1; i < NP; i++) {
      _p.x += corrector[i] * _u[i - 1];
      _p.y += corrector[i] * _v[i - 1];
      _p.z += corrector[i] * _w[i - 1];
    }

    if (!_E -> locate (_p.x, _p.y, _r, _s, _Work)) {
      _E = 0;
      for (i = 0; i < _NEL; i++) {
	_r = _s = 0.0;
	if (_Dom -> elmt[i] -> locate (_p.x, _p.y, _r, _s, _Work, guess)) {
	  _E = _Dom -> elmt[i];
	  break;
	}
      }
      if (!_E) return;
    }
    if   (_p.z < 0.0) _p.z = _Lz - fmod (fabs (_p.z), _Lz);
    else              _p.z = fmod (_p.z, _Lz);

    // -- Maintain multilevel storage.

    rollv (_u, _TORD);
    rollv (_v, _TORD);
    rollv (_w, _TORD);
  }
}

