///////////////////////////////////////////////////////////////////////////////
// statistics.C: routines for statistical analysis of AuxFields.
//
// Copyright (c) 1994 <--> $Date: 2015/04/20 11:14:18 $, Hugh Blackburn
//
// The collection of statistics is controlled by the setting of the
// AVERAGE token. Legal values are 0 (default), 1, 2, 3. The routines
// here do not control how often statistics are updated: that happens
// in Analyser class methods.
//
// All collected statistics are given 1-character names. Potentially
// these may clash with quantities derived/used by addfield utility
// for example. Caveat emptor. In the longer term we could move to
// strings for field names.
//
// AVERAGE = 0.  No statistics.
//
// AVERAGE = 1.  Running averages of variables held by the Domain used
// for initialisation. The data are held in semi-Fourier state.
//
// AVERAGE = 2.  Additionally, correlation terms for computation of
// Reynolds stresses. (Correlations, based on products of variables,
// are computed and held in physical space.)
// 
// Naming scheme for components of the symmetric "Reynolds stresses" tensor:
//
//                      / uu uv uw \     /  A  B  D \     /  A  B \
//                      | .  vv vw |  =  |  .  C  E |  =  \  .  C /  -- if 2C
//                      \ .  .  ww /     \  .  .  F /
//
// What is computed are the running average of the products uu, uv,
// etc, which are NOT the actual Reynolds stresses: they need to have
// the products of the mean values UU, UV etc subtracted, assumed to
// occur in postprocessing. If there are only two velocity components
// present in the initialising Domain, only averages for A, B & C are
// made.
//
// AVERAGE = 3. Further additional correlations are kept for
// computation of energy equation terms. Again, the correct terms need
// to be made in post-processing.
//
// a) Scalar: 
//    i) q = 0.5 [u^2 + v^2 (+ w^2)]
//   ii) d = SijSij
//
// b) Vector: Naming:
//    i) p u_i
//                      / pu \   / m \   / m \
//                      | pv | = | n | = \ n /  -- if 2C
//                      \ pw /   \ o /
//   ii) q u_i
//                      / qu \   / r \   / r \
//                      | qv | = | s | = \ s /  -- if 2C
//                      \ qw /   \ t /
//
//   iii) Sij u_j       / SxxU + SxyV + SxzW \   / a \   / a \
//                      | SyxU + SyyV + SyzW | = | b | = \ b /  -- if 2C
//                      \ SzxU + SzyV + SzzW /   \ c /
// 
// c) Tensor: symmetric rate-of-strain tensor S_ij. Naming:
//
//                      / xx xy xz \     /  G  H  J \     /  G  H \
//                      | .  yy yz |  =  |  .  I  K |  =  \  .  I /  -- if 2C
//                      \ .  .  zz /     \  .  .  L /
//
// NB: The veracity of the energy equation terms has not been checked
// for cylindrical coordinates. They are probably OK provided the
// domain is invariant in the axial direction (e.g. a straight tube).
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

static char RCS[] = "$Id: statistics.cpp,v 8.1 2015/04/20 11:14:18 hmb Exp $";

#include <sem.h>


Statistics::Statistics (Domain* D) :
// ---------------------------------------------------------------------------
// Store averages for all Domain Fields, and correlations.
// ---------------------------------------------------------------------------
  _name (D -> name),
  _base (D),
  _iavg (Femlib::ivalue ("AVERAGE")),
  _nraw (_base -> nField()),
  _nvel (_nraw - 1),
  _nrey ((_iavg > 1) ? ((_nvel+1)*_nvel)/2 : 0),
  _neng (0)
{
  if (_iavg == 0) return;
  if ((_iavg  < 0) || (_iavg > 3))
    message ("Statistics::Statistics", "AVERAGE token out of [0,3]", ERROR);
					 
  int_t       i;
  const int_t nz   = Geometry::nZProc();
  const int_t ntot = Geometry::nTotProc();

  // -- Set pointers, allocate storage.

  for (i = 0; i < _nraw; i++)	// -- Local pointers to raw variables.
    _raw[_base -> u[i] -> name()] = (AuxField*) _base -> u[i];

  if (_iavg > 0) // -- Set up buffers for averages of raw variables.
    for (i = 0; i < _nraw; i++)
      _avg[_base -> u[i] -> name()] =
	new AuxField (new real_t[ntot],nz,_base->elmt,_base->u[i]->name());

  if (_iavg > 1) // -- Set up buffers for Reynolds stress correlations.
    for (i = 0; i < _nrey; i++)
      _avg['A' + i] = new AuxField (new real_t[ntot],nz,_base->elmt,'A'+i);

  if (_iavg > 2) { // -- Set up addtional buffers for energy correlations.

    // -- Scalar.

    _avg['q'] = new AuxField (new real_t[ntot],nz,_base->elmt,'q'); ++_neng;
    _avg['d'] = new AuxField (new real_t[ntot],nz,_base->elmt,'d'); ++_neng;

    // -- Vector.

    _avg['m'] = new AuxField (new real_t[ntot],nz,_base->elmt,'m'); ++_neng;
    _avg['n'] = new AuxField (new real_t[ntot],nz,_base->elmt,'n'); ++_neng;
    if (_nvel == 3) {
      _avg['o'] = new AuxField (new real_t[ntot],nz,_base->elmt,'o'); ++_neng;
    }

    _avg['r'] = new AuxField (new real_t[ntot],nz,_base->elmt,'r'); ++_neng;
    _avg['s'] = new AuxField (new real_t[ntot],nz,_base->elmt,'s'); ++_neng;
    if (_nvel == 3) {
      _avg['t'] = new AuxField (new real_t[ntot],nz,_base->elmt,'t'); ++_neng;
    }

    _avg['a'] = new AuxField (new real_t[ntot],nz,_base->elmt,'a'); ++_neng;
    _avg['b'] = new AuxField (new real_t[ntot],nz,_base->elmt,'b'); ++_neng;
    if (_nvel == 3) {
      _avg['c'] = new AuxField (new real_t[ntot],nz,_base->elmt,'c'); ++_neng;
    }
      
    // -- Tensor.
    
    _avg['G'] = new AuxField (new real_t[ntot],nz,_base->elmt,'G'); ++_neng;
    _avg['H'] = new AuxField (new real_t[ntot],nz,_base->elmt,'H'); ++_neng;
    _avg['I'] = new AuxField (new real_t[ntot],nz,_base->elmt,'I'); ++_neng;

    if (_nvel == 3) {
      _avg['J'] = new AuxField (new real_t[ntot],nz,_base->elmt,'J'); ++_neng;
      _avg['K'] = new AuxField (new real_t[ntot],nz,_base->elmt,'K'); ++_neng;
      _avg['L'] = new AuxField (new real_t[ntot],nz,_base->elmt,'L'); ++_neng;
    }
  }
}


void Statistics::initialise (const char* filename)
// ---------------------------------------------------------------------------
// This is for standard running averages. Try to initialize from file
// filename (e.g. "session.avg"), failing that set all buffers to
// zero.  Number of fields in file should be same as Statistics::_avg
// buffer.
// ---------------------------------------------------------------------------
{
  ROOTONLY cout << "-- Initialising averaging  : ";  

  ifstream file (filename);
  map<char, AuxField*>::iterator k;

  if (file) {
    ROOTONLY {
      cout << "read from file " << filename;
      cout.flush();
    }
    file >> *this;
    file.close();

    // -- Fourier transform raw data components.

    for (k = _raw.begin(); k != _raw.end(); k++)
      _avg[k -> second -> name()] -> transform (FORWARD);
  
  } else {			// -- No file, set to zero.
    ROOTONLY cout << "set to zero";
    for (k = _avg.begin(); k != _avg.end(); k++) *(k -> second) = 0.0;
    _navg = 0;
  }

  ROOTONLY cout << endl;
}


void Statistics::update (AuxField** wrka,
			 AuxField** wrkb)
// ---------------------------------------------------------------------------
// Update running averages, using arrays wrka & wrkb as workspace.
// All product/correlation terms are calculated without dealiasing,
// and are held in physical space.
// ---------------------------------------------------------------------------
{
  if (_iavg < 1) return;
  
  char   key;
  int_t  i, j;
  Field* master = _base -> u[0];
  map<char, AuxField*>::iterator k;

  // -- Weight old running averages.

  for (k = _avg.begin(); k != _avg.end(); k++)
    *(k -> second) *= static_cast<real_t>(_navg);
    
  // -- Always do running averages of raw data.

  for (k = _raw.begin(); k != _raw.end(); k++)
    *_avg[k -> second-> name()] += *(k -> second);

  // -- Reynolds stress correlations.
  //    After this, wrka contains current velocity data in physical space.

  if (_iavg > 1) {
    (*wrka[0] = *_raw['u']) . transform (INVERSE);
    (*wrka[1] = *_raw['v']) . transform (INVERSE);
    if (_nvel == 3) (*wrka[2] = *_raw['w']) . transform (INVERSE);

    _avg['A'] -> timesPlus (*wrka[0], *wrka[0]);
    _avg['B'] -> timesPlus (*wrka[0], *wrka[1]);
    _avg['C'] -> timesPlus (*wrka[1], *wrka[1]);
    
    if (_nvel == 3) {
      _avg['D'] -> timesPlus (*wrka[0], *wrka[2]);
      _avg['E'] -> timesPlus (*wrka[1], *wrka[2]);
      _avg['F'] -> timesPlus (*wrka[2], *wrka[2]);
    }
  }

  // -- Additional working for energy terms.

  if (_iavg > 2) {

    // -- Pressure--velocity terms.

    (*wrkb[0] = *_raw['p']) . transform (INVERSE);
    _avg['m'] -> timesPlus (*wrkb[0], *wrka[0]);
    _avg['n'] -> timesPlus (*wrkb[0], *wrka[1]);
    if (_nvel == 3)
      _avg['o'] -> timesPlus (*wrkb[0], *wrka[2]);

    // -- q, TKE.

    wrkb[0] -> times (*wrka[0], *wrka[0]);
    for (i = 1; i < _nvel; i++)
      wrkb[0] -> timesPlus (*wrka[i], *wrka[i]);
    *wrkb[0] *= 0.5;
    *_avg['q'] += *wrkb[0];

    // -- q u_i.

    _avg['r'] -> timesPlus (*wrka[0], *wrkb[0]);
    _avg['s'] -> timesPlus (*wrka[1], *wrkb[0]);
    if (_nvel == 3)
      _avg['t'] -> timesPlus (*wrka[2], *wrkb[0]);

    // -- Strain-rate tensor (see also eddyvis.C in les-smag).

    if (Geometry::cylindrical()) { // -- see Bird Stewart & Lightfoot.

      // -- Off-diagonal terms.

      AuxField* tp1 = wrka[0];
      AuxField* tp2 = wrka[1];
  
      for (i = 0; i < _nvel; i++)
	for (j = 0; j < _nvel; j++) {
	  if (j == i) continue;
	  if (i == 2 && j == 1) {
	    (*tp1 = *_raw['w']) . gradient (1);
	    (*tp2 = *_raw['w']) . divY();
	    *tp1 -= *tp2;
	  } else {
	    (*tp1 = *_raw['u' + i]) . gradient (j);
	    if (j == 2) tp1 -> divY();
	  }
	  if   (j > i) *wrkb[i + j - 1]  = *tp1;
	  else         *wrkb[i + j - 1] += *tp1;
	}
  
      for (i = 0; i < _nvel; i++) *wrkb[i] *= 0.5;
      
      // -- Diagonal.
      
      for (i = 0; i < _nvel; i++) {
	(*wrka[i] = *_raw['u' + i]) . gradient (i);
	if (i == 2) (*wrka[2] += *_raw['v']) . divY();
      }

    } else {			// -- Cartesian geometry.
      
      // -- Off-diagonal terms.

      AuxField* tmp = wrka[0];

      for (i = 0; i < _nvel; i++)
	for (j = 0; j < _nvel; j++) {
	  if (j == i) continue;
	  (*tmp = *_raw['u' + i]) . gradient (j);
	  if   (j > i) *wrkb[i + j - 1]  = *tmp;
	  else         *wrkb[i + j - 1] += *tmp;
	}
      
      for (i = 0; i < _nvel; i++) *wrkb[i] *= 0.5;

      // -- Diagonal.

      for (i = 0; i < _nvel; i++)
	(*wrka[i] = *_raw['u' + i]) . gradient (i);
    }

    // -- Bring strain rate tensor components into physical space.

    wrka[0] -> transform (INVERSE);
    wrkb[0] -> transform (INVERSE);
    wrka[1] -> transform (INVERSE);
    if (_nvel == 3) {
      wrkb[1] -> transform (INVERSE);
      wrkb[2] -> transform (INVERSE);
      wrka[2] -> transform (INVERSE);
    }

    // -- Add strain rate tensor components into running averages.

    *_avg['G'] += *wrka[0];
    *_avg['H'] += *wrkb[0];
    *_avg['I'] += *wrka[1];
    if (_nvel == 3) {
      *_avg['J'] += *wrkb[1];
      *_avg['K'] += *wrkb[2];
      *_avg['L'] += *wrka[2];
    }

    // -- Compute the strain-velocity correlation.

    _raw['u'] -> transform (INVERSE);
    _raw['v'] -> transform (INVERSE);
    if (_nvel == 3)
      _raw['w'] -> transform (INVERSE);

    _avg['a'] -> timesPlus (*wrka[0], *_raw['u']);
    _avg['a'] -> timesPlus (*wrkb[0], *_raw['v']);
    if (_nvel == 3)
      _avg['a'] -> timesPlus (*wrkb[1], *_raw['w']);

    _avg['b'] -> timesPlus (*wrkb[0], *_raw['u']);
    _avg['b'] -> timesPlus (*wrka[1], *_raw['v']);
    if (_nvel == 3)
      _avg['b'] -> timesPlus (*wrkb[2], *_raw['w']);

    if (_nvel == 3) {
      _avg['c'] -> timesPlus (*wrkb[1], *_raw['u']);
      _avg['c'] -> timesPlus (*wrkb[2], *_raw['v']);
      _avg['c'] -> timesPlus (*wrka[2], *_raw['w']);
    }

    _raw['u'] -> transform (FORWARD);
    _raw['v'] -> transform (FORWARD);
    if (_nvel == 3)
      _raw['w'] -> transform (FORWARD);

    // -- Compute fully-contracted strain rate scalar, d.

    *wrkb[0] *= sqrt (2.0);
    if (_nvel == 3) {
      *wrkb[1] *= sqrt (2.0);
      *wrkb[2] *= sqrt (2.0);
    }

    _avg['d'] -> timesPlus (*wrka[0], *wrka[0]);
    _avg['d'] -> timesPlus (*wrkb[0], *wrkb[0]);
    _avg['d'] -> timesPlus (*wrka[1], *wrka[1]);
    if (_nvel == 3) {
      _avg['d'] -> timesPlus (*wrkb[1], *wrkb[1]);
      _avg['d'] -> timesPlus (*wrkb[2], *wrkb[2]);
      _avg['d'] -> timesPlus (*wrka[2], *wrka[2]);
    }
  }

  // -- Normalise and smooth running averages.

  for (k = _avg.begin(); k != _avg.end(); k++) {
    master -> smooth (k -> second);
    *(k -> second) /= static_cast<real_t>(_navg + 1);
  }

  _navg++;
}


void Statistics::dump (const char* filename)
// ---------------------------------------------------------------------------
// Similar to Domain::dump.
//
// As of 24/11/2004, we deleted the checkpointing that used to happen:
// all dumping now happens to file named on input.
//
// We also smooth all the outputs with the mass matrix.
// ---------------------------------------------------------------------------
{
  const int_t step     = _base -> step;
  const bool  periodic = !(step %  Femlib::ivalue ("IO_FLD"));
  const bool  initial  =   step == Femlib::ivalue ("IO_FLD");
  const bool  final    =   step == Femlib::ivalue ("N_STEP");

  if (!(periodic || final)) return;

  ofstream    output;
  int_t       i;
  map<char, AuxField*>::iterator k;

  ROOTONLY {
    const char routine[] = "Statistics::dump";
    const bool verbose   = static_cast<bool> (Femlib::ivalue ("VERBOSE"));

    output.open (filename);
    if (!output) message (routine, "can't open dump file", ERROR);
    if (verbose) message (routine, ": writing field dump", REMARK);
  }
  
  // -- All terms are written out in physical space but some are
  //    held internally in Fourier space.

  for (k = _raw.begin(); k != _raw.end(); k++)
    _avg[k -> second -> name()] -> transform (INVERSE);

  if (_neng) {
    _avg['G'] -> transform (INVERSE);
    _avg['H'] -> transform (INVERSE);
    _avg['I'] -> transform (INVERSE);
    if (_nvel == 3) {
      _avg['J'] -> transform (INVERSE);
      _avg['K'] -> transform (INVERSE);
      _avg['L'] -> transform (INVERSE);
    }
  }

  output << *this;

  for (k = _raw.begin(); k != _raw.end(); k++)
    _avg[k -> second -> name()] -> transform (FORWARD);

  if (_neng) {
    _avg['G'] -> transform (FORWARD);
    _avg['H'] -> transform (FORWARD);
    _avg['I'] -> transform (FORWARD);
    if (_nvel == 3) {
      _avg['J'] -> transform (FORWARD);
      _avg['K'] -> transform (FORWARD);
      _avg['L'] -> transform (FORWARD);
    }
  }

  ROOTONLY output.close();
}


ofstream& operator << (ofstream&   strm,
		       Statistics& src )
// ---------------------------------------------------------------------------
// Output Statistics class to file.  Like similar Domain routine.
// ---------------------------------------------------------------------------
{
  int_t             i;
  const int_t       N = src._avg.size();
  vector<AuxField*> field (N);
  map<char, AuxField*>::iterator k;

  for (i = 0, k = src._avg.begin(); i < N; i++, k++) field[i] = k -> second;

  writeField (strm, src._name, src._navg, src._base -> time, field);

  return strm;
}


ifstream& operator >> (ifstream&   strm,
		       Statistics& tgt )
// ---------------------------------------------------------------------------
// Input Statistics class from file.  Like similar Domain routine.
// ---------------------------------------------------------------------------
{
  const char routine[] = "strm>>Statistics";
  int_t i, j, np, nz, nel, ntot, nfields;
  int_t npchk,  nzchk, nelchk;
  char  s[StrMax], f[StrMax], err[StrMax], fields[StrMax];
  bool  swap = false;
  map<char, AuxField*>::iterator k;

  if (strm.getline(s, StrMax).eof()) return strm;
  
  strm.getline (s, StrMax) . getline (s, StrMax);

  string ss(s);
  istringstream sss (ss);
  sss >> np >> np >> nz >> nel;
 
  tgt._avg.begin()->second->describe (f);
  sss.clear ();
  sss.str   (ss = f);
  sss >> npchk >> npchk >> nzchk >> nelchk;

  if (np  != npchk ) message (routine, "element size mismatch",       ERROR);
  if (nz  != nzchk ) message (routine, "number of z planes mismatch", ERROR);
  if (nel != nelchk) message (routine, "number of elements mismatch", ERROR);
  
  ntot = np * np * nz * nel;
  if (ntot != Geometry::nTot())
    message (routine, "declared sizes mismatch", ERROR);

  strm.getline (s, StrMax);

  sss.clear();
  sss.str  (ss = s);
  sss >> tgt._navg;
    
  strm.getline (s, StrMax) . getline (s, StrMax);
  strm.getline (s, StrMax) . getline (s, StrMax) . getline (s, StrMax);
    
  nfields = 0;
  while (isalpha (s[nfields])) {
    fields[nfields] = s[nfields];
    nfields++;
  }
  fields[nfields] = '\0';
  if (nfields != tgt._avg.size()) {
    sprintf (err, "strm: %1d fields, avg: %1d", 
	     nfields, static_cast<int_t>(tgt._avg.size()));
    message (routine, err, ERROR);
  }

  for (i = 0, k = tgt._avg.begin(); k != tgt._avg.end(); k++, i++)
    if (!strchr (fields, k -> second -> name())) {
      sprintf (err, "field %c not present in avg", fields[i]);
      message (routine, err, ERROR);
    }

  strm.getline (s, StrMax);
  Veclib::describeFormat (f);

  if (!strstr (s, "binary"))
    message (routine, "input field strm not in binary format", ERROR);
  
  if (!strstr (s, "endian"))
    message (routine, "input field strm in unknown binary format", WARNING);
  else {
    swap = ((strstr (s, "big") && strstr (f, "little")) ||
	    (strstr (f, "big") && strstr (s, "little")) );
    ROOTONLY {
      if (swap) cout << " (byte-swapping)";
      cout.flush();
    }
  }

  for (j = 0; j < nfields; j++) {
    strm >> *tgt._avg[fields[j]];
    if (swap) tgt._avg[fields[j]] -> reverse();
  }
  
  ROOTONLY if (strm.bad())
    message (routine, "failed reading average file", ERROR);

  return strm;
}


void Statistics::phaseUpdate (const int_t j   ,
			      AuxField**  wrka,
			      AuxField**  wrkb)
// ---------------------------------------------------------------------------
// Phase updates are running updates, like those for standard
// statistics, but they are only computed at (a presumed very limited)
// number of instants per run. And so there would be a number of
// averaging files, to be called e.g. session.0.phs, session.1.phs
// ... session.(N-1).phs. Instead of reserving enough memory for all
// these buffers, we just keep workspace reserved, and at every phase
// point the appropriate file (number j) is read in, updated, and
// written back out. If the file does not exist, it is created.
//
// NB: the number of averages computed (needed for the running
// averaging) is only updated when j == 0, which will happen as the
// last of a set of phase averages.
// ---------------------------------------------------------------------------
{
  char filename[StrMax];
  sprintf (filename, "%s.%1d.phs", _name, j);

  this -> initialise (filename);
  this -> update     (wrka, wrkb);
  this -> dump       (filename);
}
