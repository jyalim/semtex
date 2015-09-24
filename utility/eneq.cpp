//////////////////////////////////////////////////////////////////////////////
// eneq.C: from a field file containing correlations, compute terms of
// either the fluctuating flow energy equation, Tennekes & Lumley
// (3.2.1) which is the default action, or the mean flow energy
// equation, T&L (3.1.11). Output values have dimensions of kinetic
// energy per unit mass per unit time, i.e. power per unit mass.

//
// Copyright (c) 2004 <--> $Date: 2015/04/20 11:14:18 $, Hugh Blackburn
//
// NB: the input field file is assumed to contain only velocity and
// appropriate correlation data. The naming conventions employed for
// example in addfield.cpp are broken (names are here used for different
// variables).  The appropriate averages/correlations are collected by
// setting token AVERAGE=3.
//
// The program works for both 2C (uvp) 3C (uvwp) input files.
//
// Also: the reduction has only been confirmed to work for Cartesian,
// not cylindrical, coordinate systems.
//
// Usage:
// -----
// eneq [options] session session.avg
//   options:
//   -h        ... print this message
//   -m        ... output terms for mean flow energy instead of TKE
//
// Output to standard output.
//
// Input field names:
// ------------------
//
// u -- x velocity component (cylindrical: axial)
// v -- y velocity component (cylindrical: radial)
// w -- z velocity component (cylindrical: azimuthal)
// p -- pressure/density
//
// Names for components of the symmetric "Reynolds stress" correlations:
//
//                      / uu uv uw \     /  A  B  D \     /  A  B \
//                      | .  vv vw |  =  |  .  C  E |  =  \  .  C /  -- if 2C
//                      \ .  .  ww /     \  .  .  F /
//
// Names for correlations specific to the energy equation:
// 
// a) Scalar: 
//    i) q = 0.5 [u^2 + v^2 (+ w^2)]
//   ii) d = sijsij
//
// b) Vector:
//    i) p u_i
//                      / pu \   / m \     / m \
//                      | pv | = | n |  =  \ n /  -- if 2C
//                      \ pw /   \ o /
//   ii) q u_i
//                      / qu \   / r \     / r \
//                      | qv | = | s |  =  \ s /  -- if 2C
//                      \ qw /   \ t /
//
//   iii) Sij u_j       / SxxU + SxyV + SxzW \   / a \   / a \
//                      | SyxU + SyyV + SxzW | = | b | = \ b /  -- if 2C
//                      \ SzxU + SzyV + SzzW /   \ c /
// 
// c) Tensor: symmetric rate-of-strain tensor sij:
//
//                      / xx xy xz \     /  G  H  J \     /  G  H \
//                      | .  yy yz |  =  |  .  I  K |  =  \  .  I /  -- if 2C
//                      \ .  .  zz /     \  .  .  L /
//
// Names for (output) terms in the fluctutating energy equation (TKE):
// -------------------------------------------------------------------
// In the case that we consider non-Newtonian flow, there are ten
// terms to compute.  For Newtonian flow we only need six.  We will
// use arabic numeric characters to name these terms.
//
// '0': ui.uj.Sij                (turbulent) Production
// '1': Uj dq/dxj                Mean flow transport
// '2': d(q.uj)/dxj              Turbulent transport
// '3': d(1/rho p.uj)/dxj        Pressure work
// '4': -2 kinvis d(sij.ui)/dxj  Viscous transport
// '7': 2 kinvis sij.sij         Dissipation
//
// Note: the "missing" names 5,6,8,9 are present for non-Newtonian
// flow.  The sum of terms 1,2,3,4,7,0 should be 0 for stationary
// turbulence, so we also write out
//
// S: 0 + 1 + 2 + 3 + 4 + 7
//
// NEW OUTPUT NAMES (to get around tecplot's problems with numeric names);
//
// S: a + b + c + d + e + h
//
// which should approach zero pointwise, and when integrated over domain.
//
// Names for (output) terms in the MEAN energy equation (MKE):
// -----------------------------------------------------------
//
// '0': -ui.uj.Sij                (turbulent) Production
// '1': Uj dQ/dxj                 Mean flow transport, Q = 0.5*UiUi
// '2': d(ui.uj.Ui)/dxj           Reynolds stress transport                  
// '3': 1/rho*d(P.Uj)/dxj         Pressure work
// '4': -2 kinvis d(Sij.Ui)/dxj   Viscous transport
// '7': 2 kinvis Sij.Sij          Mean flow dissipation
//
// Again, the sum of terms 1,2,3,4,7,0 should equal 0 for stationary
// turbulence, so we also write out
//
// S: 0 + 1 + 2 + 3 + 4 + 7
//
// NEW OUTPUT NAMES (to get around tecplot's problems with numeric names);
//
// S: a + b + c + d + e + h
//
// which should approach zero pointwise, and when integrated over
// domain. For steady flow, terms 2 and 0 are also zero.
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
// 02110-1301 USA
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: eneq.cpp,v 8.1 2015/04/20 11:14:18 hmb Exp $";

#include <sem.h>

static char prog[] = "eneq";
static void getargs (int,char**, bool&, const char*&, const char*&);
static void getMesh (const char*,vector<Element*>&);
static void makeBuf (map<char,AuxField*>&,vector<AuxField*>&,
		     vector<Element*>&);
static bool getDump (ifstream&,map<char, AuxField*>&,vector<Element*>&);
static bool doSwap  (const char*);
static void covary  (map<char,AuxField*>&,map<char,AuxField*>&,
		     vector<AuxField*>&);
static void meanflo (map<char,AuxField*>&,map<char,AuxField*>&,
		     vector<AuxField*>&);

static const char* fieldNames(map<char, AuxField*>&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver -- adapted from probe.C.
// ---------------------------------------------------------------------------
{
  const char           *session, *dump;
  ifstream             avgfile;
  vector<Element*>     elmt;
  map<char, AuxField*> input, output;
  vector<AuxField*>    work,  outbuf;
  bool                 MKE = false;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, MKE, session, dump);

  avgfile.open (dump, ios::in);
  if (!avgfile) message (prog, "no field file", ERROR);
  
  getMesh (session, elmt);
  makeBuf (output, work, elmt);

  // -- Need to link outbuf to output so we can use writeField.
  //    Maybe in the longer term we should overload writeField.

  outbuf.resize (output.size()); int_t i = 0;
  for (map<char,AuxField*>::iterator k = output.begin();
       k != output.end(); k++, i++) outbuf[i] = k -> second;
  
  while (getDump (avgfile, input, elmt)) {
    if (MKE) meanflo (input, output, work);
    else     covary  (input, output, work);
    writeField (cout, session, 0, 0.0, outbuf);
  }
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int          argc   ,
		     char**       argv   ,
		     bool&        MKE    ,
		     const char*& session,
		     const char*& dump   )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] =
    "Usage: %s [options] session dump.avg\n"
    "options:\n"
    "  -h ... print this message \n"
    "  -m ... output terms for mean flow energy instead of TKE\n";
              
  char buf[StrMax];
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog); cerr << buf; exit (EXIT_SUCCESS);
      break;
    case 'm':
      MKE = true;
      break;
    default:
      sprintf (buf, usage, prog); cerr << buf; exit (EXIT_FAILURE);
      break;
    }

  if (argc != 2) {
    sprintf (buf, usage, prog); cerr << buf; exit (EXIT_FAILURE);
  } else {
    --argc; session = *argv;
    --argc; dump    = *++argv;
  }
}


static void getMesh (const char*       session,
		     vector<Element*>& elmt   )
// ---------------------------------------------------------------------------
// Set up 2D mesh information. Note that parser tokens and Geometry
// are set here, too.
// ---------------------------------------------------------------------------
{
  FEML* F = new FEML (session);
  Mesh* M = new Mesh (F);
  
  const int_t nel = M -> nEl();  
  const int_t np  = Femlib::ivalue ("N_P");
  const int_t nz  = Femlib::ivalue ("N_Z");

  Geometry::CoordSys space = (Femlib::ivalue ("CYLINDRICAL")) ?
    Geometry::Cylindrical : Geometry::Cartesian;
  
  Geometry::set (np, nz, nel, space);
  elmt.resize   (nel);

  for (int_t k = 0; k < nel; k++) elmt[k] = new Element (k, np, M);
}


static void makeBuf (map<char, AuxField*>& output,
		     vector<AuxField*>&    work  ,
		     vector<Element*>&     elmt  )
// ---------------------------------------------------------------------------
// Note that we only set up the output and work buffers here. The
// input buffers get created in getDump.
// 
// NB: the numeric term names seem to give Tecplot some problems.  Fix?
// Or maybe it's preplot?
// ---------------------------------------------------------------------------
{
  const int_t nz   = Geometry::nZ();
  const int_t ntot = Geometry::nTotal();

  work.resize (2);

  output['0'] = new AuxField (new real_t[ntot], nz, elmt, 'a');
  output['1'] = new AuxField (new real_t[ntot], nz, elmt, 'b');
  output['2'] = new AuxField (new real_t[ntot], nz, elmt, 'c');
  output['3'] = new AuxField (new real_t[ntot], nz, elmt, 'd');
  output['4'] = new AuxField (new real_t[ntot], nz, elmt, 'e');
  output['7'] = new AuxField (new real_t[ntot], nz, elmt, 'h');
  output['S'] = new AuxField (new real_t[ntot], nz, elmt, 'S');

  work[0]     = new AuxField (new real_t[ntot], nz, elmt, '\0');
  work[1]     = new AuxField (new real_t[ntot], nz, elmt, '\0');
}


static const char* fieldNames (map<char, AuxField*>& u)
// ---------------------------------------------------------------------------
// Return string containing single-character names of fields.
// ---------------------------------------------------------------------------
{
  int_t i;
  static char buf[StrMax];
  map<char, AuxField*>::iterator k;

  for (i = 0, k = u.begin(); k != u.end(); k++, i++) buf[i] = k -> first;
  buf[i] = '\0';

  return buf;
}

static bool getDump (ifstream&             file,
		     map<char, AuxField*>& u   ,
		     vector<Element*>&     elmt)
// ---------------------------------------------------------------------------
// Load data from field dump, with byte-swapping if required.
// If there is more than one dump in file, it is required that the
// structure of each dump is the same as the first.
// ---------------------------------------------------------------------------
{
  const int_t ntot = Geometry::nTotal();
  char        buf[StrMax], fields[StrMax];
  int_t       i, nf, np, nz, nel;
  bool        swab;

  if (file.getline(buf, StrMax).eof()) return false;
  
  if (!strstr (buf, "Session")) message (prog, "not a field file", ERROR);
  file.getline (buf, StrMax);

  // -- Input numerical description of field sizes.

  file >> np >> nz >> nz >> nel;
  file.getline (buf, StrMax);
  
  if (np != Geometry::nP() || nz != Geometry::nZ() || nel != Geometry::nElmt())
    message (prog, "size of dump mismatch with session file", ERROR);

  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);

  // -- Input field names, assumed to be written without intervening spaces.

  file >> fields;
  nf = strlen  (fields);
  file.getline (buf, StrMax);

  // -- Arrange for byte-swapping if required.

  file.getline  (buf, StrMax);
  swab = doSwap (buf);

  // -- Create AuxFields on first pass.

  if (u.size() == 0) {
    for (i = 0; i < nf; i++)
      u[fields[i]] = new AuxField (new real_t[ntot], nz, elmt, fields[i]);
  } else if (strcmp (fieldNames (u), fields) != 0)
    message (prog, "fields mismatch with first dump in file", ERROR);

  // -- Read binary field data.

  for (map<char, AuxField*>::iterator k = u.begin(); k != u.end(); k++) {
    file >> *(k -> second);
    if (swab) k -> second -> reverse();
  }

  return file.good();
}


static bool doSwap (const char* ffmt)
// ---------------------------------------------------------------------------
// Figure out if byte-swapping is required to make sense of binary input.
// ---------------------------------------------------------------------------
{
  char mfmt[StrMax];

  Veclib::describeFormat (mfmt);   

  if (!strstr (ffmt, "binary"))
    message (prog, "input field file not in binary format", ERROR);
  else if (!strstr (ffmt, "endian"))
    message (prog, "input field file in unknown binary format", WARNING);

  return (strstr (ffmt, "big") && strstr (mfmt, "little")) || 
         (strstr (mfmt, "big") && strstr (ffmt, "little"));
}


static void meanflo (map<char, AuxField*>& in  ,
		     map<char, AuxField*>& out ,
		     vector<AuxField*>&    work)
// ---------------------------------------------------------------------------
// This does the actual work of building MKE equation terms.
// ---------------------------------------------------------------------------
{
  const char   list2d[] = "ABCGHIabdmnpqrsuv";
  const char   list3d[] = "ABCDEFGHIJKLabcdmnopqrstuvw";
  const char*  names    = fieldNames (in);
  const int_t  nvel     = (strchr (names, 'w')) ? 3 : 2;
  const real_t kinvis   = Femlib::value ("KINVIS");
  char         err[StrMax];

  if (nvel == 2) {
    if (strcmp (names, list2d) != 0) {
      sprintf (err,"list of names should match %s: have %s", list2d, names);
      message (prog, err, ERROR);
    } 
  } else if (strcmp (names, list3d) != 0) {
    sprintf (err,"list of names should match %s: have %s", list3d, names);
    message (prog, err, ERROR);
  }

  // -- Turn ABC... into Reynolds stresses.

  in['A'] -> timesMinus (*in['u'], *in['u']);
  in['B'] -> timesMinus (*in['u'], *in['v']);
  in['C'] -> timesMinus (*in['v'], *in['v']);
  if (nvel == 3) {
    in['D'] -> timesMinus (*in['u'], *in['w']);
    in['E'] -> timesMinus (*in['v'], *in['w']);
    in['F'] -> timesMinus (*in['w'], *in['w']);
  }

  // -- At this stage, in['q'] still holds the total kinetic energy.
  //    Replace it with the mean flow kinetic energy.
 
  in['q'] -> times     (*in['u'], *in['u']);
  in['q'] -> timesPlus (*in['v'], *in['v']);
  if (nvel == 3)
    in['q'] -> timesPlus (*in['w'], *in['w']);
  *in['q'] *= 0.5;

  // -- Next we build the mean rate-of-strain tensor for further work.
  //    The '2d' components get stored in 'm', 'n', 'r' while 
  //    the '3d' ones, if needed,   go in 'o', 's', 't'.

  (*in['m']  = *in['u']) . gradient (0);
  (*in['r']  = *in['v']) . gradient (1);
  (*in['n']  = *in['u']) . gradient (1);
  (*work[0]  = *in['v']) . gradient (0);
  (*in['n'] += *work[0]) *= 0.5;

  if (nvel == 3) {
    (*in['o']  = *in['u']).transform(FORWARD).gradient(2).transform(INVERSE);
    (*work[0]  = *in['w']).gradient(0);
    (*in['o'] += *work[0]) *= 0.5;
    (*in['s']  = *in['v']).transform(FORWARD).gradient(2).transform(INVERSE);
    (*work[0]  = *in['w']).gradient(1);
    (*in['s'] += *work[0]) *= 0.5;
    (*in['t']  = *in['w']).transform(FORWARD).gradient(2).transform(INVERSE);
  }

  // -- So let's get started by making term '1' (advection):

  *work[0]  = *in['q']; work[0]   -> gradient (0);
  *work[0] *= *in['u']; *out['1']  = *work[0];
  *work[0]  = *in['q']; work[0]   -> gradient (1);
  *work[0] *= *in['v']; *out['1'] += *work[0];
  if (nvel == 3) {
   (*work[0]  = *in['q']).transform(FORWARD).gradient(2).transform(INVERSE);
    *work[0] *= *in['w']; *out['1'] += *work[0];
  }

  // -- Term '2' (turbulent transport, zero for steady flow). 
  //    Destroy 'a', 'b', 'c' along the way:

  in['a'] -> times     (*in['A'], *in['u']);
  in['a'] -> timesPlus (*in['B'], *in['v']);
  in['b'] -> times     (*in['B'], *in['u']);
  in['b'] -> timesPlus (*in['C'], *in['v']);
  if (nvel == 3) {
    in['a'] -> timesPlus (*in['D'], *in['w']);
    in['b'] -> timesPlus (*in['E'], *in['w']);
    in['c'] -> times     (*in['D'], *in['u']);
    in['c'] -> timesPlus (*in['E'], *in['v']);
    in['c'] -> timesPlus (*in['F'], *in['w']);
  }

  in['a'] -> gradient (0); *out['2']  = *in['a'];
  in['b'] -> gradient (1); *out['2'] += *in['b'];
  if (nvel == 3) {
    (in['c']  -> transform (FORWARD)) . gradient (2) . transform (INVERSE);
    *out['2'] += *in['c'];
  }

#if 0
  *out['0'] = *in['A'];
  *out['1'] = *in['B'];
  *out['2'] = *in['C'];
  *out['3'] = *in['a'];
  *out['4'] = *in['b'];
#else

  // -- Term '3' (pressure work).

  out['3'] -> times (*in['p'], *in['u']);
  out['3'] -> gradient (0);
  work[0]  -> times (*in['p'], *in['v']);
  work[0]  -> gradient (1); *out['3'] += *work[0];
  if (nvel == 3) {
    work[0]  -> times (*in['p'], *in['w']);
    (work[0] ->  transform (FORWARD)) . gradient (2) . transform (INVERSE);
    *out['3'] += *work[0];
  }

  // -- Compute term '4' (viscous transport), somewhat like term '2'.

  in['a'] -> times     (*in['m'], *in['u']);
  in['a'] -> timesPlus (*in['n'], *in['v']);
  in['b'] -> times     (*in['n'], *in['u']);
  in['b'] -> timesPlus (*in['r'], *in['v']);
  if (nvel == 3) {
    in['a'] -> timesPlus (*in['o'], *in['w']);
    in['b'] -> timesPlus (*in['s'], *in['w']);
    in['c'] -> times     (*in['o'], *in['u']);
    in['c'] -> timesPlus (*in['s'], *in['v']);
    in['c'] -> timesPlus (*in['t'], *in['w']);
  }

  in['a'] -> gradient (0); *out['4']  = *in['a'];
  in['b'] -> gradient (1); *out['4'] += *in['b'];
  if (nvel == 3) {
    (in['c']  -> transform (FORWARD)) . gradient (2) . transform (INVERSE);
    *out['4'] += *in['c'];
  }
  
  *out['4'] *= -2.0 * kinvis;
  
  // -- Compute term '7' (mean flow dissipation):

  *in['n'] *= sqrt (2.0);
  if (nvel == 3) {
    *in['o'] *= sqrt (2.0);
    *in['s'] *= sqrt (2.0);
  }

  out['7'] -> times     (*in['m'], *in['m']);
  out['7'] -> timesPlus (*in['n'], *in['n']);
  out['7'] -> timesPlus (*in['r'], *in['r']);
  if (nvel == 3) {
    out['7'] -> timesPlus (*in['o'], *in['o']);
    out['7'] -> timesPlus (*in['s'], *in['s']);
    out['7'] -> timesPlus (*in['t'], *in['t']);
  }
  
  *out['7'] *= 2.0 * kinvis;

  // -- Compute term '0' (turbulence production, zero for steady flow):

  *in['n'] *= sqrt (2.0);
  if (nvel == 3) {
    *in['r'] *= sqrt (2.0);
    *in['s'] *= sqrt (2.0);
  }
 
  out['0'] -> times     (*in['A'], *in['m']);
  out['0'] -> timesPlus (*in['B'], *in['n']);
  out['0'] -> timesPlus (*in['C'], *in['r']);
  if (nvel == 3) {
    out['0'] -> timesPlus (*in['D'], *in['o']);
    out['0'] -> timesPlus (*in['E'], *in['s']);
    out['0'] -> timesPlus (*in['F'], *in['t']);
  }

  *out['0'] *= -1.0;

  // -- Finally, compute the sum, 'S' (should converge to zero):

 (*out['S']  = *out['1']);
  *out['S'] += *out['2'];
  *out['S'] += *out['3'];
  *out['S'] += *out['4'];
  *out['S'] += *out['7'];
  *out['S'] += *out['0'];
#endif
}


static void covary  (map<char, AuxField*>& in  ,
		     map<char, AuxField*>& out ,
		     vector<AuxField*>&    work)
// ---------------------------------------------------------------------------
// Build TKE equation terms.
// ---------------------------------------------------------------------------
{
  const char   list2d[] = "ABCGHIabdmnpqrsuv";
  const char   list3d[] = "ABCDEFGHIJKLabcdmnopqrstuvw";
  const char*  names    = fieldNames (in);
  const int_t  nvel     = (strchr (names, 'w')) ? 3 : 2;
  const real_t kinvis   = Femlib::value ("KINVIS");
  char         err[StrMax];

  if (nvel == 2) {
    if (strcmp (names, list2d) != 0) {
      sprintf (err,"list of names should match %s: have %s", list2d, names);
      message (prog, err, ERROR);
    } 
  } else if (strcmp (names, list3d) != 0) {
    sprintf (err,"list of names should match %s: have %s", list3d, names);
    message (prog, err, ERROR);
  }

  // -- Deal with all the correlations first, before any differentiation.

  // -- Turn ABC... into Reynolds stresses.

  in['A'] -> timesMinus (*in['u'], *in['u']);
  in['B'] -> timesMinus (*in['u'], *in['v']);
  in['C'] -> timesMinus (*in['v'], *in['v']);
  if (nvel == 3) {
    in['D'] -> timesMinus (*in['u'], *in['w']);
    in['E'] -> timesMinus (*in['v'], *in['w']);
    in['F'] -> timesMinus (*in['w'], *in['w']);
  }

  // -- At this stage, in['q'] still holds the total kinetic energy.

  in['r'] -> timesMinus (*in['A'], *in['u']);
  in['r'] -> timesMinus (*in['B'], *in['v']);
  in['r'] -> timesMinus (*in['q'], *in['u']);
  in['s'] -> timesMinus (*in['B'], *in['u']);
  in['s'] -> timesMinus (*in['C'], *in['v']);
  in['s'] -> timesMinus (*in['q'], *in['v']);
  if (nvel == 3) {
    in['r'] -> timesMinus (*in['D'], *in['w']);
    in['s'] -> timesMinus (*in['E'], *in['w']);
    in['t'] -> timesMinus (*in['D'], *in['u']);
    in['t'] -> timesMinus (*in['E'], *in['v']);
    in['t'] -> timesMinus (*in['F'], *in['w']);
    in['t'] -> timesMinus (*in['q'], *in['w']);
  }

  // -- Rework in['q'] so it holds the fluctuating KE.

  *in['q'] = *in['A']; *in['q'] += *in['C'];
  if (nvel == 3) *in['q'] += *in['F'];
  *in['q'] *= 0.5;

  // -- The pressure-velocity covariances: m, n (o).

  in['m'] -> timesMinus (*in['p'], *in['u']);
  in['n'] -> timesMinus (*in['p'], *in['v']);
  if (nvel == 3) in['o'] -> timesMinus (*in['p'], *in['w']);

  // -- Everything from now on involves derivatives in one way or another.

  // -- So let's get started by making term '1' (advection):

  *work[0]  = *in['q']; work[0]   -> gradient (0);
  *work[0] *= *in['u']; *out['1']  = *work[0];
  *work[0]  = *in['q']; work[0]   -> gradient (1);
  *work[0] *= *in['v']; *out['1'] += *work[0];
  if (nvel == 3) {
   (*work[0]  = *in['q']).transform(FORWARD).gradient(2).transform(INVERSE);
    *work[0] *= *in['w']; *out['1'] += *work[0];
  }

  // -- Term '2' (turbulent transport). Destroy 'r', 's', 't' along the way:
  
  in['r'] -> gradient (0); *out['2']  = *in['r'];
  in['s'] -> gradient (1); *out['2'] += *in['s'];
  if (nvel == 3) {
    (in['t']  -> transform (FORWARD)) . gradient (2) . transform (INVERSE);
    *out['2'] += *in['t'];
  }

  // -- Term '3' (pressure work), destroying 'm', 'n', 'o' as we go:

  in['m'] -> gradient (0); *out['3']  = *in['m'];
  in['n'] -> gradient (1); *out['3'] += *in['n'];
  if (nvel == 3) {
    (in['o'] ->  transform (FORWARD)) . gradient (2) . transform (INVERSE);
    *out['3'] += *in['o'];
  }

  // -- Next we build the mean rate-of-strain tensor for further work.
  //    The '2d' components get stored in 'm', 'n', 'r' while 
  //    the '3d' ones, if needed,   go in 'o', 's', 't'.

  (*in['m']  = *in['u']) . gradient (0);
  (*in['r']  = *in['v']) . gradient (1);
  (*in['n']  = *in['u']) . gradient (1);
  (*work[0]  = *in['v']) . gradient (0);
  (*in['n'] += *work[0]) *= 0.5;

  if (nvel == 3) {
    (*in['o']  = *in['u']).transform(FORWARD).gradient(2).transform(INVERSE);
    (*work[0]  = *in['w']).gradient(0);
    (*in['o'] += *work[0]) *= 0.5;
    (*in['s']  = *in['v']).transform(FORWARD).gradient(2).transform(INVERSE);
    (*work[0]  = *in['w']).gradient(1);
    (*in['s'] += *work[0]) *= 0.5;
    (*in['t']  = *in['w']).transform(FORWARD).gradient(2).transform(INVERSE);
  }

  // -- Reduce the rate-of-strain correlations to covariances.

  *in['G'] -= *in['m'];
  *in['H'] -= *in['n'];
  *in['I'] -= *in['r'];

  if (nvel == 3) {
    *in['J'] -= *in['o'];
    *in['K'] -= *in['s'];
    *in['L'] -= *in['t'];
  }

  // -- Compute term '4' (viscous transport):

  work[0] -> times (*in['m'], *in['u']);
  work[1] -> times (*in['n'], *in['v']);
  *work[0] += *work[1];
  if (nvel == 3) {
    work[1] -> times (*in['o'], *in['w']);
    *work[0] += *work[1];
  }
  (*in['a'] -= *work[0]) . gradient (0);
  *out['4']  = *in['a'];

  work[0] -> times (*in['n'], *in['u']);
  work[1] -> times (*in['r'], *in['v']);
  *work[0] += *work[1];
  if (nvel == 3) {
    work[1] -> times (*in['s'], *in['w']);
    *work[0] += *work[1];
  }
  (*in['b'] -= *work[0]) . gradient (1);
  *out['4'] += *in['b'];

  if (nvel == 3) {
    work[0] -> times (*in['o'], *in['u']);
    work[1] -> times (*in['s'], *in['v']);
    *work[0] += *work[1];
    work[1] -> times (*in['t'], *in['w']);
    *work[0] += *work[1];
    (*in['c'] -= *work[0]).transform(FORWARD).gradient(2).transform(INVERSE);
    *out['4'] += *in['c'];
  }
  
  *out['4'] *= -2.0 * kinvis;
  
  // -- Compute term '7' (dissipation):

  *in['n'] *= sqrt (2.0);
  if (nvel == 3) {
    *in['o'] *= sqrt (2.0);
    *in['s'] *= sqrt (2.0);
  }

  *out['7'] = *in['d'];
  out['7'] -> timesMinus (*in['m'], *in['m']);
  out['7'] -> timesMinus (*in['n'], *in['n']);
  out['7'] -> timesMinus (*in['r'], *in['r']);
  if (nvel == 3) {
    out['7'] -> timesMinus (*in['o'], *in['o']);
    out['7'] -> timesMinus (*in['s'], *in['s']);
    out['7'] -> timesMinus (*in['t'], *in['t']);
  }
  *out['7'] *= 2.0 * kinvis;

  // -- Compute term '0' (production):

  *in['n'] *= sqrt (2.0);
  if (nvel == 3) {
    *in['r'] *= sqrt (2.0);
    *in['s'] *= sqrt (2.0);
  }
 
  out['0'] -> times     (*in['A'], *in['m']);
  out['0'] -> timesPlus (*in['B'], *in['n']);
  out['0'] -> timesPlus (*in['C'], *in['r']);
  if (nvel == 3) {
    out['0'] -> timesPlus (*in['D'], *in['o']);
    out['0'] -> timesPlus (*in['E'], *in['s']);
    out['0'] -> timesPlus (*in['F'], *in['t']);
  }

  // -- Finally, compute the sum, 'S' (should converge to zero):

 (*out['S']  = *out['1']);
  *out['S'] += *out['2'];
  *out['S'] += *out['3'];
  *out['S'] += *out['4'];
  *out['S'] += *out['7'];
  *out['S'] += *out['0'];
}
