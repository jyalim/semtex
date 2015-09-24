///////////////////////////////////////////////////////////////////////////////
// analysis.C: implement Analyser class for NS-type solvers.
//
// Copyright (c) 1994 <--> $Date: 2015/04/20 11:14:17 $, Hugh Blackburn
//
// This deals with output of runtime information such as step numbers,
// CFL estimation, modal energies, etc. If set, also output history
// point and particle track information, and adminster update of
// standard and phase-averaged field statistics.
//
// It is assumed that the first 2 or 3 (for 3D) entries in the Domain
// u vector are velocity fields.
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

static char RCS[] = "$Id: analysis.cpp,v 8.1 2015/04/20 11:14:17 hmb Exp $";

#include <sem.h>


Analyser::Analyser (Domain* D   ,
		    FEML*   file) :
// ---------------------------------------------------------------------------
// Set up.
//
// Fluid particle files (session.par) are dealt with here.
// Each line is of form
//
// #     tag  time  ctime  x     y      z
//       1    0.0   21.0   1.0   10.0   0.5.
//
// Output is of the same form, called session.trk. Time is the time at
// which the information was dumped, ctime the time at which the
// particle was created.
//
// NB: Particle tracking does not work for multiprocessor runs.
//
// History points are also set up here.  They are nominated in the
// optional HISTORY section of the session file.  Output is to
// session.his.
// ---------------------------------------------------------------------------
  _src (D)
{
  const char   routine[] = "Analyser::Analyser";
  char         str[StrMax];
  const real_t x_shft = Femlib::value ("X_SHIFT");
  const real_t y_shft = Femlib::value ("Y_SHIFT");
  const real_t x_scal = Femlib::value ("X_SCALE");
  const real_t y_scal = Femlib::value ("Y_SCALE");

  cout << setprecision (6);

  // -- Set up for particle tracking.

  ROOTONLY {

    // -- Open particle track file, create particles.

    ifstream pfile (strcat (strcpy (str, _src -> name), ".par"));

    if (!pfile.fail()) {
      const int_t    add = Femlib::ivalue ("SPAWN");
      int_t          id, i = 0;
      Point          P, *I;
      FluidParticle* F;

      _par_strm.open (strcat (strcpy (str, _src -> name), ".trk"));
      _par_strm.setf (ios::scientific, ios::floatfield);
      _par_strm.precision (6);

      while (pfile >> id >> P.x >> P.x >> P.x >> P.y >> P.z) {
        P.x = (P.x + x_shft) * x_scal;  // -- shft and scal default to 0 and 1
        P.y = (P.y + y_shft) * y_scal;
	F = new FluidParticle (_src, ++i, P);
	if (!(F -> inMesh())) {
	  sprintf (str, "Particle at (%f, %f, %f) not in mesh", P.x, P.y, P.z);
	  message (routine, str, WARNING);
	} else
	  _particle.push_back (F);
	if (add && F -> inMesh()) {
	  I = new Point;
	  I -> x = P.x; I -> y = P.y; I -> z = P.z;
	  _initial.insert (_initial.end(), I);
	}
      }
    }
  }

  // -- Set up for history points: open files, create points.

  if (file -> seek ("HISTORY")) {
    int_t          i, id, num = 0;
    const int_t    NH = file -> attribute ("HISTORY", "NUMBER");
    const Element* E;
    HistoryPoint*  H;
    real_t         r, s, x, y, z;

    for (i = 0; i < NH; i++) {
      file -> stream() >> id >> x >> y >> z;
      x = (x + x_shft) * x_scal;  // -- shft and scal default to 0 and 1
      y = (y + y_shft) * y_scal;
      if (E = HistoryPoint::locate (x, y, D -> elmt, r, s)) {
	H = new HistoryPoint (id, E, r, s, x, y, z);
	_history.insert (_history.end(), H);
      } else {
	sprintf (str, "History point at (%f, %f, %f) not in mesh", x, y, z);
	message (routine, str, WARNING);
      }
    }

    ROOTONLY {
      _his_strm.open (strcat (strcpy (str, _src -> name), ".his"));
      _his_strm.setf (ios::scientific, ios::floatfield);
      _his_strm.precision (6);
      if (!_his_strm) message (routine, "can't open history file", ERROR);
    }
  }

  // -- Initialise standard averaging.

  if (Femlib::ivalue ("AVERAGE")) {
    char     filename[StrMax];
    ifstream file (strcat (strcpy (filename, _src -> name), ".avg"));
    (_stats = new Statistics (D)) -> initialise (filename);

  } else
    _stats = 0;

  // -- Initialise phase averaging by setting N_PHASE > 0, as well.

  if (Femlib::ivalue ("N_PHASE") > 0) {

    if (!Femlib::ivalue ("AVERAGE"))
      message (routine, "if N_PHASE is set, AVERAGE > 0 also required", ERROR);

    // -- Must also have defined tokens STEPS_P (steps per period)
    //    and N_PHASE (number of phase points per period)
    //    and STEPS_P modulo N_PHASE must be 0
    //    and N_STEP  modulo N_PHASE must be 0
    //    and IO_FLD = STEPS_P / N_PHASE.

    if (!Femlib::ivalue ("STEPS_P"))
      message (routine, "phase averaging is on but STEPS_P not set", ERROR);

    if ( Femlib::ivalue ("STEPS_P") % Femlib::ivalue ("N_PHASE") )
      message (routine, "STEPS_P / N_PHASE non-integer", ERROR);

    if ( Femlib::ivalue ("N_STEP")  % Femlib::ivalue ("N_PHASE") )
      message (routine, "N_STEP / N_PHASE non-integer", ERROR);

    if (Femlib::ivalue("IO_FLD") != Femlib::ivalue("STEPS_P / N_PHASE"))
      message (routine, "phase averaging: IO_FLD != STEPS_P / N_PHASE", ERROR);

    _ph_stats = new Statistics (D);

  } else
    _ph_stats = 0;

  // -- Set up for output of modal energies if toggled.

  if (Femlib::ivalue ("IO_MDL")) {
    strcat (strcpy (str, _src -> name), ".mdl");
    ROOTONLY {
      _mdl_strm.open (str, ios::out);
      _mdl_strm << "#     Time Mode         Energy" << endl
		<< "# ----------------------------" << endl;
    }
  }
}


void Analyser::analyse (AuxField** work0,
			AuxField** work1)
// ---------------------------------------------------------------------------
// Step-by-step processing.  If SPAWN was set, add more particles at
// original absolute positions.
// ---------------------------------------------------------------------------
{
  const int_t verbose = Femlib::ivalue ("VERBOSE");
  const int_t cflstep = Femlib::ivalue ("IO_CFL" );
  const bool  add     = Femlib::ivalue ("SPAWN"  ) &&
    ! (_src -> step   % Femlib::ivalue ("SPAWN"  ));

  list<FluidParticle*>::iterator p;

  // -- Step-by-step updates.

  ROOTONLY {

    // -- Run information update.

    cout << "Step: " << _src -> step << "  Time: " << _src -> time << endl;

    // -- Track particles.

    if (add) {
      FluidParticle *F;
      Point          P, *I;
      int_t          i = 0, j;
      const int_t    N = _initial.size();

      for (j = 0; j < N; j++) {
	I = _initial[j];
	P.x = I -> x; P.y = I -> y; P.z = I -> z;
	_particle.push_back (F = new FluidParticle (_src, ++i, P));
      }
    }

    for (p = _particle.begin(); p != _particle.end(); p++) (*p) -> integrate();
  }

  // -- CFL, divergence information.

  if (cflstep && !(_src -> step % cflstep))
    ROOTONLY { this -> estimateCFL(); this -> divergence (work0); }

  // -- Phase averaging.

  if (_ph_stats) {
    const int_t nPeriod = Femlib::ivalue ("STEPS_P");
    const int_t nPhase  = Femlib::ivalue ("STEPS_P / N_PHASE");
    const bool  update  = !(_src -> step % nPhase);
    const int_t iPhase  =  (_src -> step % nPeriod) / nPhase;

    if (update) _ph_stats -> phaseUpdate (iPhase, work0, work1);
  }

  // -- Periodic dumps and global information.

  // -- Note that history points and particle positions are guaranteed
  //    to be available whenever a field dump is written as well as every
  //    IO_HIS steps. But note you may thereby end up with more history
  //    data than expected.

  const bool periodic = !(_src -> step %  Femlib::ivalue ("IO_HIS")) ||
                        !(_src -> step %  Femlib::ivalue ("IO_FLD")) ;
  const bool final    =   _src -> step == Femlib::ivalue ("N_STEP");
  const bool state    = periodic || final;

  if (state) {

    ROOTONLY {			// -- Output particle locations.
      Point          P;
      FluidParticle* F;

      for (p = _particle.begin(); p != _particle.end(); p++) {
	F = *p;
	if (F -> inMesh()) {
	  P = F -> location();
	  _par_strm
	    << setw (6) << F -> ID()
	    << setw(14) << _src -> time
	    << setw(14) << F -> ctime()
	    << setw(14) << P.x
	    << setw(14) << P.y
	    << setw(14) << P.z
	    << endl;
	}
      }
    }

    // -- Output history point data.

    int_t             i, j;
    const int_t       NH = _history.size();
    const int_t       NF = _src-> u.size();
    HistoryPoint*     H;
    vector<real_t>    tmp (NF);
    vector<AuxField*> u   (NF);

    for (i = 0; i < NF; i++) u[i] = _src -> u[i];

    for (i = 0; i < NH; i++) {
      H = _history[i];

      H -> extract (u, &tmp[0]);

      ROOTONLY {
	_his_strm << setw(4) << H->ID()
		  << setprecision(8) << setw(15)
		  << _src->time
		  << setprecision(6);
	for (j = 0; j < NF-1; j++) _his_strm << setw(14) << tmp[j];
	_his_strm<< setprecision(11)<< setw(19)<< tmp[NF-1]<< setprecision(6);
	_his_strm << endl;
      }
    }
  }

  // -- Statistical analysis, updated every IO_HIS steps.

  if (_stats && !(_src -> step % Femlib::ivalue ("IO_HIS")))
    _stats -> update (work0, work1);

  // -- Modal energies, written every IO_MDL steps.

  if (!(_src -> step % Femlib::ivalue ("IO_MDL")))
    this -> modalEnergy();

  // -- Field and statistical dumps.

  _src -> dump ();

  if (_stats) {
    char filename[StrMax];
    strcat (strcpy (filename, _src -> name), ".avg");
    _stats -> dump (filename);
  }
}


void Analyser::modalEnergy ()
// ---------------------------------------------------------------------------
// Print out modal energies per unit area, output by root processor.
// ---------------------------------------------------------------------------
{
  const int_t    DIM   = Geometry::nDim();
  const int_t    N     = Geometry::nModeProc();
  const int_t    base  = Geometry::baseMode();
  const int_t    nProc = Geometry::nProc();
  register int_t i, m;
  vector<real_t> ek (N);

  for (m = 0; m < N; m++) {
    ek[m] = 0.0;
    for (i = 0; i < DIM; i++) ek[m] += _src -> u[i] -> mode_L2 (m);
  }

  if (nProc > 1) {

    ROOTONLY {
      for (m = 0; m < N; m++)
	_mdl_strm << setw(10) << _src -> time
		  << setw( 5) << m
		  << setw(15) << ek[m]
		  << endl;

      for (i = 1; i < nProc; i++) {
	Femlib::recv (&ek[0], N, i);
	for (m = 0; m < N; m++)
	  _mdl_strm << setw(10) << _src -> time
		    << setw( 5) << m + i * N
		    << setw(15) << ek[m]
		    << endl;
      }

      _mdl_strm.flush();

    } else
      Femlib::send (&ek[0], N, 0);

  } else
    for (m = 0; m < N; m++)
      _mdl_strm << setw(10) << _src -> time
		<< setw( 5) << m
		<< setw(15) << ek[m]
		<< endl;
}


void Analyser::divergence (AuxField** Us) const
// ---------------------------------------------------------------------------
// Print out the velocity field's divergence energy per unit area.  Us
// is used as work area.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Analyser::divergence";
  const Geometry::CoordSys space = Geometry::system();

  const int_t  DIM = Geometry::nDim();
  const int_t  N   = Geometry::nModeProc();
  const real_t Lz  = Femlib::value ("TWOPI / BETA");
  int_t        i, m;
  real_t       L2 = 0.0;

  if (Geometry::cylindrical()) {
    for (i = 0; i < DIM; i++) *Us[i] = *_src -> u[i];
    Us[1] -> mulY();
    for (i = 0; i < DIM; i++)  Us[i] -> gradient (i);
    Us[1] -> divY();
    if (DIM == 3) Us[2] -> divY();
  } else {
    for (i = 0; i < DIM; i++) {
      *Us[i] = *_src -> u[i];
      Us[i] -> gradient (i);
    }
  }

  for (i = 1; i < DIM; i++) *Us[0] += *Us[i];

  for (m = 0; m < N; m++) L2 += Us[0] -> mode_L2 (m);

  L2 /= Lz;

  cout << "------- Divergence Energy: " << L2 << endl;

  // -- Crash stop.

  // This looks like it should always be true, but it's false if L2 is a NaN.

  if (L2 != L2) message (routine, "forcing termination on NaN.", ERROR);
}


void Analyser::estimateCFL () const
// ---------------------------------------------------------------------------
// Estimate and print the peak CFL number, based on zero-mode velocities.
// ---------------------------------------------------------------------------
{
  const real_t CFL_max = 0.7;	// -- Approximate maximum for scheme.
  const real_t SAFETY  = 0.9;	// -- Saftey factor.
  const real_t dt      = Femlib::value ("D_T");
  real_t       CFL_dt, dt_max;
  int_t        percent;

  CFL_dt = max (_src -> u[0] -> CFL (0), _src -> u[1] -> CFL (1));
  if (_src -> nField() > 3) CFL_dt = max (CFL_dt, _src -> u[2] -> CFL (2));

  dt_max  = SAFETY * CFL_max / CFL_dt;
  percent = static_cast<int_t>(100.0 * dt / dt_max);

  cout << "-- CFL: "     << CFL_dt * dt;
  cout << ", dt (max): " << dt_max;
  cout << ", dt (set): " << dt;
  cout << " ("           << percent << "%)" << endl;
}
