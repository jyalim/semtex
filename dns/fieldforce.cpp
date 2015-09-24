#include <sem.h>
#include <fieldforce.h>
#include <feml.h>

//////////////////////////////////////////////////////////////////////////////
// This code was originally contributed by Thomas Albrecht.
//
// This implements a somewhat generic body forcing interface. Basically,
// there are two (types of) classes involved here:
//
// 1. the 'FieldForce' class, which provides an interface to nonlinear.C;
//
// 2. several forcing subclasses, which actually compute a specific
//    type of forcing. Those subclasses are derived from the base
//    class 'VirtualForce'.
//
// The FieldForce instance holds the storage for the final force that
// is eventually added to the nonlinear term. To compute the final
// force, FieldForce calls all forcing subclasses and sums up their
// contribution to the final force.
//
// There are actually two storages, one for physical, and one for
// Fourier space.  A specific forcing subclass implements whichever
// suits best.
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: fieldforce.cpp,v 8.1 2015/04/20 11:14:13 hmb Exp $";

static int_t NCOM;

FieldForce::FieldForce (Domain* D   ,
                        FEML*   file)
// ---------------------------------------------------------------------------
// This constructor deals with <FORCING> section of FEML file
// It maintains a list of forcing subclasses, initialized here.
//
// On initialisation, create concrete forcing type subclasses, which
// are derived from abstract class VirtualForce.  Subclasses then read
// their respective data from the session file.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "FieldForce::FieldForce";
  const int_t nTotP     = Geometry::nTotProc();
  const int_t nzP       = Geometry::nZProc();
  char        s[StrMax];
  const       int_t verbose = Femlib::ivalue ("VERBOSE");
  int_t       i;
  
  NCOM = D -> nField() - 1;	// -- Number of velocity components.

  _D = D;

  // -- Check for FORCE section.

  if (!file -> seek (secForce)) {
    VERBOSE cout << "FORCE section not found. Disabling forcing." << endl;
    _enabled = false;
    return;
  }
  _enabled = true;

  // -- init classes
  //    NB: Sponge must be first in list.

  _classes.push_back(new SpongeForce         (_D, file));
  _classes.push_back(new CoriolisForce       (_D, file));
  _classes.push_back(new ConstForce          (_D, file));
  _classes.push_back(new WhiteNoiseForce     (_D, file));
  _classes.push_back(new SteadyForce         (_D, file));
  _classes.push_back(new ModulatedForce      (_D, file));
  _classes.push_back(new SpatioTemporalForce (_D, file));
  _classes.push_back(new DragForce           (_D, file));
  _classes.push_back(new SFDForce            (_D, file));
}


void FieldForce::addPhysical (AuxField*         Ni ,
                              AuxField*         wrk,
                              const int         com,
                              vector<AuxField*> U  )
// ---------------------------------------------------------------------------
// When called from nonlinear() to apply forcing, each subclass'
// update method is called subsequently to sum up the force, which is
// then applied to the nonlinear term.
//
// Input value com is the velocity component index: 0 <==> u; 1 <==>
// v; 2 <==> w. Routine is called component-by-component. 
// ---------------------------------------------------------------------------
{
  const char  routine[] = "FieldForce::computePhysical";

  if (!_enabled) return;

  // -- clear workspace.
  //    FIXME: should not be neccessary since SpongeForce overwrites
  //           check here if sponge is enabled

  *wrk = 0.;

  vector<VirtualForce*>::iterator p;
  for (p = _classes.begin(); p != _classes.end(); p++)
    (*p) -> physical (wrk, com, U);

  if (Geometry::cylindrical() && (com <  2)) wrk -> mulY ();
  *Ni += *wrk;

#if 0
  if (com == NCOM - 1) dump();
#endif
}


void FieldForce::addFourier (AuxField*         Ni ,
                             AuxField*         wrk,
                             const int         com,
                             vector<AuxField*> U  )
// ---------------------------------------------------------------------------
// As above, Fourier space.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "FieldForce::addFourier";
  const int_t NCOM      = _D -> nField() - 1;	// -- Number of vel. components

  if (!_enabled) return;

  *wrk = 0.;

  vector<VirtualForce*>::iterator p;
  for (p = _classes.begin(); p != _classes.end(); p++)
    (*p) -> fourier(wrk, com, U);

  if (Geometry::cylindrical() && (com <  2)) wrk -> mulY ();
  *Ni += *wrk;

#if 0
  if (com == NCOM - 1) dump();
#endif
}


void FieldForce::writeAux (vector<AuxField *> N)
// ---------------------------------------------------------------------------
// some debug stuff.
// ---------------------------------------------------------------------------
{
  int_t       step     = _D -> step;
  const bool  periodic = !(step %  Femlib::ivalue ("IO_FLD"));
  const bool  initial  =   step == Femlib::ivalue ("IO_FLD");
  const bool  final    =   step == Femlib::ivalue ("N_STEP");
  char        s[StrMax];

  if (!(periodic || final)) return;

  for (int i = 0; i < NCOM; i++) N[i]->setName ( forcename + i );
  ofstream output;
  sprintf(s, "%s.f.%03i.chk", _D->name, _D->step);
  ROOTONLY output.open (s, ios::out);

  writeField(output, _D -> name, _D->step, _D->time, N);
  ROOTONLY output.close();
}


void FieldForce::dump()
// ----------------------------------------------------------------------------
// Dump internal physical and/or Fourier space storage to file.
// ----------------------------------------------------------------------------
{
  const char  routine[]  = "FieldForce::dump";
  int         i, j = 0;
  char        s[StrMax];
  int_t       step  = _D -> step;
  const bool  periodic = !(step %  Femlib::ivalue ("IO_FLD"));
  const bool  initial  =   step == Femlib::ivalue ("IO_FLD");
  const bool  final    =   step == Femlib::ivalue ("N_STEP");

  if (!_enabled || !(periodic || final)) return;

/* -- FIXME: broken, since we no longer have all components for forcing
             available. Fix by calling component-wise.

  vector<AuxField*> vf(NCOM);
  for (i = 0; i < NCOM; i++, j++) vf[j] = _fip[i];  // dump physical space forcing
  //for (i = 0; i < NCOM; i++, j++) vf[j] = _fif[i];  // dump Fourier space forcing

  ofstream output;
  sprintf(s, "%s.ff.%03i.chk", _D->name, _D->step);
  ROOTONLY {
    output.open (s, ios::out);
    if (!output) message (routine, "can't open force file", ERROR);
  }
  writeField(output, _D -> name, _D->step, _D->time, vf);
  ROOTONLY output.close();
*/
  return;
}


// ---------------------------------------------------------------------------
// VirtualForce is the virtual base class for different types of
// forcing subclasses, providing a uniform interface (e. g., init,
// update in physical and Fourier space).
// ---------------------------------------------------------------------------


AuxField* VirtualForce::allocAuxField (Domain *D   ,
				       char    type)
{
  const char  routine[] = "VirtualForce::allocAuxField";
  const int_t nTotP     = Geometry::nTotProc();
  const int_t nzP       = Geometry::nZProc();

  return new AuxField (new real_t [(size_t)nTotP], nzP, D->elmt, type);
}


void VirtualForce::readSteadyFromFile(char*             fname, 
				      vector<AuxField*> a    )
// ---------------------------------------------------------------------------
// Read steady forcing from file, store in a.
// ---------------------------------------------------------------------------
{
  const char routine[] = "VirtualForce::readSteadyFromFile";
  ifstream   input;

  input.open(fname);
  if (!input) {
    char s[StrMax];
    sprintf(s, "can't open '%s' file", fname);
    message (routine, s, ERROR);
  }
  readField(input, a);
  input.close();
}


// ===========================================================================
// Concrete forcing subclasses follow.
// ===========================================================================


ConstForce::ConstForce (Domain* D   ,
			FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.  A force constant in both space in time, applied in
// Fourier space.
// ---------------------------------------------------------------------------
{
  const char  routine[]          = "ConstForce::ConstForce";
  const int_t verbose            = Femlib::ivalue ("VERBOSE");
  const char  tok[3][StrMax]     = {"CONST_X", 
				    "CONST_Y", 
				    "CONST_Z"};
  const char  tok_old[3][StrMax] = {"FFX", 
				    "FFY", 
				    "FFZ"};
  _D = D;

  for (int i = 0; i < NCOM; i++) {
    _v[i] = 0.;	// -- default
    if (file -> valueFromSection (&_v[i], secForce, tok[i])) {
      VERBOSE cout << "    " << tok[i]     << " = " << _v[i] << endl;
    } 
  }
}


void ConstForce::fourier (AuxField*         ff ,
			  const int         com,
			  vector<AuxField*> U  )
// ---------------------------------------------------------------------------
// Applicator.
// ---------------------------------------------------------------------------
{
  ROOTONLY if (fabs (_v[com]) > EPSDP) ff -> addToPlane (0, _v[com]);
}


SteadyForce::SteadyForce (Domain* D   ,
			  FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.  A steady, spatially varying force, computed during
// pre-processing.  Applied in physical space.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "SteadyForce::SteadyForce";
  const int_t nTotP     = Geometry::nTotProc();
  const int_t nzP       = Geometry::nZProc();
  const int_t verbose   = Femlib::ivalue ("VERBOSE");

  const char tok[3][StrMax] = {"STEADY_X", 
			       "STEADY_Y", 
			       "STEADY_Z"};
  const char tok_file[]     =  "STEADY_FILE";
  char       fname[StrMax];
  int        i;

  VERBOSE cout << "  " << routine << endl;
  _enabled = false;
  _D = D;

  _a.resize (NCOM);

  // -- If given a filename, read steady force from file.

  if (file -> valueFromSection (fname, secForce, tok_file)) {
    _enabled = true;
    for (i = 0; i < NCOM; i++) _a[i]  = allocAuxField(D, forcename + i);
    VERBOSE cout << "    reading file " << fname << endl;
    readSteadyFromFile(fname, _a);
  }
  // -- Otherwise, try to read force components from session file.
  //    If found, allocate storage.
  else for (i = 0; i < NCOM; i++) {
    char a[StrMax];
    sprintf(a, "0");	// -- default
    if (file -> valueFromSection (a, secForce, tok[i])) {
      _enabled = true;
      *(_a[i]  = allocAuxField(D, forcename + i)) = a;
      VERBOSE cout << "    " << tok[i] << " = " << a << endl;
    }
    else _a[i] = NULL;
  }
}


void SteadyForce::physical (AuxField*         ff ,
			    const int         com,
			    vector<AuxField*> U  )
// ---------------------------------------------------------------------------
// Applicator.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "SteadyForce::physical";
  const int_t verbose   = Femlib::ivalue ("VERBOSE");

  if (_a[com]) *ff += (*_a[com]);
}


WhiteNoiseForce::WhiteNoiseForce (Domain* D   ,
				  FEML*   file)
// ---------------------------------------------------------------------------
// Constructor. WhiteNoiseForce adds white noise in given
// direction(s), to all or a given mode, every _apply_step'th step
// ---------------------------------------------------------------------------
{
  const char  routine[]         = "WhiteNoiseForce::WhiteNoiseForce";
  const char  tok[3][StrMax]    = {"WHITE_EPS_X", 
				   "WHITE_EPS_Y", 
				   "WHITE_EPS_Z"};
  const char  tok_mode[StrMax]  = "WHITE_MODE";
  const char  tok_apply[StrMax] = "WHITE_APPLY_STEP";
  const int_t verbose           = Femlib::ivalue ("VERBOSE");

  VERBOSE cout << "  " << routine << endl;

  for (int i = 0; i < NCOM; i++) {
    _eps[i] = 0.;	// -- default
    if (file -> valueFromSection (&_eps[i], secForce, tok[i]))
      VERBOSE cout << "    " << tok[i] << " = " << _eps[i] << endl;
  }
  _mode = -1;
  if (file -> valueFromSection (&_mode, secForce, tok_mode)) {
    if (_mode < PERTURB_UNSET)
      message(routine, "WHITE_MODE must be >= -1", ERROR);
    VERBOSE cout << "    " << tok_mode << " = " << _mode << endl;
  }

  _apply_step = 1;
  if ((file -> valueFromSection (&_apply_step, secForce, tok_apply)))
    VERBOSE cout <<  "  Applied every " << _apply_step << ". step." << endl;
}


void WhiteNoiseForce::fourier (AuxField*         ff ,
			       const int         com, 
			       vector<AuxField*> U  )
// ---------------------------------------------------------------------------
// Applicator.
// ---------------------------------------------------------------------------
{
  if ((fabs(_eps[com]) > EPSDP) && (_D->step % _apply_step == 0))
    ff -> perturb(_mode, _eps[com]);
}


ModulatedForce::ModulatedForce (Domain* D   ,
				FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "ModulatedForce::ModulatedForce";
  const int_t verbose   = Femlib::ivalue ("VERBOSE");

  const char tok_a[3][StrMax]     = {"MOD_A_X",    
				     "MOD_A_Y",
				     "MOD_A_Z"};
  const char tok_alpha[3][StrMax] = {"MOD_ALPHA_X",
				     "MOD_ALPHA_Y",
				     "MOD_ALPHA_Z"};
  const char tok_file[]            = "MOD_A_FILE";
  char       a[StrMax], fname[StrMax];
  int        i;

  VERBOSE cout << "  " << routine << endl;
  _enabled = false;
  _D = D;

  _a.resize (NCOM);

  bool alpha_found = false;
  for (i = 0; i < NCOM; i++) {
    // -- try and read time-varying function alpha.
    sprintf(_alpha[i], "0");
    if (file -> valueFromSection (_alpha[i], secForce, tok_alpha[i])) {
      alpha_found = true;
      VERBOSE cout << "    " << tok_alpha[i] << " = " << _alpha[i] << endl;
    }
  }

  // -- disable if no alpha found.
  if (!alpha_found) return;

  // -- spatially-varying function a.
  //    If given a filename, read from file...

  if (file -> valueFromSection (fname, secForce, tok_file)) {
    _enabled = true;
    for (i = 0; i < NCOM; i++) _a[i]  = allocAuxField(D, forcename + i);
    VERBOSE cout << "    reading file " << fname << endl;
    readSteadyFromFile(fname, _a);
  }
  // -- ... otherwise try and read force components from session file.
  //    Allocate storage only if needed.
  else for (i = 0; i < NCOM; i++) {
    sprintf(a, "0");	// -- defaults
    if (file -> valueFromSection (a, secForce, tok_a[i])) {
      _enabled = true;
      _a[i]  = allocAuxField(D, forcename + i);
      *_a[i] = a;
      VERBOSE cout << "    " << tok_a[i] << " = " << a << endl;
    }
    else _a[i] = NULL;
  }
}


void ModulatedForce::physical (AuxField*         ff,
			       const int         com, 
			       vector<AuxField*> U  )
// ---------------------------------------------------------------------------
// Applicator.
// ---------------------------------------------------------------------------
{
  const char routine[] = "ModulatedForce::physical";
  real_t     alpha     = Femlib::value (_alpha[com]);

  if (_a[com] && fabs(alpha) > EPSDP) ff -> axpy (alpha, *_a[com]);
}


SpatioTemporalForce::SpatioTemporalForce (Domain* D   ,
					  FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.  SpatioTemporalForce -- f = alpha(x, t) WARNING! We
// evaluate alpha on full field EACH TIME STEP.  This severely degrades
// performance!
// ---------------------------------------------------------------------------
{
  const char   routine[] = "SpatioTemporalForce::SpatioTemporalForce";
  const int_t verbose    = Femlib::ivalue ("VERBOSE");
  const char  tok_alpha[3][StrMax] = {"SPATIOTEMP_ALPHA_X",
				      "SPATIOTEMP_ALPHA_Y",
				      "SPATIOTEMP_ALPHA_Z"};
  char        a[StrMax], fname[StrMax];
  int         i;

  VERBOSE cout << "  " << routine << endl;
  _enabled = false;
  _D = D;

  _a.resize (NCOM);

  for (i = 0; i < NCOM; i++) {
    // -- try and read spatio-temporally-varying function alpha.
    sprintf(_alpha[i], "0");
    if (file -> valueFromSection (_alpha[i], secForce, tok_alpha[i])) {
      _enabled = true;
      VERBOSE cout << "    " << tok_alpha[i] << " = " << _alpha[i] << endl;
      _a[i]  = allocAuxField(D, forcename + i);
    }
    else _a[i] = NULL;
  }
}


void SpatioTemporalForce::physical (AuxField*         ff ,
				    const int         com,
				    vector<AuxField*> U  )
// ---------------------------------------------------------------------------
// Applicator.
// ---------------------------------------------------------------------------
{
  const char routine[] = "SpatioTemporalForce::physical";

  if (_a[com]) {
    *_a[com] = _alpha[com];
    *ff += *_a[com];
  }
}


SpongeForce::SpongeForce (Domain* D   ,
			  FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "SpongeForce::SpongeForce";
  const int_t verbose   = Femlib::ivalue ("VERBOSE");
  char        s[StrMax];
  const char  tok_ref[3][StrMax] = {"SPONGE_U", 
				    "SPONGE_V", 
				    "SPONGE_W"};
  _enabled = false;
  _D = D;
  _update = 0;

  // -- setup and read sponge mask from session file

  _mask = allocAuxField(D, forcename);
  if (!(file -> valueFromSection (s, secForce, "SPONGE_M"))) {
    VERBOSE cout <<  "  SPONGE_M not found. Disabling sponge layer." << endl;
    return;
  }

  _enabled = true;

  VERBOSE cout << "    SPONGE_M = " << s << endl;

  // -- time-depended mask?

  if ((file -> valueFromSection (&_update, secForce, "SPONGE_UPDATE"))) {
    VERBOSE cout <<  "  Updating every " << _update << ". step." << endl;
    strcpy(_mask_func, s);
  }
  else
    *_mask = s;

  // -- read reference velocity from session file

  _Uref.resize(3);
  for (int i = 0; i < NCOM; i++) {
    char s[StrMax];
    sprintf(s, "0");	// -- default
    file -> valueFromSection (s, secForce, tok_ref[i]);
    _Uref[i] = allocAuxField(_D, 'r' + i);
    *_Uref[i]  = s;

  }
}


void SpongeForce::physical (AuxField*         ff ,
			    const int         com,
			    vector<AuxField*> U  )
// ---------------------------------------------------------------------------
// Applicator.
// -- Since this overwrites ff, SpongeForce must be first in list!!
// ---------------------------------------------------------------------------
{
  const char   routine[] = "SpongeForce::physical";
  const int_t verbose    = Femlib::ivalue ("VERBOSE");

  if (!_enabled) return;

  if (_update && ((_D->step % _update) == 0)) *_mask = _mask_func;

  ff -> vvmvt (*_Uref[com], *U[com], *_mask);   // ff = (uref - u) * mask
}


DragForce::DragForce (Domain* D   ,
		      FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.
//
// Drag force: |F| = - m(X) * U/|U| * |U|^2 (i.e., F and U are coplanar)
// componentwise, this reads:
//   Fx = -m(X) |U| u
//   Fy = -m(X) |U| v
//   Fz = -m(X) |U| w
// m(X) is a shape function.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "DragForce::DragForce";
  const int_t verbose   = Femlib::ivalue ("VERBOSE");
  char        s[StrMax];

  _enabled = false;

  // -- setup and read mask from session file
  if (!(file -> valueFromSection (s, secForce, "DRAG_M"))) {
    VERBOSE cout <<  "  DRAG_M not found. Disabling drag force." << endl;
    return;
  }
  _D = D;
  _mask = allocAuxField (_D, 0);
  _umag = allocAuxField (_D, 0);

  _enabled = true;
  VERBOSE cout << "  DRAG_M = " << s << endl;
  *_mask = s;
}


void DragForce::physical (AuxField*         ff ,
			  const int         com,
			  vector<AuxField*> U  )
// ---------------------------------------------------------------------------
// Applicator.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "DragForce::physical";
  const int_t verbose   = Femlib::ivalue ("VERBOSE");

  if (!_enabled) return;

  // -- compute velocity magnitude
  if (com == 0) {
    _umag  -> mag(U);
    *_umag *= *_mask;
  }

  ff -> timesMinus (*_umag, *U[com]);
}


CoriolisForce::CoriolisForce (Domain* D   ,
			      FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.
// ---------------------------------------------------------------------------
{
  const char routine[] = "CoriolisForce::CoriolisForce";
  const int_t verbose = Femlib::ivalue ("VERBOSE");
  const char  tokOmega[3][StrMax]    = {"CORIOLIS_OMEGA_X",
					"CORIOLIS_OMEGA_Y",
					"CORIOLIS_OMEGA_Z"};
  const char  tokDomegaDt[3][StrMax] = {"CORIOLIS_DOMEGA_X_DT",
                                        "CORIOLIS_DOMEGA_Y_DT",
                                        "CORIOLIS_DOMEGA_Z_DT"};
  const char  tokUnsteady[StrMax]    = {"CORIOLIS_UNSTEADY"};
  char        s[StrMax];

  VERBOSE cout << "  " << routine << endl;

  _enabled = false;
  _unsteady = 0; // int_t, bec. there's no bool version of valueFromSection()
  _D = D;
  _o.resize(3);
  _minus_2o.resize(3);

  _minus_o.resize(3);
  _DoDt.resize(3);

  file -> valueFromSection (&_unsteady, secForce, tokUnsteady);

  // -- try and read omega's components from session file

  for (int i = 0; i < 3; i++) {

    sprintf(_omega[i], "0");	// -- default
    sprintf(_DomegaDt[i], "0");
    if (NCOM == 3 || i == 2) {
      if (file -> valueFromSection (_omega[i], secForce, tokOmega[i])) {
        _enabled = true;
        VERBOSE cout << "    " << tokOmega[i] << " = " << _omega[i] << endl;
      }

      if (_unsteady && 
	  (file -> valueFromSection (_DomegaDt[i], secForce, tokDomegaDt[i])))
        VERBOSE cout << "    " 
		     << tokDomegaDt[i] << " = " << _DomegaDt[i] << endl;
    }
  }

  if (_enabled && !_unsteady) {
    VERBOSE cout << "    Steady Coriolis," 
      " ignoring possible CORIOLIS_DOMEGA_[XYZ]_DT."     << endl;
    VERBOSE cout << "    Also, make sure" 
      " you include centrifugal force via STEADY_[XYZ]." << endl;
    // FIXME: Could we also evaluate centrifugal force during pre-processing,
    //        then add it to SteadyForce somehow?
  }

  //  if (NCOM == 2) 
  //    VERBOSE cout << "    2-D: ignoring possible CORIOLIS_[D]OMEGA_[XY][_DT]."
  //		 << endl;

  if (!_unsteady)
    for (int i = 0; i < 3; i++) _minus_2o[i] = -2. * Femlib::value (_omega[i]);

  _a.resize (NCOM);
  for (int i = 0; i < NCOM; i++) _a[i] = allocAuxField(D, forcename + i);
}


void CoriolisForce::physical (AuxField*               ff ,
			      const int               com,
			      const vector<AuxField*> U  )
// ---------------------------------------------------------------------------
// Applicator.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "CoriolisForce::physical";
  const int_t verbose   = Femlib::ivalue ("VERBOSE");

  if (!_enabled) return;

  if (NCOM == 2 && Geometry::cylindrical())
    message (routine, "2-D, cylindrical not implemented yet.", ERROR);
  if (_unsteady) {
    if (com == 0) 
      for (int i = 0; i < 3; i++) {
	if (NCOM == 2) i = 2;
	_o[i] = Femlib::value (_omega[i]);

        // -- we'll need several negative values later on.

	_minus_o[i]  = - _o[i];
	_minus_2o[i] = - 2. * _o[i];
	_DoDt[i]     = - Femlib::value (_DomegaDt[i]);
      }

    ff -> crossXPlus (com, _DoDt);            // - DOmega/Dt x X

    // - Omega x Omega x X  (actually, we compute  + Omega x ((- Omega) x X))

    if (com == 0)
      for (int i = 0; i < NCOM; i++) {
	*_a[i] = 0.;
	_a[i] -> crossXPlus (i, _minus_o);
      }
    ff -> crossProductPlus (com, _o, _a);
  }

  ff -> crossProductPlus (com, _minus_2o, U); // - 2 Omega x U
}


SFDForce::SFDForce (Domain* D   ,
		    FEML*   file)
// ---------------------------------------------------------------------------
// Constructor.
// 
// Internal storage _a plays the role of qbar in Akervik at al.'s
// description. 
// ---------------------------------------------------------------------------
{
  const char  routine[] = "SFDForce::SFDForce";
  const int_t verbose   = Femlib::ivalue ("VERBOSE");

  _enabled = false;

  if (!(file -> valueFromSection (&_SFD_DELTA, secForce, "SFD_DELTA") ||
	file -> valueFromSection (&_SFD_CHI,   secForce, "SFD_CHI"))) return;


  if (!(file -> valueFromSection (&_SFD_DELTA, secForce, "SFD_DELTA") &&
	file -> valueFromSection (&_SFD_CHI,   secForce, "SFD_CHI"))) {
    message (routine, "SFD_DELTA and SFD_CHI must be paired.", ERROR);
  }

  _D = D;
  _enabled = true;

  _a.resize (NCOM);
  for (int i = 0; i < NCOM; i++)
    *(_a[i] = allocAuxField (D, forcename + i)) = 0.0;

  VERBOSE {
    cout << "  SFD_DELTA = " << _SFD_DELTA << endl;
    cout << "  SFD_CHI   = " << _SFD_CHI   << endl;
  }
}


void SFDForce::physical (AuxField*         ff ,
			 const int         com,
			 vector<AuxField*> U  )
// ---------------------------------------------------------------------------
// Applicator.  Forcing -SFD_CHI * (u - qbar) is added to the
// (negative of the) nonlinear terms. Then qbar is updated using
// forwards Euler.
// ---------------------------------------------------------------------------
{
  const char   routine[] = "SFDForce::physical";
  const int_t  verbose   = Femlib::ivalue ("VERBOSE");
  const real_t dt        = Femlib::value  ("D_T");
  static int   step      = 0;  // -- Flag restart.

  if (!_enabled) return;

  // -- Restart qbar from previous velocity field.
  
  if (step < NCOM) *_a[com] = *U[com];

  // -- Subtract CHI*(u-ubar) from -nonlinear terms.

  ff -> axpy (-_SFD_CHI, * U[com]);
  ff -> axpy (+_SFD_CHI, *_a[com]);

  // -- Explicit update for ubar = qbar = _a.

  *_a[com] *= 1.0 - dt / _SFD_DELTA;
  _a[com]  -> axpy (dt / _SFD_DELTA, *U[com]);

  step++;
}
