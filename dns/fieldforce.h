#ifndef FIELDFORCE_H
#define FIELDFORCE_H

const char secForce[] = "FORCE";
const char forcename  = 'u';		// forcing fields are uvw


class VirtualForce
// ---------------------------------------------------------------------------
// Virtual base class for various types of forcing.
// ---------------------------------------------------------------------------
{
public:
  void      allocStorage       (Domain*);
  AuxField* allocAuxField      (Domain*, char);
  void      readSteadyFromFile (char*, vector<AuxField*>);

  virtual void physical (AuxField*, const int, vector<AuxField*>) {};
  virtual void fourier  (AuxField*, const int, vector<AuxField*>) {};

  vector<AuxField*> _a;		// -- storage for pre-processed part

protected:
  Domain* _D;
  bool    _enabled;
};


class FieldForce
// ---------------------------------------------------------------------------
// Provides external access for applications. See e.g. calls in nonlinear.C.
// ---------------------------------------------------------------------------
{
public:
  FieldForce            (Domain*, FEML*);
  void addPhysical      (AuxField*, AuxField*, const int, vector<AuxField*>);
  void addFourier       (AuxField*, AuxField*, const int, vector<AuxField*>);
  void dump             ();
  void writeAux		(vector<AuxField*>);
protected:
  bool			_enabled;
  vector<VirtualForce*> _classes;   // -- vector of concrete forcing classes
  Domain*		_D;
  vector<AuxField*>	_u;         // -- storage for physical space velocity
};


class ConstForce : public VirtualForce
// ---------------------------------------------------------------------------
// A force constant in both space in time, applied in Fourier space.
// ---------------------------------------------------------------------------
{
public:
  ConstForce            (Domain*, FEML*);
  void fourier          (AuxField*, const int, vector<AuxField*>);
protected:
  real_t                _v[3];	// Force components
};


class SteadyForce : public VirtualForce
// ---------------------------------------------------------------------------
// Constant in time but a function of space, applied in physical space.
// ---------------------------------------------------------------------------
{
public:
  SteadyForce           (Domain*, FEML*);
  void physical         (AuxField*, const int, vector<AuxField*>);
protected:
};


class WhiteNoiseForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Forcing stochastic in space and time, applied in Fourier space.
// ---------------------------------------------------------------------------
{
public:
  WhiteNoiseForce       (Domain*, FEML*);
  void fourier		(AuxField*, const int, vector<AuxField*>);
protected:
  real_t                _eps[3];
  int_t                 _mode;
  int_t                 _apply_step; // apply force every _apply_step'th step.
};


class ModulatedForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Forcing which is a constant function of space (may be read from
// file) modulated by a function of time.  Applied in physical space.
// ---------------------------------------------------------------------------
{
public:
  ModulatedForce        (Domain*, FEML*);
  void physical         (AuxField*, const int, vector<AuxField*>);
protected:
  char                  _alpha[3][StrMax]; // -- temporally varying part.
};


class SpatioTemporalForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Forcing that is an arbitrary function of space-time. Physical space.
// ---------------------------------------------------------------------------
{
public:
  SpatioTemporalForce   (Domain*, FEML*);
  void physical         (AuxField*, const int, vector<AuxField*>);
protected:
  char                  _alpha[3][StrMax]; // -- spatio-temporal-varying part.
};


class SpongeForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Forcing penalises difference between velocity field and a given
// function of spatial position. Physical space.
// ---------------------------------------------------------------------------
{
public:
  SpongeForce           (Domain*, FEML*);
  void physical         (AuxField*, const int, vector<AuxField*>);
protected:
  vector<AuxField*>     _Uref;
  AuxField*             _mask;
  char                  _mask_func[StrMax]; // -- mask function, f(x,y,z,t)
  int                   _update; // mask update frequency
};


class DragForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Forcing acts against velocity field according to its magnitude. Physical.
// ---------------------------------------------------------------------------
{
public:
  DragForce             (Domain*, FEML*);
  void physical         (AuxField*, const int, vector<AuxField*>);
protected:
  AuxField		*_mask;
  AuxField		*_umag;
};


class CoriolisForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Forcing appopriate to solution in a rotating frame of reference.
// ---------------------------------------------------------------------------
{
public:
  CoriolisForce         (Domain*, FEML*);
  void physical         (AuxField*, const int, vector<AuxField*>);
  void OmegaTimesOmegaTimesX();
protected:
  char                  _omega[3][StrMax];    // -- angular velocity = f(t) ..
  char                  _DomegaDt[3][StrMax]; // -- and its time derivative
  vector<real_t>        _o;		      // -- evaluated at current time
  vector<real_t>        _minus_o;	      // -- - omega
  vector<real_t>        _minus_2o;	      // -- - 2 * omega
  vector<real_t>        _DoDt;		      // -- evaluated at current time
  int_t                 _unsteady;            // -- 1 if omega is unsteady
};


class SFDForce : virtual public VirtualForce
// ---------------------------------------------------------------------------
// Selective frequency damping: add forcing designed to achieve steady
// solution to NSE.  Physical space.  Parameters SFD_DELTA and SFD_CHI.
//
// Reference: Akervik et al., (2006), Steady solutions to the
// Navier--Stokes equations by selective frequency damping, Phys
// Fluids 18: 068102.
// ---------------------------------------------------------------------------
{
public:
  SFDForce              (Domain*, FEML*);
  void physical         (AuxField*, const int, vector<AuxField*>);
protected:
  real_t                _SFD_DELTA, _SFD_CHI;
};


#endif
