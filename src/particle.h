#ifndef PARTICLE_H
#define PARTICLE_H

class FluidParticle
// ===========================================================================
// Class used to locate and integrate positions of massless particles.
// ===========================================================================
{
public:
  FluidParticle (Domain*, const int_t, Point&);
  void            integrate (); 
  int_t           ID        () const { return _id;     }
  real_t          ctime     () const { return _ctime;  }
  const  Element* inMesh    () const { return _E;      }
  const  Point&   location  () const { return _p;      } 
  static int_t    IDMax     ()       { return _ID_MAX; }

private:
  const int_t    _id   ;	// Numeric tag.
  const real_t   _ctime;	// Time of initialisation.
  int_t          _step ;	// Number of integration steps.
  const Element* _E    ;        // Pointer to the element particle is in.
  real_t         _r    ;        // Corresponding "r" location within element.
  real_t         _s    ;	// likewise for "s".
  Point          _p    ;	// Physical space location.
  real_t*        _u    ;	// Multilevel "x" velocity storage.
  real_t*        _v    ;	// Multilevel "y" velocity storage.
  real_t*        _w    ;	// Multilevel "z" velocity storage.

  static Domain* _Dom    ;	// Velocity fields and class functions.
  static int_t   _NCOM   ;	// Number of velocity components.
  static int_t   _NEL    ;	// Number of elements in mesh.
  static int_t   _NZ     ;	// Number of z planes.
  static int_t   _TORD   ;	// Order of N--S timestepping.
  static int_t   _ID_MAX ;	// Highest issued id.
  static real_t* _P_coeff;	// Integration (predictor) coefficients.
  static real_t* _C_coeff;	// Integration (corrector) coefficients.
  static real_t* _Work   ;	// Work area for point location routines.
  static real_t  _DT     ;	// Time step.
  static real_t  _Lz     ;	// Periodic length.
};

#endif
