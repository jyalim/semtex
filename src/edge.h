#ifndef EDGE_H
#define EDGE_H

class Edge
// ===========================================================================
// Element-edge data and operators.
// ===========================================================================
{
public:
  Edge (const char*, const Element*, const int_t); 

  int_t  bOff  () const { return _elmt -> ID() * Geometry::nExtElmt(); }
  int_t  dOff  () const { return _doffset; }
  int_t  dSkip () const { return _dskip; }

  const real_t* nx () const { return _nx; }
  const real_t* ny () const { return _ny; }

  void   get      (const real_t*,real_t*)                                const;
  void   geometry (real_t*,real_t*,real_t* = 0,real_t* = 0,real_t* = 0)  const;
  bool   isCurved () const { return _curved; }

  void   sideGrad (const real_t* u, real_t* ux, real_t* uy, real_t* wk) const 
                       { _elmt -> sideGrad (_side, u + _eoffset, ux, uy, wk); }
  void   curlCurl (const int_t,
		   const real_t*,const real_t*,const real_t*,
		   const real_t*,const real_t*,const real_t*,
		   real_t*,real_t*,real_t*,real_t*,real_t*)              const;

  void   mulY     (real_t*)                                              const;
  void   divY     (real_t*)                                              const;
  void   sideEvaluate (const char*, real_t*)                             const;

  bool   inGroup      (const char* grp) const { return !(strcmp (grp,_group)); }
  void   addForGroup  (const char*, const real_t,  real_t*)              const;
  void   setForGroup  (const char*, const real_t,  real_t*)              const;
  void   dotInForGroup(const char*, const Vector&, real_t*)              const;

  real_t vectorFlux (const char*,const real_t*,const real_t*,real_t*)    const;
  real_t scalarFlux (const char*,const real_t*,real_t*)                  const;

  Vector normTraction (const char*,const real_t*,real_t*)                const;
  Vector tangTraction (const char*,const real_t*,const real_t*,real_t*)  const;

  void traction (const int_t,const real_t,
		 const real_t*,const real_t*,const real_t*,const real_t*,
		 const real_t*,const real_t*,const real_t*,const real_t*,
		 real_t*,real_t*,real_t*,real_t*,real_t*,real_t*,real_t*)const;

protected:
  int_t          _np     ;	// Matches Geometry::nP().
  char*          _group  ;	// Group string.

  const Element* _elmt   ;	// Corresponding element.
  int_t          _side   ;	// Corresponding side.

  int_t          _eoffset;      // Offset of corresponding element in Field.
  int_t          _doffset;	// Offset in Field data plane (matches BLAS).
  int_t          _dskip  ;	// Skip   in Field data plane (matches BLAS).

  real_t*        _x      ;	// Node locations.
  real_t*        _y      ;	//
  real_t*        _nx     ;	// Unit outward normal components at nodes.
  real_t*        _ny     ;	// 
  real_t*        _area   ;	// Weighted multiplier for parametric mapping.

  bool           _curved ;      // True if the edge is curved.
};

#endif
