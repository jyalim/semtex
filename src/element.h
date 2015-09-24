#ifndef ELEMENT_H
#define ELEMENT_H

class Element
// ===========================================================================
// Virtual 2D quadrilateral element class, equal order in each direction.
//
//         s                                y ^     +4
//         ^                                  |    / \
//     4   |                                  |  1+   +3
//     +------+                               |    \  |
//     |   |  |      <== Logical  space       |     \ |
//     |   +--|-> r      Physical space ==>   |      \|
//     |      |                               |       +2
//     +------+                               +-----------> x
//     1      2
//
// Master element coordinates: -1 < r,s < 1; edge traverses CCW.
//
// Elements in this class will use Gauss--Lobatto--Legendre
// meshes, integration, and nodal basis functions.
//
// All 2D storage is row-major.
// ===========================================================================
{
public:
  Element (const int_t,const int_t,const Mesh*);
  ~Element();
  
  int_t ID () const { return _id; }

  // -- Elemental Helmholtz matrix constructor, operator.

  void HelmholtzSC   (const real_t,const real_t,real_t*,real_t*,
		      real_t*,real_t*,real_t*,int_t*)                    const;
  void HelmholtzDiag (const real_t,const real_t,real_t*,real_t*)         const;
  void HelmholtzKern (const real_t,const real_t,
		      real_t*,real_t*,real_t*,real_t*)                   const;
  void HelmholtzOp   (const real_t,const real_t,real_t*,real_t*,real_t*) const;

  // -- Local/global projectors.

  void local2global      (const real_t*,const int_t*,real_t*,real_t*)    const;
  void local2globalSum   (const real_t*,const int_t*,real_t*,real_t*)    const;
  void local2globalSumSC (real_t*,const int_t*,real_t*,const real_t*,
			  real_t*)                                       const;
  void global2local      (real_t*,const int_t*,const real_t*,
			  const real_t*)                                 const;
  void global2localSC    (const real_t*,const int_t*,real_t*,real_t*,
			  const real_t*,const real_t*,real_t*)           const;

  // -- Project from one interpolation order to another.

  void project (const int_t,const real_t*,const int_t,real_t*,real_t*)   const;
  
  // -- Element-boundary operators.

  void bndryDsSum  (const int_t*,const real_t*,real_t*)                  const;
  void bndryMask   (const int_t*,real_t*,const real_t*,const int_t*)     const;
  void bndryInsert (const int_t*,const real_t*,real_t*)                  const;

  // -- Element-side operators.

  void sideGeom  (const int_t,real_t*,real_t*,real_t*,real_t*,real_t*)   const;
  void sideEval  (const int_t,real_t*,const char*)                       const;
  void sideGrad  (const int_t,const real_t*,real_t*,real_t*,real_t*)     const;
  void sideGet   (const int_t,const real_t*,real_t*)                     const;
  void sideGetY  (const int_t,real_t*)                                   const;
  void sideMulY  (const int_t,const real_t*,real_t*)                     const;
  void sideDivY  (const int_t,const real_t*,real_t*)                     const;
  void sideDivY2 (const int_t,const real_t*,real_t*)                     const;

  // -- Element-operator functions.

  void grad  (real_t*,real_t*,real_t*)             const;
  void gradX (const real_t*,const real_t*,real_t*) const;
  void gradY (const real_t*,const real_t*,real_t*) const;

  void divY (real_t*) const;
  void mulY (real_t*) const;
  void mulX (real_t*) const;
  void crossXPlus (const int, const real_t,const vector<real_t>&,real_t*) const;
  
  void evaluate (const char*,real_t*) const;

  real_t integral (const char*,  real_t*) const;
  real_t integral (const real_t*,real_t*) const;
  real_t momentX  (const char*,  real_t*) const;
  real_t momentY  (const char*,  real_t*) const;
  real_t momentX  (const real_t*,real_t*) const;
  real_t momentY  (const real_t*,real_t*) const;
  real_t area     ()                      const;
  void   weight   (real_t*)               const;

  void   lengthScale (real_t*)                                           const;
  real_t CFL         (const real_t,const real_t*,const real_t*,real_t*)  const;

  real_t norm_inf (const real_t*) const;
  real_t norm_L2  (const real_t*) const;
  real_t norm_H1  (const real_t*) const;

  // -- Probe functions.

  bool   locate (const real_t,const real_t,real_t&,real_t&,
		 real_t*,const bool = false)                      const;
  real_t probe  (const real_t,const real_t,const real_t*,real_t*) const;

  // -- Debugging/informational routines.

  void printMesh  ()              const;
  void printBndry (const real_t*) const;
  void printMatSC (const real_t*,const real_t*,const real_t*)            const;
  void Helmholtz  (const real_t,const real_t,real_t*,real_t*,real_t*)    const;

protected:

  const int_t   _id   ;		// Element identifier.
  const int_t   _np   ;		// Number of points on an edge.
  const int_t   _npnp ;		// Total number = np * np.
  const int_t   _next ;		// Number of points on periphery.
  const int_t   _nint ;		// Number of internal points.
  const bool    _cyl  ;         // Cylindrical coordinate problem.

  const real_t* _zr   ;		// Master element mesh points on [-1, 1], r.
  const real_t* _zs   ;		// Master element mesh points on [-1, 1], s.
  const real_t* _wr   ;		// Master element quadrature weights, r.
  const real_t* _ws   ;		// Master element quadrature weights, s.
  const real_t* _DVr  ;		// Master element derivative operator, r.
  const real_t* _DTr  ;		// Transpose.
  const real_t* _DVs  ;		// Master element derivative operator, s.
  const real_t* _DTs  ;		// Transpose.
  const real_t* _SDVr ;		// As above, but with SVV modifications.
  const real_t* _SDTr ;
  const real_t* _SDVs ;
  const real_t* _SDTs ;

  int_t*        _emap ;		// Indices of edges in nodal matrices.
  int_t*        _pmap ;		// Inversion of emap (pmap[emap[i]] = i).

  real_t*       _xmesh;		// Physical space mesh.
  real_t*       _ymesh;		// 2D row-major store.

  real_t*       _drdx ;		// Partial derivatives (r, s) --> (x, y),
  real_t*       _dsdx ;		//   evaluated at quadrature points.
  real_t*       _drdy ;		//   (2D row-major storage.)
  real_t*       _dsdy ;		//

  real_t*       _delta;		// Local length scale.

  real_t*       _Q1   ;		// Geometric factor 1 at quadrature points.
  real_t*       _Q2   ;		// Geometric factor 2 at quadrature points.
  real_t*       _Q3   ;		// Geometric factor 3 at quadrature points.
  real_t*       _Q4   ;		// Geometric factor 4 at quadrature points.
  real_t*       _Q8   ;	        // Like _Q4 but without possible factor of y.

  // -- Geometric and quadrature-specific internal functions.

  void mapping      ();
  void HelmholtzRow (const real_t,const real_t,const int_t,
		     const int_t,real_t*,real_t*) const;

  // -- BLAS-conforming edge offsets & skips for element-edge traverses.

  void terminal (const int_t side,int_t& start,int_t& skip) const
  { switch (side) {
  case 0: start = 0;             skip  = 1;    break;
  case 1: start = _np - 1;       skip  = _np;  break;
  case 2: start = _np*(_np - 1); skip  = -1;   break;
  case 3: start = 0;             skip  = -_np; break;
  } }
};

#endif
