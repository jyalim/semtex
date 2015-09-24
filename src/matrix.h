#ifndef MATRIX_H
#define MATRIX_H


class MatrixSys
// ===========================================================================
// System of global and local Helmholtz matrices and partitions.
// Matrix factorisations are Cholesky, use LAPACK storage schemes.
// ===========================================================================
{
friend class Field;
//friend ostream& operator << (ostream&, MatrixSys&);
//friend istream& operator >> (istream&, MatrixSys&);
public:
  MatrixSys     (const real_t, const real_t, const int_t,
		 const vector<Element*>&, const BoundarySys*,
		 const SolverKind);
 ~MatrixSys     ();
  bool match    (const real_t, const real_t, const NumberSys*,
		 const SolverKind) const;

private:
  real_t  _HelmholtzConstant;	// Same for all modes.
  real_t  _FourierConstant  ;	// Varies with mode number.
 
  const vector<Boundary*>& _BC;	// Internal copy of Boundary conditions.
  const NumberSys*         _NS;	// Internal copy of NumberSys.

  int_t    _nel     ;		// Number of elemental matrices.
  int_t    _nglobal ;		// Number of unique element-boundary nodes.
  int_t    _singular;		// If system is potentially singular.
  int_t    _nsolve  ;		// System-specific number of global unknowns.
  SolverKind _method  ;		// Flag specifies direct or iterative solver.

  // -- For _method == DIRECT:

  int_t    _nband ;		// Bandwidth of global matrix (incl. diagonal).
  int_t    _npack ;		// Number of real_ts for global matrix.
  real_t*  _H     ;		// (Factored) packed global Helmholtz matrix.
  real_t** _hbi   ;		// Element external-internal coupling matrices.
  real_t** _hii   ;		// (Factored) internal-internal matrices.
  int_t*   _bipack;		// Size of hbi for each element.
  int_t*   _iipack;		// Size of hii for each element.

  // -- For _method == JACPCG:

  int_t    _npts;		// Total number of unique meshpoints.
  real_t*  _PC  ;		// Diagonal preconditioner matrix.
};

ostream& operator << (ostream&, MatrixSys&);
istream& operator >> (istream&, MatrixSys&);


class ModalMatrixSys
// ===========================================================================
// A way of organising MatrixSys*'s by Fourier mode.
// ===========================================================================
{
public:
  ModalMatrixSys (const real_t, const real_t, const int_t, const int_t,
		  const vector<Element*>&, const BoundarySys*, 
		  const SolverKind);
 ~ModalMatrixSys ();

  const MatrixSys* operator [] (const int_t i) const { return _Msys[i]; }

private:
  char*              _fields;	// Character field tags for this system.
  vector<MatrixSys*> _Msys  ;	// One MatrixSys for each Fourier mode.
};

#endif
