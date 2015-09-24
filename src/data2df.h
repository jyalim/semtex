#ifndef DATA2DF_H
#define DATA2DF_H

///////////////////////////////////////////////////////////////////////////////
// Define simple routines to handle quad spectral element x Fourier
// data (class Data2DF), plus simple header data I/O class (class
// Header) both with public data.
///////////////////////////////////////////////////////////////////////////////


class Data2DF
// ============================================================================
// Canonical field class, each np X np element is defined on [-1,1] X [-1, 1].
// Data are arranged element-ordered in 2D planes to create a 3D scalar field.
// ============================================================================
{
friend istream& operator >> (istream&, Data2DF&);
friend ostream& operator << (ostream&, Data2DF&);

public:
  Data2DF  (const int_t nP, const int_t nZ, const int_t nEl,
	    const char Name='\0');
  ~Data2DF () { delete [] _data; delete [] _plane; }

  char getName () { return _name; }
  Data2DF& reverse ();

  Data2DF& operator  =   (const real_t);
  Data2DF& operator *=   (const real_t);
  Data2DF& operator  =   (const Data2DF&);
  Data2DF& operator +=   (const Data2DF&);
  Data2DF& operator -=   (const Data2DF&);
  Data2DF& operator *=   (const Data2DF&);

  Data2DF& DFT1D         (const int_t);
  Data2DF& DPT2D         (const int_t, const char);

  Data2DF& filter1D      (const real_t, const int_t);
  Data2DF& filter2D      (const real_t, const int_t);

  Data2DF& F_conjugate   (const bool);
  Data2DF& F_symmetrize  (const bool);
  Data2DF& F_shift       (const real_t alpha, const bool);

  Data2DF& reflect2D     (vector<int_t>&, vector<int_t>&);

#if 0 // -- Sorry: wanted the convenience of public data.
protected:
#endif
  const char  _name;
  const int_t _np, _nz, _nel, _np2;
  int_t       _nplane, _ntot;
  real_t*     _data;
  real_t**    _plane;
};


class Header
// ===========================================================================
// Nekton/Prism/Semtex-compatible header struct + I/O routines +
// public data descriptors.  No array data storage.
// ===========================================================================
{
public:
  Header();
 ~Header() { delete [] sess; delete [] sesd; delete [] flds; delete [] frmt; }

  bool  swab    () const;
  int_t nFields () const { return strlen (flds); }

  char*  sess;
  char*  sesd;
  int_t  nr  ;
  int_t  ns  ;
  int_t  nz  ;
  int_t  nel ;
  int_t  step;
  real_t time;
  real_t dt  ;
  real_t visc;
  real_t beta;
  char*  flds;
  char*  frmt;
};
istream& operator >> (istream&, Header&);
ostream& operator << (ostream&, Header&);

#endif
