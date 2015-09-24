///////////////////////////////////////////////////////////////////////////////
// data2df.C: simple 2D x Fourier data class; an AuxField without
// geometric information.  Also included in this file are Header class
// definitions for IO & maintenance of input file information.
//
// Copyright (c) 2004 <--> $Date: 2015/04/20 11:14:17 $, Hugh Blackburn
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

static char RCS[] = "$Id: data2df.cpp,v 8.1 2015/04/20 11:14:17 hmb Exp $";

#include <sem.h>
#include <data2df.h>


Data2DF::Data2DF (const int_t nP  ,
		  const int_t nZ  ,
		  const int_t nEl ,
		  const char  Name) :
  _name (Name),
  _np   (nP  ),
  _nz   (nZ  ),
  _nel  (nEl ),
  _np2  (nP * nP)
// ---------------------------------------------------------------------------
// Data2DF constructor. Storage area set to zero.
// ---------------------------------------------------------------------------
{
  int_t i;
  
  _nplane = _np * _np * _nel;
  if (_nplane & 1) _nplane++;
  _ntot   = _nplane * _nz;

  _data  = new real_t  [_ntot];
  _plane = new real_t* [_nz];

  for (i = 0; i < _nz; i++) _plane[i] = _data + i*_nplane;
  Veclib::zero (_ntot, _data, 1);
}


Data2DF& Data2DF::DFT1D (const int_t sign)
// ---------------------------------------------------------------------------
// Carry out discrete Fourier transformation in z direction.
// ---------------------------------------------------------------------------
{
  if (_nz > 2) Femlib::DFTr (_data, _nz, _nplane, sign);

  return *this;
}


Data2DF& Data2DF::DPT2D (const int_t sign, 
			 const char  basis)
// ---------------------------------------------------------------------------
// Carry out 2D discrete polynomial transform (element-by-element) on planes.
//
// Use basis = 'l' to transform to Legendre polynomial basis, 'm' for
// modal polynomial basis.
// ---------------------------------------------------------------------------
{
  int_t          i;
  vector<real_t> work (_nplane);
  const real_t   *Fu, *Ft, *Bu, *Bt;

  if (basis == 'l')
    Femlib::legTran (_np, &Fu, &Ft, &Bu, &Bt, 0, 0);
  else
    Femlib::modTran (_np, &Fu, &Ft, &Bu, &Bt, 0, 0);

  if (sign == FORWARD)
    for (i = 0; i < _nz; i++)
      Femlib::tpr2d (_plane[i], _plane[i], &work[0], Fu, Ft, _np, _np, _nel);
  else
    for (i = 0; i < _nz; i++)
      Femlib::tpr2d (_plane[i], _plane[i], &work[0], Bu, Bt, _np, _np, _nel);

  return *this;
}


Data2DF& Data2DF::F_conjugate (const bool zero)
// ---------------------------------------------------------------------------
// Take complex conjugate in the Fourier coordinate direction. If zero
// is true, assume that mode zero is complex instead of two real_t
// modes packed together.
// ---------------------------------------------------------------------------
{
  int_t       i;
  const int_t first = (zero) ? 1 : 3;

  if (_nz > 1)
    for (i = first; i < _nz; i += 2)
      Veclib::neg (_nplane, _plane[i], 1);

  return *this;
} 


Data2DF& Data2DF::F_symmetrize (const bool zero)
// ---------------------------------------------------------------------------
// Symmetrize velocity field component in Fourier coordinate
// direction. This enforces a reflection symmetry of the velocity and
// pressure fields, as follows.
// 
// 'u': mode_k.Im = 0, k > 0
// 'v': mode_k.Im = 0, k > 0
// 'w': mode_k.Re = 0, k > 0
// 'p': mode_k.Im = 0, k > 0
//
// If zero is true, assume that mode zero is complex instead of two
// real modes packed together.
// ---------------------------------------------------------------------------
{
  int_t i, first;

  switch (_name) {
  case 'u': case 'v': case 'p': first = (zero) ? 1 : 3; break;
  case 'w':                     first = (zero) ? 0 : 2; break;
  default:                      return *this;           break;
  }

  if (_nz > 1)
    for (i = first; i < _nz; i += 2)
      Veclib::zero (_nplane, _plane[i], 1);

  return *this;
} 


Data2DF& Data2DF::F_shift (const real_t alpha,
			   const bool   zero )
// ---------------------------------------------------------------------------
// Use the shift-rotation duality of the Fourier transform to shift
// the data a proportion alpha of the fundamental length in the
// Fourier coordinate direction.  Data are assumed to be in
// Fourier-transformed state on input.
// ---------------------------------------------------------------------------
{
  const int_t    N = _nz >> 1;
  const int_t    first = (zero) ? 0 : 1;
  register int_t i;
  int_t          k;
  real_t         cosA, sinA, tmp;
  real_t         *Re, *Im;

  for (k = first; k < N; k++) {
    Re   = _plane[2*k];
    Im   = _plane[2*k + 1];
    cosA = cos(TWOPI * (k-first+1) * alpha);
    sinA = sin(TWOPI * (k-first+1) * alpha);
    for (i = 0; i < _nplane; i++) {
      tmp = Re[i];
      Re[i] = Re[i] * cosA - Im[i] * sinA;
      Im[i] = tmp   * sinA + Im[i] * cosA;
    }
  }

  return *this;
}


Data2DF& Data2DF::operator = (const Data2DF& rhs)
// ---------------------------------------------------------------------------
// If the two fields conform, copy rhs's data storage to lhs.
//
// Otherwise perform projection/interpolation of rhs's data area to lhs.
// Interpolation ASSUMES THAT FOURIER TRANSFORMATION HAS ALREADY OCCURRED
// in z direction if rhs is 3D.  Truncation of Fourier modes occurs if this
// Data2DF has less modes than rhs (to avoid aliasing).
// ---------------------------------------------------------------------------
{
  if (rhs._nel != _nel)
    message ("Data2DF::operator =", "fields can't conform", ERROR);

  if (rhs._np == _np && rhs._nz == _nz)
    Veclib::copy (_ntot, rhs._data, 1, _data, 1);

  else {			// -- Perform projection.

    register int_t  i, k;
    register real_t *LHS, *RHS;
    const real_t    *IN,  *IT;
    const int_t     nzm = min (rhs._nz, _nz);
    vector<real_t>  work (rhs._np * _np);
    real_t*         tmp = &work[0];

    Femlib::projection (&IN, &IT, rhs._np, GLJ, 0.0, 0.0, _np, GLJ, 0.0, 0.0);

    for (k = 0; k < nzm; k++) {	// -- 2D planar projections.
      LHS = _plane[k];
      RHS = rhs._plane[k];

      if (rhs._np == _np)
	Veclib::copy (_nplane, RHS, 1, LHS, 1);
      else
	for (i = 0; i < _nel; i++, LHS += _np2, RHS += rhs._np2) {
	  Blas::mxm (IN,  _np, RHS, rhs._np, tmp, rhs._np);
	  Blas::mxm (tmp, _np, IT,  rhs._np, LHS,     _np);
	}
    }

    if ((i = _nz - rhs._nz) > 0) // -- Zero pad for Fourier projections.
      Veclib::zero (i * _nplane, _data + rhs._ntot, 1);
  }

  return *this;
}


Data2DF& Data2DF::operator += (const Data2DF& rhs)
// ---------------------------------------------------------------------------
// If two fields conform, add rhs's data storage pointwise into lhs.
// ---------------------------------------------------------------------------
{
  if (rhs._nel != _nel || rhs._np != _np || rhs._nz != _nz)
    message ("Data2DF::operator +=", "fields don't conform", ERROR);

  Veclib::vadd (_ntot, rhs._data, 1, _data, 1, _data, 1);

  return *this;
}


Data2DF& Data2DF::operator -= (const Data2DF& rhs)
// ---------------------------------------------------------------------------
// If two fields conform, subtarct rhs's data storage pointwise from lhs.
// ---------------------------------------------------------------------------
{
  if (rhs._nel != _nel || rhs._np != _np || rhs._nz != _nz)
    message ("Data2DF::operator +=", "fields don't conform", ERROR);

  Veclib::vsub (_ntot, rhs._data, 1, _data, 1, _data, 1);

  return *this;
}


Data2DF& Data2DF::operator *= (const Data2DF& rhs)
// ---------------------------------------------------------------------------
// If two fields conform, multiply rhs's data storage pointwise into lhs.
// ---------------------------------------------------------------------------
{
  if (rhs._nel != _nel || rhs._np != _np || rhs._nz != _nz)
    message ("Data2DF::operator *=", "fields don't conform", ERROR);

  Veclib::vmul (_ntot, rhs._data, 1, _data, 1, _data, 1);

  return *this;
}


Data2DF& Data2DF::operator = (const real_t val)
// ---------------------------------------------------------------------------
// Set field storage area to val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (_ntot,      _data, 1);
  else              Veclib::fill (_ntot, val, _data, 1);

  return *this;
}


Data2DF& Data2DF::operator *= (const real_t val)
// ---------------------------------------------------------------------------
// Multiply field storage area by val.
// ---------------------------------------------------------------------------
{
  Blas::scal (_ntot, val, _data, 1);

  return *this;
}


Data2DF& Data2DF::reverse ()
// ---------------------------------------------------------------------------
// Reverse order of bytes within each word of data.
// ---------------------------------------------------------------------------
{
  Veclib::brev (_ntot, _data, 1, _data, 1);

  return *this;
}


ostream& operator << (ostream&  strm,
		      Data2DF& F    )
// ---------------------------------------------------------------------------
// Binary write of F's data area.
// ---------------------------------------------------------------------------
{
  int_t i;
  
  for (i = 0; i < F._nz; i++)
    strm.write (reinterpret_cast<char*> (F._plane[i]),
		F._np * F._np * F._nel * sizeof (real_t));

  return strm;
}


istream& operator >> (istream&  strm,
		      Data2DF& F    )
// ---------------------------------------------------------------------------
// Binary read of F's data area.
// ---------------------------------------------------------------------------
{
  int_t i;
  
  for (i = 0; i < F._nz; i++)
    strm.read (reinterpret_cast<char*> (F._plane[i]),
	       F._np * F._np * F._nel * sizeof (real_t));

  return strm;
}


Data2DF& Data2DF::filter1D (const real_t roll ,
			    const int_t  order)
// ---------------------------------------------------------------------------
// Carry out filtering in the Fourier domain, with a Boyd--VanDeven
// filter (see Leven, Iskandarani and Haidvogel JCP 137). The input
// parameter roll [0, 1] specifies where the filter begins to roll
// over, while parameter order, 2 or larger, specifies the shape of the
// filter.
//
// It is assumed that the data are already Fourier transformed on input.
// ---------------------------------------------------------------------------
{
  int_t          i;
  const int_t    nh = _nz >> 1, ord = max (2, order);
  vector<real_t> filter (nh + 1), mask (_nz);
  const real_t   lag = clamp (roll, 0.0, 1.0);

  Femlib::erfcFilter (nh, ord, lag * nh, 1.0, &filter[0]);
  mask[0] = filter[ 0];
  mask[1] = filter[nh];
  for (i = 1; i < nh; i++) mask[2*i] = mask[2*i+1] = filter[i];

  for (i = 0; i < _nz; i++)
    Veclib::smul  (_nplane, mask[i], _plane[i], 1, _plane[i], 1);

  return *this;
}


Data2DF& Data2DF::filter2D (const real_t roll ,
			    const int_t  order)
// ---------------------------------------------------------------------------
// Carry out 2D filtering element-by-element in the modal polynomial
// domain, with a Boyd--VanDeven filter (see Leven, Iskandarani and
// Haidvogel JCP 137). The input parameter roll [0, 1] specifies where
// the filter begins to roll over, while parameter order, 2 or larger,
// specifies the shape of the filter.
//
// Doing the filtering in the modal polynomial space ensures that the field
// retains C0 continuity across element boundaries, see Blackburn & Schmidt
// JCP 186.
// ---------------------------------------------------------------------------
{
  int_t          i, j;
  const real_t   lag = clamp (roll, 0.0, 1.0);
  const int_t    ord = max (2, order);
  vector<real_t> filter (_np);
  vector<real_t> work (_nplane);
  vector<real_t> Iu (_np2), It (_np2);
  const real_t   *Du, *Dt, *dpt;
  real_t         *dataplane;

  Femlib::erfcFilter (_np - 1, ord, lag, 1.0, &filter[0]);

  // -- Polynomial transform+filter (Du, Dt) and inversion (Iu, It) matrices.

  Femlib::modTran (_np, &Du, &Dt, 0, &dpt, 0, 0);

  for (i = 0; i < _np; i++)
    Veclib::smul (_np, filter[i], dpt + i*_np, 1, &It[i*_np], 1);

  Blas::mxm    (Dt, _np, &It[0], _np, &Iu[0], _np);
  Veclib::copy (_np2,    &Iu[0], 1,   &It[0], 1);

  for (i = 0; i < _np; i++)
    for (j = 0; j < _np; j++)
      Iu[Veclib::row_major (j, i, _np)] = It[Veclib::row_major (i, j, _np)];

  for (i = 0; i < _nz; i++)
    Femlib::tpr2d (_plane[i],_plane[i],&work[0],&Iu[0],&It[0],_np,_np,_nel);

  return *this;
}


Data2DF& Data2DF::reflect2D (vector<int_t>& pos,
			     vector<int_t>& neg)
// ---------------------------------------------------------------------------
// Use gathr_scatr maps to carry out a 2D spatial reflection of data.
// Note that no sign change occurs (data is scalar) -- if needed that
// must be done separately.
// ---------------------------------------------------------------------------
{
  int_t          k;
  vector<real_t> tmp (_nplane);

  for (k = 0; k < _nz; k++) {
    Veclib::copy (_nplane, _plane[k], 1, &tmp[0], 1);
    Veclib::gathr_scatr (pos.size(), &tmp[0], &neg[0], &pos[0], _plane[k]);
  }

  return *this;
}


//////////////////////////////////////////////////////////////////////////////
// Header class definitions are also in this file.
//////////////////////////////////////////////////////////////////////////////


Header::Header ()
// ---------------------------------------------------------------------------
// Allocate storage for Header.
// ---------------------------------------------------------------------------
{
  sess = new char [StrMax];
  sesd = new char [StrMax];
  flds = new char [StrMax];
  frmt = new char [StrMax];

  sess[0] = sesd[0] = flds[0] = '\0';
  sprintf (frmt, "binary "); Veclib::describeFormat (frmt + strlen (frmt));
  nr = ns = nz = nel = step = 0;
  time = dt = visc = beta = 0.0;
}


istream& operator >> (istream& file,
		      Header&  hdr )
// ---------------------------------------------------------------------------
// Insert data into Header struct from file.
// ---------------------------------------------------------------------------
{
  char routine[] = "operator: istream >> Header";
  char s[StrMax];

  if (file.get(hdr.sess, 25).eof()) return file; file.getline(s, StrMax);
  file.get(hdr.sesd, 25);                        file.getline(s, StrMax);
  file >> hdr.nr >> hdr.ns >> hdr.nz >> hdr.nel; file.getline(s, StrMax);
  file >> hdr.step;                              file.getline(s, StrMax);
  file >> hdr.time;                              file.getline(s, StrMax);
  file >> hdr.dt;                                file.getline(s, StrMax);
  file >> hdr.visc;                              file.getline(s, StrMax);
  file >> hdr.beta;                              file.getline(s, StrMax);
  file >> hdr.flds;                              file.getline(s, StrMax);
  file.get(hdr.frmt, 26);                        file.getline(s, StrMax);

  if (!file) message (routine, "failed reading header information", ERROR);
  return file;
}


ostream& operator << (ostream& file,
		      Header&  hdr )
// ---------------------------------------------------------------------------
// Put data from Header struct onto file. Use current time info.  If
// the data are in binary format, over-ride the stated input format
// with whatever the current machine supports.
// ---------------------------------------------------------------------------
{
  const char routine [] = "operator: ofstream << Header";
  const char *hdr_fmt[] = { 
    "%-25s "                  "Session\n",
    "%-25s "                  "Created\n",
    "%-5d %-5d %-5d %-5d   "  "Nr, Ns, Nz, Elements\n",
    "%-25d "                  "Step\n",
    "%-25.6g "                "Time\n",
    "%-25.6g "                "Time step\n",
    "%-25.6g "                "Kinvis\n",
    "%-25.6g "                "Beta\n",
    "%-25s "                  "Fields written\n",
    "%-25s "                  "Format\n"
  };

  char   s1[StrMax], s2[StrMax];
  time_t tp (time (0));

  strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));

  sprintf  (s1, hdr_fmt[0], hdr.sess);                        file << s1;
  sprintf  (s1, hdr_fmt[1], s2);                              file << s1;
  sprintf  (s1, hdr_fmt[2], hdr.nr, hdr.ns, hdr.nz, hdr.nel); file << s1;
  sprintf  (s1, hdr_fmt[3], hdr.step);                        file << s1;
  sprintf  (s1, hdr_fmt[4], Femlib::value ("t"));             file << s1;
  sprintf  (s1, hdr_fmt[5], Femlib::value ("D_T"));           file << s1;
  sprintf  (s1, hdr_fmt[6], Femlib::value ("KINVIS"));        file << s1;
  sprintf  (s1, hdr_fmt[7], Femlib::value ("BETA"));          file << s1;
  sprintf  (s1, hdr_fmt[8], hdr.flds);                        file << s1;

  sprintf  (s2, "%s", "binary "); Veclib::describeFormat (s2 + strlen (s2));

  sprintf  (s1, hdr_fmt[9], s2);                              file << s1;

  if (!file) message (routine, "failed writing field file header", ERROR);
  file << flush;

  return file;
}


bool Header::swab() const
// ---------------------------------------------------------------------------
// Return true if coding of binary information in *this conflicts with 
// that for the machine (indicating byte swapping is required).
// ---------------------------------------------------------------------------
{
  char routine[] = "Header::swab";
  char machine[StrMax];
  bool swap = false;

  Veclib::describeFormat (machine);

  if (!strstr (frmt, "binary"))
    message (routine, "input field file not in binary format", ERROR);
  
  if (!strstr (frmt, "endian"))
    message (routine, "input field file in unknown binary format", WARNING);
  else
    swap = ((strstr (machine, "big") && strstr (frmt,    "little")) ||
	    (strstr (frmt,    "big") && strstr (machine, "little")) );
  
  return swap;
}
