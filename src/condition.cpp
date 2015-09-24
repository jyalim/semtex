//////////////////////////////////////////////////////////////////////////////
// condition.cpp: functions used to evaluate & apply BCs.
//
// Copyright (c) 1994 <--> $Date: 2015/04/20 11:14:17 $, Hugh Blackburn
//
// All classes inherit the (semi-abstract) base class Condition.
// Owing to the different behaviour of essential and natural BCs, a
// number of routines take no action; due also to the generality of
// the base class function calls, many routines do not use all
// parameters.
//
// The two basic kinds of conditions are essential (Dirichlet) and
// natural (Neumann).  Essential conditions are set/imposed directly,
// while natural conditions are applied by summation of integral
// terms, owing to the fact that they derive from integration by
// parts.  So function "set" is only used by essential-type conditions
// while function "sum" is only used by natural-type conditions.
//
// The HOPBC condition are a special case of natural conditions used
// for the pressure Poisson equation: they are derived from the
// momentum equations and computed using an extrapolative process
// described in [1].
//
// A further general class of available condition is the Robin or
// mixed type, of form dc/dn + K( c - C ) = 0. See [2]. 

// NB: if you invent a new BC which is of essential/Dirichlet type,
// you also need to edit mesh.cpp/buildMask() so that enumerate sets the
// global numbering mask to 1 for this type. Otherwise the BC won't be
// correctly applied.  In this, Robin/mixed count as natural/Neumann.
//
// REFERENCES
// ----------
// [1] Karniadakis, Israeli & Orszag (1991) "High-order splitting methods
//     for the incompressible Navier--Stokes equations", JCP 91(2).
// [2] Blackburn (2001) "Dispersion and diffusion in coated tubes of 
//     arbitrary cross-section", Comput Chem Eng 25(2/3).
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

static char RCS[] = "$Id: condition.cpp,v 8.1 2015/04/20 11:14:17 hmb Exp $";

#include <sem.h>


void Essential::evaluate (const Field*   src    ,
			  const int_t    id     , 
			  const int_t    plane  , 
			  const Element* elmt   ,
			  const int_t    side   , 
			  const int_t    step   ,
			  const bool     Fourier,
			  real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external value storage area tgt with installed constant.
// Physical space.
// ---------------------------------------------------------------------------
{
  if (!Fourier) Veclib::fill (Geometry::nP(), _value, tgt, 1);
}


void Essential::set (const int_t   side,
		     const int_t*  bmap,
		     const real_t* src ,
		     real_t*       tgt ) const
// ---------------------------------------------------------------------------
// Scatter external value storage area src into globally-numbered tgt. 
// ---------------------------------------------------------------------------
{
  const int_t  nm    = Geometry::nP() - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::scatr (nm, src, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] = src[nm];
  else             tgt[start[nm]] = src[nm];  
}


void Essential::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "essential:\t%g", _value);
}


EssentialFunction::EssentialFunction (const char* f)
// ---------------------------------------------------------------------------
// Essential condition that sets value by interpreting function
// string, which can be an explicit function of x, y, z & t as well as
// any installed token values.  Otherwise the same as Essential class.
// ---------------------------------------------------------------------------
{
  strcpy ((_function = new char [strlen (f) + 1]), f);
}


void EssentialFunction::evaluate (const Field*   src    ,
				  const int_t    id     , 
				  const int_t    plane  , 
				  const Element* elmt   ,
				  const int_t    side   , 
				  const int_t    step   ,
				  const bool     Fourier,
				  real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external tgt by interpreting function.
// Physical space.
// ---------------------------------------------------------------------------
{
  if (!Fourier) elmt -> sideEval (side, tgt, _function);
}


void EssentialFunction::set (const int_t   side,
			     const int_t*  bmap,
			     const real_t* src ,
			     real_t*       tgt ) const
// ---------------------------------------------------------------------------
// Scatter external value storage area src into globally-numbered tgt
// (as for Essential class).
// ---------------------------------------------------------------------------
{
  const int_t  nm    = Geometry::nP() - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::scatr (nm, src, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] = src[nm];
  else             tgt[start[nm]] = src[nm];  
}


void EssentialFunction::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "essential-function:\t%s", _function);
}


void Natural::evaluate (const Field*   src    ,
			const int_t    id     , 
			const int_t    plane  , 
			const Element* elmt   ,
			const int_t    side   , 
			const int_t    step   ,
			const bool     Fourier,
			real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external value storage area tgt with installed constant.
// Physical space.
// ---------------------------------------------------------------------------
{
  if (!Fourier) Veclib::fill (Geometry::nP(), _value, tgt, 1);
}


void Natural::sum (const int_t   side  ,
		   const int_t*  bmap  ,
		   const real_t* src   ,
		   const real_t* weight,
		   real_t*       work  ,
		   real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms (from src) into globally-numbered tgt.
// Work vector must be np long.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void Natural::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "natural:\t%g", _value);
}


NaturalFunction::NaturalFunction (const char* f)
// ---------------------------------------------------------------------------
// Natural condition that sets value by interpreting function string,
// which can be an explicit function of x, y, z & t as well as any
// installed token values.  Otherwise the same as Natural class.
// ---------------------------------------------------------------------------
{
  strcpy ((_function = new char [strlen (f) + 1]), f);
}


void NaturalFunction::evaluate (const Field*   src    ,
				const int_t    id     , 
				const int_t    plane  , 
				const Element* elmt   ,
				const int_t    side   , 
				const int_t    step   ,
				const bool     Fourier,
				real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external tgt by interpreting function.
// Physical space.
// ---------------------------------------------------------------------------
{
  if (!Fourier) elmt -> sideEval (side, tgt, _function);
}


void NaturalFunction::sum (const int_t   side  ,
			   const int_t*  bmap  ,
			   const real_t* src   ,
			   const real_t* weight,
			   real_t*       work  ,
			   real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void NaturalFunction::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "natural-function:\t%s", _function);
}


Mixed::Mixed (const char* v)
// ---------------------------------------------------------------------------
// The format for a Mixed BC is: "field = mulvalue;refvalue", where in
// dc/dn + K( c - C ) = 0, K is mulvalue and C is refvalue.  The
// separator must be ';' or ',' (and without white space). Each of the
// two supplied values is expected to evaluate to a real_t constant.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Mixed::Mixed";
  const char sep[] = ";,";
  char       buf[StrMax], *tok;

  strcpy (buf, v);
  _K_ = Femlib::value (tok = strtok (buf,  sep));
  _C_ = Femlib::value (tok = strtok (NULL, sep));
}


void Mixed::sum (const int_t   side  ,
		 const int_t*  bmap  ,
		 const real_t* src   ,
		 const real_t* weight,
		 real_t*       work  ,
		 real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.  This is
// used to add K*C to RHS forcing.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::smul (np, _K_ * _C_, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];
}


void Mixed::augmentSC (const int_t   side  ,
		       const int_t   nband ,
		       const int_t   nsolve,
		       const int_t*  bmap  ,
		       const real_t* area  ,
		       real_t*       work  ,
		       real_t*       tgt   )  const
// ---------------------------------------------------------------------------
// Add <K, w> terms to (banded LAPACK) matrix tgt.
// ---------------------------------------------------------------------------
{
  const int_t    np    = Geometry::nP();
  const int_t    nm    = np - 1;
  const int_t*   start = bmap;
  register int_t i, k;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::smul (np, _K_, area, 1, work, 1);

  for (i = 0; i < nm; i++)
    if ((k = start[i]) < nsolve)
      tgt[Lapack::band_addr (k, k, nband)] += work[i];

  i = (side == 3) ? bmap[0] : start[nm];
  if (i < nsolve) tgt[Lapack::band_addr (i, i, nband)] += work[nm];
}


void Mixed::augmentOp (const int_t   side, 
		       const int_t*  bmap,
		       const real_t* area,
		       const real_t* src ,
		       real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are mixed BCs.  Add in diagonal terms
// <K*src, w> to tgt.  Both src and tgt are globally-numbered vectors.
// ---------------------------------------------------------------------------
{
  const int_t    np    = Geometry::nP();
  const int_t    nm    = np - 1;
  const int_t*   start = bmap;
  register int_t i;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += _K_ * area[i] * src[start[i]];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += _K_ * area[nm] * src[i];
}


void Mixed::augmentDg (const int_t   side, 
		       const int_t*  bmap,
		       const real_t* area,
		       real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  const int_t    np    = Geometry::nP();
  const int_t    nm    = np - 1;
  const int_t*   start = bmap;
  register int_t i;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += _K_ * area[i];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += _K_ * area[nm];
}


void Mixed::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "mixed:\t%g\t%g", _K_, _C_);
}

//////////////////////////////////////////////////////////////////////////////
// Internally computed BC types follow.
//
// These are all evaluated in Fourier space.
//////////////////////////////////////////////////////////////////////////////

void NaturalCBCp::evaluate (const Field*   src    ,
			    const int_t    id     , 
			    const int_t    plane  , 
			    const Element* elmt   ,
			    const int_t    side   , 
			    const int_t    step   ,
			    const bool     Fourier,
			    real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external via a call to CBCmgr to compute terms.
// ---------------------------------------------------------------------------
{
  if (Fourier) _BCmgr -> evaluateCNBCp (id, plane, step, tgt); 
}


void NaturalCBCp::sum (const int_t   side  ,
		       const int_t*  bmap  ,
		       const real_t* src   ,
		       const real_t* weight,
		       real_t*       work  ,
		       real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void NaturalCBCp::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "computed-natural-pressure, see KIO91");
}


void EssentialCBCp::evaluate (const Field*   src    ,
			      const int_t    id     , 
			      const int_t    plane  , 
			      const Element* elmt   ,
			      const int_t    side   , 
			      const int_t    step   ,
			      const bool     Fourier,
			      real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external via a call to CBCmgr to compute terms.
// ---------------------------------------------------------------------------
{
  if (Fourier) _BCmgr -> evaluateCEBCp (src, id, plane, step, tgt); 
}


void EssentialCBCp::set (const int_t   side,
			 const int_t*  bmap,
			 const real_t* src ,
			 real_t*       tgt ) const
// ---------------------------------------------------------------------------
// Scatter external value storage area src into globally-numbered tgt. 
// ---------------------------------------------------------------------------
{
  const int_t  nm    = Geometry::nP() - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::scatr (nm, src, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] = src[nm];
  else             tgt[start[nm]] = src[nm];  
}


void EssentialCBCp::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "computed-essential-pressure, see DKC14");
}


void NaturalCBCu::evaluate (const Field*   src    ,
			    const int_t    id     , 
			    const int_t    plane  , 
			    const Element* elmt   ,
			    const int_t    side   , 
			    const int_t    step   ,
			    const bool     Fourier,
			    real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external via a call to CBCmgr to compute terms.
// ---------------------------------------------------------------------------
{
  if (Fourier) _BCmgr -> evaluateCNBCu (src, id, plane, step, 'u', tgt); 
}


void NaturalCBCu::sum (const int_t   side  ,
		       const int_t*  bmap  ,
		       const real_t* src   ,
		       const real_t* weight,
		       real_t*       work  ,
		       real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void NaturalCBCu::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "computed-natural-velocity (u component) see DKC14 (8c)");
}


void NaturalCBCv::evaluate (const Field*   src    ,
			    const int_t    id     , 
			    const int_t    plane  , 
			    const Element* elmt   ,
			    const int_t    side   , 
			    const int_t    step   ,
			    const bool     Fourier,
			    real_t*        tgt    ) const
// ---------------------------------------------------------------------------
// Load external via a call to CBCmgr to compute terms.
// ---------------------------------------------------------------------------
{
  if (Fourier) _BCmgr -> evaluateCNBCu (src, id, plane, step, 'v', tgt); 
}


void NaturalCBCv::sum (const int_t   side  ,
		       const int_t*  bmap  ,
		       const real_t* src   ,
		       const real_t* weight,
		       real_t*       work  ,
		       real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void NaturalCBCv::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "computed-natural-velocity (v component) see DKC14 (8c)");
}
