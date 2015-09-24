#ifndef CONDITION_H
#define CONDITION_H

class Condition
// ===========================================================================
// Virtual base class for boundary condition application.
//
// Each concrete class is derived from the virtual base class Condition:
// 1. Essential         : essential BC with constant, supplied, value.
// 2. EssentialFunction : essential BC, value obtained by parsing a function.
// 3. Natural           : natural BC with constant, supplied, value.
// 4. NaturalFunction   : natural BC, value obtained by parsing a function.
// 5. Mixed             : transfer coefficient type, 2 supplied values.
// 6. NaturalCBCp       : "high-order" pressure BC, natural, computed value.
// 7. EssentialCBCp     : Outflow-specific pressure BC, essential, computed.
// 8. NaturalCBCu       : Outflow-specific u velocity BC, natural, computed.
// 9. NaturalCBCv       : Outflow-specific v velocity BC, natural, computed.
//
// Note that for supplied-value BC types, the value is a physical-space
// description of the BC, and is now set once at run-time (cannot be reset).
// This is not true for those obtained by parsing a function: this can
// be re-parsed every timestep (if so flagged by dns command-line option
// -t or -tt).
//
// Also, each condition class derived from the base has to define all the
// pure virtual functions listed below (except the destructor) but
// some of these will just be stubs that do nothing for any particular
// type.  Those stubs are indicated in the present header
// file with the function body "{ };".
//
// For essential/Dirichlet BCs, the method "set" must be defined;
// For natural/Neumann BCs, the method "sum" must be defined, while
// For mixed/Robin BCs, the three methods "augmentSC", "augmentOp", 
// and "augmentDG" must be defined.
//
// See also boundary.h, edge.h, bcmgr.cpp mesh.cpp.
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
// ===========================================================================
{
public:
  virtual void evaluate  (const Field*   src    ,
			  const int_t    id     , 
			  const int_t    plane  , 
			  const Element* elmt   ,
			  const int_t    side   , 
			  const int_t    step   ,
			  const bool     Fourier,
			  real_t*        tgt    )                      const=0;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                      const=0;
  virtual void sum       (const int_t, const int_t*,
		          const real_t*,const real_t*,real_t*,real_t*) const=0;
  virtual void augmentSC (const int_t,  const int_t, const int_t,
			  const int_t*,const real_t*,real_t*,real_t*)  const=0;
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)       const=0;
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                      const=0;
  virtual void describe  (char* tgt)                                   const=0;

  virtual ~Condition()   { }
};


class Essential : public Condition
// ===========================================================================
// Essential BC applicator.  This one is for plain (constant value)
// Dirichlet/essential BCs.
// ===========================================================================
{
public:
  Essential              (const char* v) : _value (strtod (v, 0)) { }
  virtual void evaluate  (const Field*, const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const;
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const
    { };
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void describe  (char*)                                         const;
private:
  real_t _value;
};


class EssentialFunction : public Condition
// ===========================================================================
// Essential BC applicator for specified function Dirichlet/essential BCs.
// This uses parser to set the BC values. 
// ===========================================================================
{
public:
  EssentialFunction      (const char*);
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const;
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const
    { };
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void describe  (char*)                                         const;
private:
  char* _function;
};


class Natural : public Condition
// ===========================================================================
// Natural BC applicator.  This one is for plain (constant value)
// Neumann/natural BCs.
// ===========================================================================
{
public:
  Natural                (const char* v) : _value (strtod (v, 0)) { }
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void describe  (char*)                                         const;
private:
  real_t _value;
};


class NaturalFunction : public Condition
// ===========================================================================
// Natural BC applicator for specified function Neumann/natural BCs.
// This uses parser to set the BC values. 
// ===========================================================================
{
public:
  NaturalFunction        (const char*);
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void describe  (char*)                                         const;
private:
  char* _function;
};


class Mixed : public Condition
// ===========================================================================
// Boundary condition class for mixed (a.k.a. Robin) type BCs of form
//     dc/dn + K(c - C) = 0.
// Mixed BCs affect problem Helmholtz matrices, but only on the diagonal, 
// and element-boundary, terms. Syntax in session file is
//     <M> c = K, C </M>  or 
//     <M> c = K; C </M> 
// where 'c' is a field name and K and C can be evaluated as constants
// (perhaps using defined TOKENS). White space following the separators
// above (',', ';') preceding C is optional.
// ===========================================================================
{
public:
  Mixed                  (const char*);
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const
    { };
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*, const real_t*,
			  const real_t*, real_t*, real_t*)               const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const;
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const;
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const;
  virtual void describe  (char*)                                         const;
private:
  real_t _K_;		// -- This is "K" above.
  real_t _C_;		// -- This is "C" above.
};


// -- Classes with internally computed BC ("CBC") types follow.
//    Evaluate (and other) functions dealt with in bcmgr.cpp.


class NaturalCBCp : public Condition
// ===========================================================================
// Computed Neumann BC for pressure, typical of wall boundaries.
//
// Karniadakis, Israeli & Orszag JCP 97, (1991).
// ===========================================================================
{
public:
  NaturalCBCp            (BCmgr* B) : _BCmgr (B) { }
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void describe  (char*)                                         const;
private:
  BCmgr* _BCmgr; 
};


class EssentialCBCp : public Condition
// ===========================================================================
// Computed Dirichlet BC for pressure on outflow boundaries.
//
// Dong, Karniadakis & Chryssostomides, JCP 261 (2014, DKC14).
// ===========================================================================
{
public:
  EssentialCBCp          (BCmgr* B) : _BCmgr (B) { }
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const;
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const
    { };
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void describe  (char*)                                         const;
private:
  BCmgr* _BCmgr; 
};


class NaturalCBCu : public Condition
// ===========================================================================
// Computed Neumann BC for velocity component 'u' on outflow boundaries.
//
// Dong, Karniadakis & Chryssostomides, JCP 261 (2014, DCK14).
// ===========================================================================
{
public:
  NaturalCBCu            (BCmgr* B) : _BCmgr (B) { }
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void describe  (char*)                                         const;
private:
  BCmgr* _BCmgr; 
};


class NaturalCBCv : public Condition
// ===========================================================================
// Computed Neumann BC for velocity component 'v' on outflow boundaries.
// Evaluated in Fourier space.
//
// Dong, Karniadakis & Chryssostomides, JCP 261 (2014, DKC14).
// ===========================================================================
{
public:
  NaturalCBCv            (BCmgr* B) : _BCmgr (B) { }
  virtual void evaluate  (const Field*,   const int_t, const int_t,
                          const Element*, const int_t, const int_t,
			  const bool, real_t*)                           const;
  virtual void set       (const int_t, const int_t*,
			  const real_t*, real_t*)                        const
    { };
  virtual void sum       (const int_t, const int_t*,
			  const real_t*,const real_t*,real_t*,real_t*)   const;
  virtual void augmentSC (const int_t, const int_t, const int_t,
			  const int_t*, const real_t*, real_t*, real_t*) const
    { };
  virtual void augmentDg (const int_t, const int_t*, 
			  const real_t*, real_t*)                        const
    { };
  virtual void augmentOp (const int_t, const int_t*,
			  const real_t*, const real_t*, real_t*)         const
    { };
  virtual void describe  (char*)                                         const;
private:
  BCmgr* _BCmgr; 
};

#endif
