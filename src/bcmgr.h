#ifndef BCMGR_H
#define BCMGR_H

typedef struct bctriple { char group; int_t elmt; int_t side; } BCtriple;

class BCmgr
// ===========================================================================
// This is a factory / retrieval service for classes derived from
// Condition, and maintains GROUP descriptors.  In addition, it reads
// and returns NumberSys objects from session.num.  Also, it contains
// code for maintenance and evaluation of computed BC types.
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
  BCmgr (FEML*, vector<Element*>&);

  const char*        field        () const { return _fields; }
  const char*        groupInfo    (const char) const;
  Condition*         getCondition (const char, const char, const int_t = 0);
  NumberSys*         getNumberSys (const char, const int_t = 0);
  vector<BCtriple*>& getBCedges   () { return _elmtbc; }
  int_t              nBCedges     () const { return _elmtbc.size(); }
  int_t              nWall        (); // Should be const: OSX compiler bug?

  class CondRecd {
  public: 
    char       grp  ;
    char       fld  ;
    Condition* bcn  ;
    char*      value;
  };

  // -- Routines for maintaining and evaluating computed BCs.

  // -- Chicken & egg: build can't be in class constructor 
  //    because we need a Field to do it.  Right?

  void buildComputedBCs (const Field*);

  void maintainFourier  (const int_t, const Field*, const AuxField**,
			 const AuxField**, const bool = true);  
  void maintainPhysical (const Field*, const vector<AuxField*>&, const int_t);
  void evaluateCNBCp    (const int_t, const int_t, const int_t, real_t*);
  void evaluateCEBCp    (const Field*, const int_t, const int_t, 
			 const int_t, real_t*);
  void evaluateCNBCu    (const Field*, const int_t, const int_t, 
			 const int_t, const char, real_t*);
  void accelerate       (const Vector&, const Field*);

private:
  char*              _fields  ; // String containing field names.
  vector<char>       _group   ; // Single-character group tags.
  vector<char*>      _descript; // Group name strings.
  vector<CondRecd*>  _cond    ; // Conditions in storage.
  vector<BCtriple*>  _elmtbc  ; // Group tags for each element-side BC.
  vector<NumberSys*> _numsys  ; // Numbering schemes in storage.
  bool               _axis    ; // Session file declared an axis BC group.
  bool               _outflow ; // Session file declared an outflow BC group.

  void buildnum  (const char*, vector<Element*>&);
  void buildsurf (FEML*, vector<Element*>&);

  // -- Storage of past-time values needed for computed BCs:

  int_t     _nLine;	// Same as for Field storage.
  int_t     _nEdge;     // Number of edges with BCs attached.
  int_t     _nZ;        // Same as for Field storage.
  int_t     _nP;
  int_t     _nTime;

  real_t*** _u;          // (Physical) x velocity component.
  real_t*** _v;          // (Physical) y velocity component.
  real_t*** _w;          // (Physical) z velocity component.

  real_t*** _un;	// (Fourier)  normal velocity u . n
  real_t*** _divu;	// (Fourier)  KINVIS * div(u)
  real_t*** _gradu;	// (Fourier)  KINVIS * normal gradient of velocity 
  real_t*** _hopbc;	// (Fourier)  normal component of [N(u)+f+curlCurl(u)]

  real_t*   _work;      // Computational workspace (scratch).
  real_t*   _fbuf;      // Fourier transform buffer for KE terms.
  real_t*   _ke;	// (Physical) kinetic energy 0.5*(u^2+v^2+w^2)*.
  real_t*   _unp;       // (Physical) u*.n.

  bool      _toggle;    // Toggle switch for Fourier transform of KE.
};

#endif
