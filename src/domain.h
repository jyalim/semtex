#ifndef DOMAIN_H
#define DOMAIN_H

class Domain
// ===========================================================================
// Physical domain storage class.
//
// Since users are expected to program with entities kept at the
// Domain level, all internal storage is exposed to view (like a C
// struct).
//
// Domain fields are written/read untransformed (in physical space).
// ===========================================================================
{
friend istream& operator >> (istream&, Domain&);
friend ostream& operator << (ostream&, Domain&);
public:
  Domain (FEML*, vector<Element*>&, BCmgr*);

  char*                name;	// Session name.
  char*                field; 	// Lower-case single character field names.
  int_t                step;	// Runtime step number.
  real_t               time;	// Simulation time.
  vector<Element*>&    elmt;	// Shared for equal-order interpolations.
  vector<real_t*>      udat;	// Data storage area for solution fields.
  vector<Field*>       u   ;	// Solution fields: velocities, pressure.
  vector<BoundarySys*> b   ;	// Field boundary systems.

  int_t nField     () const { return u.size(); }
  void  report     ();
  void  restart    ();
  void  dump       ();
  void  transform  (const int_t);
};

#endif
