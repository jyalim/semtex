#ifndef NUMBERSYS_H
#define NUMBERSYS_H

class NumberSys
// ===========================================================================
// This class is a holder for Field Element-boundary numbering
// schemes.  Different Fields can potentially have unique
// NumberSyss, or they may share them with other Fields, according
// to distribution of essential BCs.
//
// The numbering system is associated with the Geometry class, but
// augments it by allowing for connectivity and boundary conditions.
//
// Bmask and emask provide a hierarchy of descriptions of element-
// boundary-node essential boundary condition masks:
// 1. bmask is a vector of ints, 4*(np-1)*nel (i.e. nbndry) in length, which
//    describe the CCW traverse of element-boundary nodes.  Values
//    set to 1 indicate that the node has an imposed (essential) BC, values
//    set to 0 indicate that the node either has a natural or no BC.
// 2. emask is the next level of the hierarchy, a vector of ints nel long.
//    Values are set to 1 if the corresponding element has any bmask
//    values set to 1, otherwise set to 0.
//
// Nglobal specifies the number of global nodes, while nsolve ( <=
// nglobal) specifies the number of global nodes which have bmask
// values of zero, i.e. where the nodal values will be solved for
// instead of set explicitly.
//
// Nglobal and nsolve also effectively provide a top level in the
// bmask/emask hierarchy, since if their difference is non-zero then
// at least one bmask value (and emask value) will be 1.
// 
// Btog gives global node numbers to Element-boundary nodes, same
// length as bmask (nbndry).  The optimization level describes the
// scheme used to build btog.
//
// Inv_mass is the inverse of the values of the global mass matrix,
// but on Element boundaries only.  Used for Field smoothing
// operations after derivative operators (if required).  Length =
// nglobal.
//
// A NumberSys is uniquely identified by the entries in btog, or
// equivalently the entries of bmask, together with optimization
// level.
// ===========================================================================
{
friend class BCmgr;
public:
 ~NumberSys () { }; 

  int_t         nGlobal () const { return _nglobal; }
  int_t         nSolve  () const { return _nsolve;  }
  int_t         nBand   () const { return _nbandw;  }

  const char*   fields  () const { return _fields;            }
  const int_t*  bmask   () const { return _bmask;             }
  const int_t*  emask   () const { return _emask;             }
  int_t         fmask   () const { return _nglobal - _nsolve; }
  const int_t*  btog    () const { return _btog;              }
  const real_t* imass   () const { return _imass;             }

private:
  int_t   _optlev ;		// Optimization level used for btog.
  int_t   _nglobal;		// Length of inv_mass.
  int_t   _nsolve ;		// Number of non-masked global nodes.
  int_t   _nbandw ;		// Bandwidth of btog (includes diagonal).

  char*   _fields;		// String with character labels for Fields.
  int_t*  _bmask ;		// 1 for essential-BC nodes, 0 otherwise.
  int_t*  _emask ;		// 1 if associated Element has any esstl set.
  int_t*  _btog  ;		// Gives numbers to all element-boundary knots.
  real_t* _imass ;		// Inverse of global mass matrix;
};

#endif
