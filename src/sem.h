#ifndef SEM_H
#define SEM_H
///////////////////////////////////////////////////////////////////////////////
// Sem.h: main header file for semtex spectral element solvers.
//
// Copyright (c) 1994<-->$Date: 2015/04/20 11:14:18 $, Hugh Blackburn
//
// Conventions: 
// 1. Arrays are 0-offset.
// 2. Internal ident numbers id/ID start at 0.
// 3. Class private variable names start with _.
//
// $Id: sem.h,v 8.1 2015/04/20 11:14:18 hmb Exp $
///////////////////////////////////////////////////////////////////////////////

#include <cstdarg>		/* System C(++) headers.  */
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cctype>
#include <cstring>
#include <climits>
#include <cfloat>
#include <cassert>
#include <ctime>

#include <fenv.h>		/* Std C floating-point error handling. */


#include <iostream>		/* System C++ headers. */
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

#include <vector>
#include <list>
#include <stack>
#include <map>

using namespace std;

#include <cfemdef.h>		/* Semtex library headers.     */
#include <blas.h>
#include <lapack.h>
#include <utility.h>
#include <veclib.h>
#include <femlib.h>

#define ROOTONLY if (Geometry::procID() == 0)
#define VERBOSE  ROOTONLY if (verbose)

#include <feml.h>		/* Semtex src headers. */
#include <geometry.h>
#include <mesh.h>
#include <element.h>

#include <svv.h>

class Boundary;
class BoundarySys;
class AuxField;
class Field;
class Domain;
class BCmgr;
class Statistics;
class HistoryPoint;
class FluidParticle;
class NumberSys;

#include <analysis.h>
#include <auxfield.h>
#include <condition.h>
#include <domain.h>
#include <edge.h>
#include <boundary.h>
#include <bsys.h>
#include <bcmgr.h>
#include <data2df.h>
#include <family.h>
#include <matrix.h>
#include <field.h>
#include <flowrate.h>
#include <history.h>
#include <integration.h>
#include <misc.h>
#include <numbersys.h>
#include <particle.h>
#include <statistics.h>


template<class T> inline void rollv (T* u, const int_t n)
// ===========================================================================
// Stack roll template.  u is an array of type T, with at least n
// elements.  Roll up by one element.
// ===========================================================================
{
  if (n < 2) return;

  T tmp(u[n - 1]);

  for (register int_t q(n - 1); q; q--)
    u[q] = u[q - 1];
  u[0] = tmp;
}


template<class T> inline void rollm (T** u, const int_t m, const int_t n)
// ===========================================================================
// Stack roll template.  u is an matrix of type T, with at least n*m
// elements.  m = number of rows, n = number of columns. Roll up by one row.
// ===========================================================================
{
  if (m < 2) return;
  int_t i, j;
  for (j = 0; j < n; j++) {
    T tmp (u[m-1][j]);
    for (i = m - 1; i; i--)
      u[i][j] = u[i-1][j];
    u[0][j] = tmp;
  }
}

#endif

