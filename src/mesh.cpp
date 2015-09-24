///////////////////////////////////////////////////////////////////////////////
// mesh.C: read information from a FEML stream, provide
// facilities for generation of mesh knots and initial connectivity.
//
// Copyright (c) 1994 <--> $Date: 2015/04/29 02:03:12 $, Hugh Blackburn
//
// Example/required parts of a FEML file:
//
// <NODES NUMBER=9>
// #	tag	x	y	z
// 	1	0.0	0.0	0.0
// 	2	2.0	0.0	0.0
// 	3	4.0	0.0	0.0
// 	4	0.0	0.5	0.0
// 	5	2.0	0.5	0.0
// 	6	4.0	0.5	0.0
// 	7	0.0	1.0	0.0
// 	8	2.0	1.0	0.0
// 	9	4.0	1.0	0.0
// </NODES>
// 
// <ELEMENTS NUMBER=4>
// #	tag	type    nodes
// 	1	<Q>	1 2 5 4	 </Q>
// 	2	<Q>	2 3 6 5  </Q>
// 	3	<Q>	4 5 8 7  </Q>
// 	4	<Q>	5 6 9 8  </Q>
// </ELEMENTS>
// 
// <SURFACES NUMBER=8>
// #	tag	elmt	face	type
//	1	1	1	<B>	w	</B>
//	2	2	1	<B>	w	</B>
//	3	2	2	<B>	o	</B>
//	4	4	2	<B>	o	</B>
//	5	4	3	<B>	w	</B>
//	6	3	3	<B>	w	</B>
//	7	3	4	<B>	v	</B>
//	8	1	4	<B>	v	</B>
// </SURFACES>
//
// Optional:
//
// <CURVES NUMBER=2>
// #    tag     elmt    face    specification
//      1       4       3       <ARC>    -1.0     </ARC>
//      2       2       1       <SPLINE> filename </SPLINE>
// </CURVES>
//
// <GROUPS NUMBER=3>
// #	tag	name	descriptor
// 	1	v	value
// 	2	w	wall
// 	3	o	exit
// </GROUPS>
// 
// <BCS NUMBER=3>
// #	tag	group	number, followed by BCs.
// 	1	v	4
// 			<D>	u = 1.0-4.0*(y-0.5)^2.0	</D>
// 			<D>	v = 0.0			</D>
// 			<D>	w = 0.0			</D>
// 			<H>	p			</H>
// 	2	w	4
// 			<D>	u = 0.0			</D>
// 			<D>	v = 0.0			</D>
// 			<D>	w = 0.0			</D>
// 			<H>	p 			</H>
// 	3	o	4
// 			<N>	u = 0.0			</N>
// 			<N>	v = 0.0			</N>
// 			<N>	w = 0.0			</N>
// 			<D>	p = 0.0			</D>
// </BCS>
// 
// NB: Node, Element and Side IDs are internally held as one less than
// input value, i.e. commence at 0 instead of 1.
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

static char RCS[] = "$Id: mesh.cpp,v 8.2 2015/04/29 02:03:12 hmb Exp $";

#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cctype>
#include <cstring>
#include <climits>
#include <cfloat>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

#include <utility.h>
#include <femlib.h>
#include <veclib.h>
#include <blas.h>
#include <mesh.h>

static inline int_t rma (int_t i, int_t j, int_t n)
// -- Row-major offsetting for 2D arrays with 0-based indexing.
{ return j + i * n; }

#define VERBOSE if (verbose)

Mesh::Mesh (FEML*      f    ,
	    const bool check) :
// ---------------------------------------------------------------------------
// Create a Mesh using information available in feml.
//
// If check is true (default value) then attempt to install all mesh
// information, including surfaces and curved sides.  If it is not,
// then only sufficient information to define the elements is loaded
// (i.e. nodes and element vertices).
// ---------------------------------------------------------------------------
  _feml (*f)
{
  const char  routine[] = "Mesh::Mesh";
  const int_t verbose = Femlib::ivalue ("VERBOSE");
  char        err[StrMax], tag[StrMax];
  int_t       i, j, k, K, Nn;
  Node*       N;
  Elmt*       E;

  _nodeTable .resize (0);
  _elmtTable .resize (0);
  _curveTable.resize (0);

  // -- Input Nodes.

  _nodeTable.resize (Nn = _feml.attribute ("NODES", "NUMBER"));

  if (Nn < 4) {
    sprintf (err, "At least 4 Nodes are needed, found %1d declared", Nn);
    message (routine, err, ERROR);
  }

  VERBOSE cout << routine << ": Reading vertices ... ";

  for (i = 0; i < Nn; i++) {

    while (_feml.stream().peek() == '#') // -- Skip comments.
      _feml.stream().ignore (StrMax, '\n');

    N = new Mesh::Node;
    _feml.stream() >> N -> ID >> N -> loc.x >> N -> loc.y >> N -> loc.z;
    N -> ID--;
    N -> gID      = UNSET;
    N -> periodic = 0;

    if (N -> ID >= Nn) {
      sprintf (err, "Node ID %1d exceeds attribution (%1d)", N -> ID + 1, Nn);
      message (routine, err, ERROR);
    } else 
      _nodeTable[N -> ID] = N;
  }

  VERBOSE cout << "done" << endl;

  // -- Input Elmt corner vertex nodes.
  //    Presently, only quad (<Q>) elements are allowed.
  
  _elmtTable.resize (K = _feml.attribute ("ELEMENTS", "NUMBER"));

  if (K < 1) {
    sprintf (err, "at least 1 element needed, %1d attributed", K);
    message (routine, err, ERROR);
  }

  VERBOSE cout << "  Reading elements ... ";

  for (i = 0; i < K; i++) {

    while (_feml.stream().peek() == '#') // -- Skip comments.
      _feml.stream().ignore (StrMax, '\n');

    E = new Mesh::Elmt;

    _feml.stream() >> E -> ID >> tag;
    E -> ID--;

    if (strcmp (tag, "<Q>") == 0) {
      E -> node.resize (4);
      E -> side.resize (4);

      for (j = 0; j < 4; j++) {
	_feml.stream() >> k;
	k--;
	if (k >= Nn) {
	  sprintf (err, "in element %1d, node tag %1d exceeds maximum (%1d)",
		   E -> ID + 1, k + 1, Nn);
	  message (routine, err, ERROR);
	} else 
	  E -> node[j] = _nodeTable[k];
      }
    } else {
      sprintf (err, "unrecognized element tag: %s", tag);
      message (routine, err, ERROR);
    }

    _feml.stream() >> tag;
    if (strcmp (tag, "</Q>") != 0) {
      sprintf (err, "closing tag </Q> missing for element %1d", E -> ID + 1);
      message (routine, err, ERROR);
    }

    if (E -> ID >= K) {
      sprintf (err, "element ID (%1d) exceeds attribution (%1d)", E -> ID+1,K);
      message (routine, err, ERROR);
    } else
      _elmtTable[E -> ID] = E;
  }

  VERBOSE cout << "done" << endl;

  VERBOSE cout << "  Setting up mesh internal connectivity ... ";

  this -> assemble ();

  VERBOSE cout << "done" << endl;

  if (check) {  
    VERBOSE cout << "  Installing mesh external surface data ... ";
    this -> surfaces ();
    VERBOSE cout << "done" << endl;

    VERBOSE cout << "  Checking mesh connectivity ... ";
    this -> checkAssembly ();
    VERBOSE cout << "done" << endl;

    VERBOSE cout << "  Installing mesh curved sides ... ";
    this -> curves ();
    VERBOSE cout << "done" << endl;
  }
}


void Mesh::assemble (const bool printVacancy)
// ---------------------------------------------------------------------------
// Traverse Elmts and fill in Side-based connectivity information.
//
// On exit, Sides that don't have mateElmt set should correspond to Surfaces.
//
// If printVacancy is true (default: false) then print to cout a
// default set of surfaces with typical wall type BC indication, e.g.
//     idtag elmt side <B> w </B>
// which can be used as the basis of subsequent editing if desired, and/or
// summarizes the elements that need surface information supplied.
// ---------------------------------------------------------------------------
{
  int_t       i, j, r, s, found;
  const int_t Ne = nEl();
  Elmt        *E, *ME;	               // -- M <==> "mate".
  Side        *S, *MS;

  // -- First, build Elmt Sides.

  for (i = 0; i < Ne; i++) {
    E = _elmtTable[i];
    const int_t Nn = E -> nNodes();
    for (j = 0; j < Nn; j++) {
      S = new Side (j, E -> node[j], E -> ccwNode(j));
      S -> thisElmt = E;
      E -> side[j]  = S;
    }
  }

  // -- Now traverse Elmts and build Side--Side connections based on Node
  //    identities.  This can't pick up periodic Nodes, not yet installed.
  //    That happens in surfaces().

  for (i = 0; i < Ne; i++) {
    E = _elmtTable[i];
    const int_t Nn = E -> nNodes();

    for (j = 0; j < Nn; j++) {
      S = E -> side[j];
      found = 0;

      for (r = 0; !found && r < Ne; r++) {
	ME = _elmtTable[r];
	const int_t Nm = ME -> nNodes();
	
	for (s = 0; !found && s < Nm; s++) {
	  MS = ME -> side[s];

	  if ((found = ( S -> startNode == MS -> endNode &&
			 MS -> startNode ==  S -> endNode ))) {
	    S -> mateElmt = ME;
	    S -> mateSide = MS;
	  }
	}
      }
    }
  }

  if (printVacancy) {
 
   // -- Count/print the number of vacancies (number of surfaces needed).
    for (r = 0, i = 0; i < Ne; i++) {
      E = _elmtTable[i];
      const int_t Nn = E -> nNodes();
      for (j = 0; j < Nn; j++) {
	S = E -> side[j];
	if (! S -> mateElmt ) r++;
      }
    }
    cout << "<SURFACES NUMBER=" << r << ">" << endl;

    // -- Now the vacant surface information.
    for (r = 0, i = 0; i < Ne; i++) {
      E = _elmtTable[i];
      const int_t Nn = E -> nNodes();
      for (j = 0; j < Nn; j++) {
	S = E -> side[j];
	if (! S -> mateElmt )
	  cout << ++r
	       << " " << E -> ID+1 << " " << S -> ID+1 
	       << " <B> w </B>" << endl;
      }
    }
    cout << "</SURFACES>" << endl;
  }
}


void Mesh::surfaces ()
// ---------------------------------------------------------------------------
// This section reads FEML "SURFACE" information and uses it to set up
// corresponding element Side storage.
//
// Surface information can either declare
//   a boundary group name <B> group </B>
// or set up
//   a periodic boundary   <P> elmt side </P>.
//
// Periodic boundaries can be set only once (one one side of the matchup).
// ---------------------------------------------------------------------------
{
  const char  routine[] = "Mesh::surfaces";
  const int_t K = _feml.attribute ("SURFACES", "NUMBER");
  char        err[StrMax], tag[StrMax];
  int_t       i, e, s, t;

  for (i = 0; i < K; i++) {

    while (_feml.stream().peek() == '#') // -- Skip comments.
      _feml.stream().ignore (StrMax, '\n');

    // -- Get element and side number information.

    _feml.stream() >> t >> e >> s;
    if (t > K) {
      sprintf (err, "Surface tag no. %1d exceeds attribution (%1d)",
	       t, K);
      message (routine, err, ERROR);
    } else if (e > nEl()) {
      sprintf (err, "Surface %1d element no. %1d too large (%1d)",
	       t, e, nEl());
      message (routine, err, ERROR);
    } else if (s > _elmtTable[e - 1] -> nNodes()) {
      sprintf (err, "Surface %1d elmt %1d side no. %1d too large (%1d)",
	       t, e, s, _elmtTable[e - 1] -> nNodes());
      message (routine, err, ERROR);
    } else if (_elmtTable[e - 1] -> side[s - 1] -> mateElmt) {
      Mesh::Elmt* ME = _elmtTable[e - 1] -> side[s - 1] -> mateElmt;
      Mesh::Side* MS = _elmtTable[e - 1] -> side[s - 1] -> mateSide;
      sprintf (err, "Surface %1d elmt %1d side %1d already set to mate "
	       "elmt %1d side %1d", t, e, s, ME -> ID + 1, MS -> ID + 1);
      message (routine, err, ERROR);
    } 
    
    // -- Set up either a boundary group or a periodic boundary.

    _feml.stream() >> tag;
    e--; s--;

    if (strcmp (tag, "<B>") == 0) {
      
      // -- Boundary group.
      //    Group information for this side should not be set already.

      if (_elmtTable[e] -> side[s] -> group) {
	sprintf (err, "Surface %1d: group already set (%c)",
		 t, _elmtTable[e] -> side[s] -> group);
	message (routine, err, ERROR);
      }
      
      _feml.stream() >> _elmtTable[e] -> side[s] -> group;

      // -- Clean up.

      _feml.stream() >> tag;
      if (strcmp (tag, "</B>") != 0) {
	sprintf (err, "Surface %1d: couldn't close tag <B> with %s", t, tag);
	message (routine, err, ERROR);
      }
      
    } else if (strcmp (tag, "<P>") == 0) {
      
      // -- Periodic.
      //    These are set two at a time (this and indicated mate).
      //    Mate information for both this side and its indicated mate
      //    should not be previously set.

      Side *S, *MS;
      int_t  me,  ms;
      _feml.stream() >> me >> ms;

      if (me < 1 || me > nEl()) {
	sprintf (err, "Surface %1d, mating elmt no. %1d out of range (1--%1d)",
		 t, me, nEl());
	message (routine, err, ERROR);
      } else if (ms < 1 || ms > _elmtTable[e] -> nNodes()) {
	sprintf (err, "Surface %1d, mating side no. %1d out of range (1--%1d)",
		 t, ms, _elmtTable[e] -> nNodes());
	message (routine, err, ERROR);
      } else if (_elmtTable[me - 1] -> side[ms - 1] -> mateElmt ||
		 _elmtTable[me - 1] -> side[ms - 1] -> group    ) {
	sprintf (err, "Surface %1d, mating elmt %1d, side %1d already set",
		 t, me, ms);
	message (routine, err, ERROR);
      }

      me--; ms--;
      
      S  = _elmtTable[e]  -> side[s];
      MS = _elmtTable[me] -> side[ms];

      S  -> mateElmt = _elmtTable[me];
      S  -> mateSide = MS;
      MS -> mateElmt = _elmtTable[e];
      MS -> mateSide = S;

      this -> chooseNode (S -> startNode, MS ->   endNode);
      this -> chooseNode (S ->   endNode, MS -> startNode);

      // -- Clean up.
      
      _feml.stream() >> tag;
      if (strcmp (tag, "</P>") != 0) {
	sprintf (err, "Surface %1d: couldn't close tag <P> with %s", t, tag);
	message (routine, err, ERROR);
      }
      
    } else {
      sprintf (err, "couldn't recognize Surface tag %s", tag);
      message (routine, err, ERROR);
    }
  }
}


void Mesh::chooseNode (Node* N1,
		       Node* N2)
// ---------------------------------------------------------------------------
// It is possible for side-to-side periodicity to miss the fact that a
// node is already periodic with something else, so we test for those
// cases, and choose the lowest Node ID out of available alternatives.
//
// There are 9 possible combinations of N1 -> periodic & N1 -> periodic
// depending on whether they are set or not and whether, if set, they
// are self-referential or not.  In 5 of these cases, a "dangling" pointer
// will be left (one that doesn't point to a self-referential periodicity),
// and these cases have to be fixed by a call to fixPeriodic.
// ---------------------------------------------------------------------------
{
  Node* PN = 0;

  if (!N1->periodic)
    if (!N2->periodic) {
      PN = (N1->ID < N2->ID) ? N1 : N2;
      N1->periodic = N2->periodic = PN;
    } else if (N2->periodic != N2) {
      PN = (N1->ID < N2->periodic->ID) ? N1 : N2->periodic;
      N1->periodic = N2->periodic = N2->periodic->periodic = PN;
    } else {
      PN = (N1->ID < N2->ID) ? N1 : N2;
      N1->periodic = N2->periodic = PN;
      this -> fixPeriodic();
    }
  else if (N1->periodic != N1)
    if (!N2->periodic) {
      PN = (N1->periodic->ID < N2->ID) ? N1->periodic : N2;
      N2->periodic = N1->periodic = N1->periodic->periodic = PN;
    } else if (N2->periodic != N2) {
      PN = (N1->periodic->ID < N2->periodic->ID) ? N1->periodic : N2->periodic;
      N1->periodic = N2->periodic =
	N1->periodic->periodic = N2->periodic->periodic = PN;
    } else {
      PN = (N1->periodic->ID < N2->ID) ? N1->periodic : N2;
      N2->periodic = N1->periodic = N1->periodic->periodic = PN;
      this -> fixPeriodic();
    }
  else {
    if (!N2 -> periodic) {
      PN = (N1->ID < N2->ID) ? N1 : N2;
      N1->periodic = N2->periodic = PN;
      fixPeriodic();
    } else if (N2 -> periodic != N2) {
      PN = (N1->ID < N2->periodic->ID) ? N1 : N2;
      N2->periodic = N1->periodic = N1->periodic->periodic = PN;
      fixPeriodic();
    } else {
      PN = (N1->ID < N2->ID) ? N1 : N2;
      N1->periodic = N2->periodic = PN;
      this -> fixPeriodic();
    }
  }
  
  if (!PN) message ("Mesh::chooseNode", "NEVER HAPPEN", ERROR);
}


void Mesh::fixPeriodic()
// ---------------------------------------------------------------------------
// Traverse all Nodes and fix up any Nodes whose periodic value doesn't
// terminate self-referentially (not recursively as this is not needed).
// ---------------------------------------------------------------------------
{
  const int_t    N = _nodeTable.size();
  register int_t i;
  Node           *np, *npp;

  for (i = 0; i < N; i++) {
    np  = _nodeTable[i];
    npp = np -> periodic;
    if (npp && npp -> periodic != npp) {
      npp = npp -> periodic;
      if (!npp) message ("Mesh::fixPeriodic", "NEVER HAPPEN", ERROR);
    }
  }
}


void Mesh::showAssembly (Mesh& m)
// ---------------------------------------------------------------------------
// Static debugging function: Print edge--edge connectivity information.
// ---------------------------------------------------------------------------
{
  int_t i, j, Ne, Nn;
  Elmt* E;
  Side* S;

  Ne = m.nEl();
  cout << "# " << Ne << " Elmts" << endl;

  for (i = 0; i < Ne; i++) {
    E  = m._elmtTable[i];
    Nn = E -> nNodes ();

    cout << "# Elmt: " << E -> ID + 1 << ", Vertices:";
    for (j = 0; j < Nn; j++)
      cout << " " << E -> node[j] -> ID + 1;

    cout << ", Mating:";
    for (j = 0; j < Nn; j++) {
      S = E -> side[j];
      cout << '\t' << S -> ID + 1;
      cout << "->";
      if   (!S -> mateElmt) cout << S -> group;
      else                  cout << S -> mateElmt -> ID + 1 << "."
				 << S -> mateSide -> ID + 1;
    }
    cout << endl;
  }
}


void Mesh::checkAssembly()
// ---------------------------------------------------------------------------
// All element sides have to either mate an adjoining element or fall on
// a boundary.  Check it out.  But surfaces() must have been called first.
// ---------------------------------------------------------------------------
{
  char           routine[] = "Mesh::checkAssembly", err[StrMax];
  const int_t    Ne = nEl();
  Elmt*          E;
  Side*          S;
  register int_t i, j;
  bool           OK = true;

  for (i = 0; i < Ne; i++) {
    E = _elmtTable[i];
    const int_t Ns = E -> nNodes();
    for (j = 0; j < Ns; j++) {
      S = E -> side[j];
      if (S -> mateSide == 0) {
	sprintf (err, "Elmt %1d Side %1d not set",
		 S -> thisElmt -> ID + 1, S -> ID + 1);
	message (routine, err, WARNING);
	OK = false;
      }
    }
  }
  
  if (!OK) message (routine, "some element edges not accounted for", ERROR);

  if (Femlib::ivalue ("VERBOSE") > 1) {
    cout << endl << "# Summary:" << endl;
    showAssembly (*this);
  }
}


void Mesh::curves ()
// ---------------------------------------------------------------------------
// Read in curved edge information and store in Mesh _curveTable.
//
// Curved edges are specified by lines like:
//   curveID  elementID  sideID <C> ... </C>
//
// Tags <C> and </C> delimit the kind of curve to be defined; presently
// the only defined kind is <ARC>.  Between the delimiters the number
// of parameters is user-defined.
// ---------------------------------------------------------------------------
{
  if (!_feml.seek ("CURVES")) return;
  
  const char routine[] = "Mesh::curves";
  char       err[StrMax], buf[StrMax];
  int_t      i, K, id, elmt, side, ns;
  Curve*     C;
  Side*      S;

  _curveTable.resize (K = _feml.attribute ("CURVES", "NUMBER"));

  for (i = 0; i < K; i++) {
    _feml.stream() >> id >> elmt >> side;

    if (id > K) {
      sprintf (err, "Curve ID %1d exceeds attribution (%1d)", id, K);
      message (routine, err, ERROR);
    } else if (elmt > nEl()) {
      sprintf (err, "Curve ID %1d, Elmt no. %1d too large (%1d)", id, elmt, K);
      message (routine, err, ERROR);
    } else if (side > (ns = _elmtTable[elmt - 1] -> nNodes())) {
      sprintf (err, "Curve ID %1d, Side no. %1d too large (%1d)", id,side, ns);
      message (routine, err, ERROR);
    }
    
    S = _elmtTable[elmt - 1] -> side[side - 1];

    _feml.stream() >> buf;

    if (strcmp (buf, "<ARC>") == 0) {
      real_t radius;
      _feml.stream() >> radius;

      C = new CircularArc (id, S, radius);

      _feml.stream() >> buf;
      if (strcmp (buf, "</ARC>") != 0) {
	sprintf (err, "Curve ID %1d, can't close <ARC> with </ARC>", id);
	message (routine, err, ERROR);
      }
    } else if (strcmp (buf, "<SPLINE>") == 0) {
      char filename[StrMax];
      _feml.stream() >> filename;

      C = new Spline (id, S, filename);

      _feml.stream() >> buf;
      if (strcmp (buf, "</SPLINE>") != 0) {
	sprintf (err, "Curve ID %1d, can't close <SPLINE> with </SPLINE>", id);
	message (routine, err, ERROR);
      }
    } else {
      sprintf (err, "Curve %1d, unknown curve kind %s", id, buf);
      message (routine, err, ERROR);
    }

    _curveTable[i] = C;
  }
}


Point Mesh::Elmt::centroid () const
// ---------------------------------------------------------------------------
// Return point that is centroid of element Node points.
// ---------------------------------------------------------------------------
{
  register int_t i;
  const    int_t K = nNodes();
  Point    C = {0.0, 0.0, 0.0};

  for (i = 0; i < K; i++) {
    Point P = node[i] -> loc;
    C.x += P.x;
    C.y += P.y;
  }

  C.x /= K;
  C.y /= K;

  return C;
}


void Mesh::meshSide (const int_t   np     ,
		     const int_t   elmt   ,
		     const int_t   side   ,
		     const real_t* spacing,
		     Point*        knot   ) const
// ---------------------------------------------------------------------------
// If a curved side can be identified for the nominated element and side,
// compute the points using appropriate routine.  Otherwise compute points
// along a straight side.
//
// Spacing gives location of knots in master coordinates [-1, 1].
// ---------------------------------------------------------------------------
{
  const char     routine[] = "Mesh::meshSide";
  const int_t    Nc = _curveTable.size();
  const int_t    Ne = _elmtTable .size();
  register int_t i;

  if (np < 2) message (routine, "must have at least two points", ERROR);

  for (i = 0; i < Nc; i++) {
    Curve* C = _curveTable[i];
    if (C -> ismatch (elmt, side)) {
      C -> compute (np, spacing, knot);
      return;
    }
  }

  // -- Fall though default: straight line.

  const Side* S  = _elmtTable[elmt] -> side[side];
  const Point P1 = S -> startNode -> loc;
  const Point P2 = S -> endNode   -> loc;
  const real_t  dx = P2.x - P1.x;
  const real_t  dy = P2.y - P1.y;

  for (i = 0; i < np; i++) {
    knot[i].x = P1.x + dx * 0.5 * (spacing[i] + 1.0);
    knot[i].y = P1.y + dy * 0.5 * (spacing[i] + 1.0);
  }
  return;
}


void Mesh::meshElmt (const int_t   ID,
		     const int_t   np,
		     const real_t* zr,
		     const real_t* zs,
		     real_t*       x ,
		     real_t*       y ) const
// ---------------------------------------------------------------------------
// This is a routine for use of Mesh by other classes (like Element).
//
// Generate mesh points for Elmt No ID (IDs begin at 0).  Generate
// element-edge points, then internal points using a Coons patch.
// Inputs zr, zs contain the spacing of edge knot points along
// interval [-1, 1], in the "r" and "s" directions.
//
// For a quad mesh, equal-order on each side, x & y have row-major ordering.
//
// Because the element edges will be mapped in CCW order, we have to
// reverse the direction of the input spacings zr and zs on sides 2
// and 3.
//
// Offset in x and y directions by X_SHIFT and Y_SHIFT TOKEN values.
// Scale  in x and y directions by X_SCALE and Y_SCALE TOKEN values.
// ---------------------------------------------------------------------------
{
  const char     routine[] = "Mesh::meshElmt";
  const int_t    nm     = np - 1;
  const int_t    ns     = _elmtTable[ID] -> nNodes();
  const real_t   x_shft = Femlib::value ("X_SHIFT");
  const real_t   y_shft = Femlib::value ("Y_SHIFT");
  const real_t   x_scal = Femlib::value ("X_SCALE");
  const real_t   y_scal = Femlib::value ("Y_SCALE");
  register int_t i, j;
  vector<Point>  P (np);
  vector<real_t> work (np);

  // -- Compute and load peripheral points.

  for (j = 0; j < ns; j++) {
    switch (j) {
    case 0:
      this -> meshSide (np, ID, j, zr, &P[0]);
      for (i = 0; i < nm; i++) {
	x[rma (     0,      i, np)] = P[i].x;
	y[rma (     0,      i, np)] = P[i].y;
      }
      break;
    case 1:
      this -> meshSide (np, ID, j, zs, &P[0]);
      for (i = 0; i < nm; i++) {
	x[rma (     i,     nm, np)] = P[i].x;
	y[rma (     i,     nm, np)] = P[i].y;
      }
      break;
    case 2:
      copy (zr, zr+np, work.begin());
      reverse (work.begin(), work.end());
      Veclib::neg (np, &work[0], 1);
      this -> meshSide (np, ID, j, &work[0], &P[0]);
      for (i = 0; i < nm; i++) {
	x[rma (    nm, nm - i, np)] = P[i].x;
	y[rma (    nm, nm - i, np)] = P[i].y;
      }
      break;
    case 3:
      copy (zs, zs+np, work.begin());
      reverse (work.begin(), work.end());
      Veclib::neg (np, &work[0], 1);
      this -> meshSide (np, ID, j, &work[0], &P[0]);
      for (i = 0; i < nm; i++) {
	x[rma (nm - i,      0, np)] = P[i].x;
	y[rma (nm - i,      0, np)] = P[i].y;
      }
      break;
    default:
      message (routine, "never happen", ERROR);
      break;
    }
  }

  // -- Coons patch on (-1, 1) X (-1, 1) to make interior points.

  for (i = 1; i < nm; i++)
    for (j = 1; j < nm; j++) {

      x[rma (i, j, np)] = 0.50 * ( (1.0 - zr[j]) * x[rma ( i,  0, np)] +
				   (1.0 + zr[j]) * x[rma ( i, nm, np)] +
				   (1.0 - zs[i]) * x[rma ( 0,  j, np)] +
				   (1.0 + zs[i]) * x[rma (nm,  j, np)] )
	     
	 - 0.25 * ( (1.0 - zr[j]) * (1.0 - zs[i]) * x[rma ( 0,  0, np)] +
	            (1.0 - zr[j]) * (1.0 + zs[i]) * x[rma (nm,  0, np)] +
		    (1.0 - zs[i]) * (1.0 + zr[j]) * x[rma ( 0, nm, np)] +
		    (1.0 + zs[i]) * (1.0 + zr[j]) * x[rma (nm, nm, np)] );

      y[rma (i, j, np)] = 0.50 * ( (1.0 - zr[j]) * y[rma ( i,  0, np)] +
				   (1.0 + zr[j]) * y[rma ( i, nm, np)] +
				   (1.0 - zs[i]) * y[rma ( 0,  j, np)] +
				   (1.0 + zs[i]) * y[rma (nm,  j, np)] )
	     
	 - 0.25 * ( (1.0 - zr[j]) * (1.0 - zs[i]) * y[rma ( 0,  0, np)] +
                    (1.0 - zr[j]) * (1.0 + zs[i]) * y[rma (nm,  0, np)] +
		    (1.0 - zs[i]) * (1.0 + zr[j]) * y[rma ( 0, nm, np)] +
		    (1.0 + zs[i]) * (1.0 + zr[j]) * y[rma (nm, nm, np)] );
    }

  // -- Carry out x--y shifting if required.

  if (fabs (x_shft) > EPSDP) Veclib::sadd (np*np, x_shft, x, 1, x, 1);
  if (fabs (y_shft) > EPSDP) Veclib::sadd (np*np, y_shft, y, 1, y, 1);

  // -- Carry out x--y scaling if required.

  if (fabs (x_scal - 1.0) > EPSDP) Blas::scal (np*np, x_scal, x, 1);
  if (fabs (y_scal - 1.0) > EPSDP) Blas::scal (np*np, y_scal, y, 1);
}


int_t Mesh::buildMap (const int_t np ,
		      int_t*      map)
// ---------------------------------------------------------------------------
// Generate connectivity (i.e. global knot numbers) for a mesh with np
// knot points (i.e. Lagrange knots) along each element side, ignoring
// internal points (i.e. generate connectivity for static-condensation form).
//
// Fill map (element-by-element storage of these global numbers) for whole
// mesh: for a mesh of quad elements, map must hold 4*(np-1)*nEl int_ts.
// Return the number of global knots (maximum global knot number + 1). 
//
// NB: np >= 2, also global numbers generated here start at 0.
// NB: this connectivity information is generated without reference to BCs.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Mesh::buildMap";

  if (np < 2) message (routine, "need at least 2 knots", ERROR);

  // -- Create element-side based gID storage, if required, & initialize gIDs.
  
  const int_t    nel = nEl(), ni = np - 2;
  register int_t i, j, k, ns;
  int_t          nGid = 0, nb = 0;
  Elmt*          E;
  Side*          S;

  // -- Allocate space, unset all knot numbers.

  for (i = 0; i < nel; i++) {
    E  = _elmtTable[i];
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S = E -> side[j];
      S -> gID.resize (ni);      
      S -> startNode -> gID = UNSET;
      S -> endNode   -> gID = UNSET;
      if (ni) Veclib::fill (ni, UNSET, &S -> gID[0], 1);
    }
  }

  // -- Generate connectivity information.

  for (i = 0; i < nel; i++) {
    E  = _elmtTable[i];
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S = E -> side[j];
      S -> connect (ni, nGid);
    }
  }

  // -- Fill map.

  for (i = 0; i < nel; i++) {
    E  = _elmtTable[i];
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S = E -> side[j];
      map[nb++] = S -> startNode -> gID;
      for (k = 0; k < ni; k++)
	map[nb++] = S -> gID[k];
    }
  }

  // -- Deallocate internal knot number storage.

  for (i = 0; i < nel; i++) {
    E  = _elmtTable[i];
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S = E -> side[j];
      S -> gID.resize (0);      
      S -> startNode -> gID = UNSET;
      S -> endNode   -> gID = UNSET;
    }
  }

  return nGid;
}


void Mesh::Side::connect (const int_t ni ,
			  int_t&      gid)
// ---------------------------------------------------------------------------
// Fill in connectivity for this element side, updating global number gid.
// ---------------------------------------------------------------------------
{
  register int_t i, k;
  register Side* otherSide;

  if (startNode -> periodic) {
    if (startNode -> periodic -> gID == UNSET)
      startNode -> periodic -> gID = gid++;
    if (startNode -> gID == UNSET)
      startNode -> gID = startNode -> periodic -> gID;
  }

  if (startNode -> gID == UNSET)
    startNode -> gID = gid++;

  if (ni) {			// -- Do side-internal gids.
    if (mateElmt) {
      otherSide = mateSide;
      if (otherSide -> gID[0] == UNSET)
	for (i = 0; i < ni; i++)
	  gID[i] = gid++;
      else
	for (i = 0, k = ni - 1; i < ni; i++, k--)
	  gID[i] = otherSide -> gID[k];
    } else
      for (i = 0; i < ni; i++)
	gID[i] = gid++;
  }

  if (endNode -> periodic) {
    if (endNode -> periodic -> gID == UNSET)
      endNode -> periodic -> gID = gid++;
    if (endNode -> gID == UNSET)
      endNode -> gID = endNode -> periodic -> gID;
  }

  if (endNode -> gID == UNSET)
    endNode -> gID = gid++;
}


void Mesh::printNek () const
// ---------------------------------------------------------------------------
// Print out mesh information in NEKTON format.
// ---------------------------------------------------------------------------
{
  const char    routine[] = "Mesh::printNek";
  char          buf[StrMax];
  string        err;
  ostringstream os(err);

  int_t      i, j, ns, nel = nEl();
  float      vbc;
  Elmt       *E, *ME;
  Side       *S;

  cout.precision (5);
  cout.setf      (ios::scientific,ios::floatfield);

  // -- Elements.

  cout << setw(14) << 10.0
       << setw(14) << 10.0
       << setw(14) << 0.0 
       << setw(14) << 0.0
       << " XFAC,YFAC,XZERO,YZERO"
       << endl;

  cout << "**MESH DATA** 1st line is X of corner 1,2,3,4. 2nd line is Y."
       << endl;

  cout << setw(10) << nel
       << setw(10) << 2 
       << setw(10) << nel
       << " NEL,NDIM,NELV"
       << endl;

  for (i = 0; i < nel; i++) {
    E = _elmtTable[i];

    cout << "ELEMENT   "
         << setw(10) << E -> ID + 1
         << " [  1A]  GROUP 0"
         << endl;

    ns = E -> nNodes();

    for (j = 0; j < ns; j++)
      cout << setw(14) << E -> node[j] -> loc.x;
    cout << endl;

    for (j = 0; j < ns; j++)
      cout << setw(14) << E -> node[j] -> loc.y;
    cout << endl;
  }

  // -- Curved sides.

  ns = _curveTable.size();
  cout << "***** CURVED SIDE DATA *****" << endl;
  cout << setw(5) << ns
       << " Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE"
       << endl;
  for (i = 0; i < ns; i++)
    _curveTable[i] -> printNek ();

  // -- Boundary conditions.

  cout << "***** BOUNDARY CONDITIONS *****" << endl;
  cout << "***** FLUID BOUNDARY CONDITIONS *****" << endl;

  for (i = 0; i < nel; i++) {
    E  = _elmtTable[i];
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S  = E -> side[j];
      ME = S -> mateElmt;
      if (ME) {
	cout << "E  "
	     << setw (5)  << E -> ID + 1
	     << setw (3)  << S -> ID + 1
	     << setw (14) << 1.0*ME -> ID + 1
	     << setw (14) << 1.0*S  -> mateSide -> ID + 1
	     << setw (14) << 1.0
	     << endl;
      } else {
	describeGrp (S -> group, buf);
	if (strstr (buf, "value")) {
	  cout << "V  "
	       << setw (5) << E -> ID + 1
	       << setw (3) << S -> ID + 1;
	  describeBC (S -> group, 'u', buf);
	  sscanf     (buf, "%*s %*s %f", &vbc);
	  cout << setw (14) << vbc;
	  describeBC (S -> group, 'v', buf);
	  sscanf     (buf, "%*s %*s %f", &vbc);
	  cout << setw (14) << vbc;
	  if (Femlib::ivalue ("N_Z") > 1) {
	    describeBC (S -> group, 'w', buf);
	    sscanf     (buf, "%*s %*s %f", &vbc);
	  } else
	    vbc = 0.0;
	  cout << setw (14) << vbc;
	  cout << endl;
	} else if  (strstr (buf, "wall")) {
	  cout << "W  "
	       << setw (5)  << E -> ID + 1
	       << setw (3)  << S -> ID + 1
	       << setw (14) << 0.0
	       << setw (14) << 0.0
	       << setw (14) << 0.0
	       << endl;
	} else if (strstr (buf, "outflow")) {
	  cout << "O  "
	       << setw (5)  << E -> ID + 1
	       << setw (3)  << S -> ID + 1
	       << setw (14) << 0.0
	       << setw (14) << 0.0
	       << setw (14) << 0.0
	       << endl;
	} else {
	  os << "Elmt " << E -> ID + 1 << " side " << S -> ID + 1
	     << " --- B.C. type "  << buf
	     << " not implemented" << ends;
	  message (routine, os.str().c_str(), WARNING);
	}
      }
    }
  }
}


void Mesh::describeGrp (char  G,
			char* S) const
// ---------------------------------------------------------------------------
// Search feml file info for string descriptor matching G, load into S.
// ---------------------------------------------------------------------------
{
  const char   routine[] = "Mesh::describeGrp";
  char         groupc, err[StrMax], buf[StrMax];
  int_t        i, id;
  bool         found = false;
  const int_t  N = _feml.attribute ("GROUPS", "NUMBER");
  
  for (i = 0; !found && i < N; i++) {
    while (_feml.stream().peek() == '#') // -- Skip comments.
      _feml.stream().ignore (StrMax, '\n');
    _feml.stream() >> id >> groupc >> buf;
    if ((found = (groupc == G))) strcpy (S, buf);
  }

  if (!found) {
    sprintf (err, "no group found to match '%c'", G);
    message (routine, err, ERROR);
  }
}


void Mesh::describeBC (char  grp,
		       char  fld, 
		       char* tgt) const
// ---------------------------------------------------------------------------
// Find BC description string matching group 'grp' and Field 'fld',
// load into tgt.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "Mesh::describeBC";
  const int_t N = _feml.attribute ("BCS", "NUMBER");
  char        eql, groupc, fieldc, err[StrMax], buf[StrMax];
  int_t       i, j, id, nbcs;
  bool        found = false;

  for (i = 0; !found && i < N; i++) {

    while (_feml.stream().peek() == '#') // -- Skip comments.
      _feml.stream().ignore (StrMax, '\n');

    _feml.stream() >> id >> groupc >> nbcs;

    for (j = 0; !found && j < nbcs; j++) {

      // -- Trash tag. (i.e. ignore type of BC. This is an error!).

      _feml.stream() >> buf;

      // -- Check for match and take appropriate action if true.

      _feml.stream() >> fieldc;

      if ((found = (groupc == grp) && (fieldc == fld))) {
	_feml.stream() >> eql;
	if (eql == '=') {
	  tgt[0] = fld;
	  tgt[1] = '\0';
	  strcat (tgt, " = ");
	  _feml.stream() >> buf;
	  strcat (tgt, buf);
	} else {
	  sprintf (err, "Group '%c', Field '%c', expected '=', got '%c",
		   grp, fld, eql);
	  message (routine, err, ERROR);
	}
      } else {
	_feml.stream().ignore (StrMax, '\n');
      }
    }
  }
      
  if (!found) {
    sprintf (err, "couldn't find BC to match Group '%c', Field '%c'",
	     grp, fld);
    message (routine, err, ERROR);
  }
}


void Mesh::buildMask (const int_t np  ,
		      const char  fld ,
		      int_t*      mask)
// ---------------------------------------------------------------------------
// This routine generates an int_t mask (0/1) vector for
// element-boundary nodes.  For any location that corresponds to a
// domain boundary with an essential boundary condition and for field
// name "fld", the corresponding mask value will be 1 -- this tags the
// corresponding node for lifting out of the field solution, since the
// field value will be set, rather than solved as part of the system
// of equations (i.e. it forms part of the RHS of a matrix system
// equation, rather than the LHS).  All other locations will be 0
// (i.e., unmasked).
//
// For quads, mask is 4 * nel * (np - 1) long, same as input for buildMap.
// Use is made of the fact that on BCs, there are no mating sides, hence
// no need to set mask on mating sides.
//
// NOTE that the default behaviour for any type of BC that is not <D>
// or <A> is for the mask to be 0.
//
// If fld is 'U', 'v, 'w', 'P' or 'C', set mask for essential BCs on
// symmetry axis.
//
// If fld is 'p' or 'P', set mask for BC of type <O> (outflow).
// ---------------------------------------------------------------------------
{
  const char routine[] = "Mesh::buildMask";

  if (np < 2) message (routine, "need at least 2 knots", ERROR);

  register int_t i, j, k, ns, nb = 0;
  const int_t    nel   = nEl(), ni = np - 2;
  const int_t    axisE = strchr ("UvwPC", fld) != 0;
  Elmt*          E;
  Side*          S;

  // -- Allocate space, unmask all gIDs.

  for (i = 0; i < nel; i++) {
    E  = _elmtTable[i];
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S = E -> side[j];
      S -> gID.resize (ni);      
      S -> startNode -> gID = 0;
      S -> endNode   -> gID = 0;
      if (ni) Veclib::fill (ni, 0, &S -> gID[0], 1);
    }
  }

  // -- Switch on gID in appropriate locations, for D, A <==> Dirichlet BCs.

  for (i = 0; i < nel; i++) {
    E  = _elmtTable[i];
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S = E -> side[j];
      if (!(S -> mateElmt)) {
	if (
	                   matchBC (S -> group, tolower (fld), 'D')  ||
	    (axisE      && matchBC (S -> group, tolower (fld), 'A')) ||
 	    ((fld == 'p' || fld == 'P') 
                        && matchBC (S -> group, tolower (fld), 'O'))
	    ) {
	  S -> startNode -> gID = 1;
	  S -> endNode   -> gID = 1;
	  if (ni) Veclib::fill (ni, 1, &S -> gID[0], 1);
	  if (S -> startNode -> periodic) S -> startNode -> periodic -> gID =1;
	  if (S -> endNode   -> periodic) S -> endNode   -> periodic -> gID =1;
	}
      }
    }
  }

  // -- Traverse mesh and load mask values.

  for (i = 0; i < nel; i++) {
    E  = _elmtTable[i];
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S = E -> side[j];
      mask[nb++] = S -> startNode -> gID;
      for (k = 0; k < ni; k++)
	mask[nb++] = S -> gID[k];
    }
  }

  // -- Deallocate internal knot number storage.

  for (i = 0; i < nel; i++) {
    E  = _elmtTable[i];
    ns = E -> nNodes();
    for (j = 0; j < ns; j++) {
      S = E -> side[j];
      S -> gID.resize (0);      
      S -> startNode -> gID = UNSET;
      S -> endNode   -> gID = UNSET;
    }
  }
}


bool Mesh::matchBC (const char grp,
		    const char fld,
		    const char bcd)
// ---------------------------------------------------------------------------
// From FEML BC information, return true if the boundary condition
// kind shown for group 'grp' and field 'fld' is of type 'bcd'.
//
// Input parameter grp is the group character set in the SURFACES
// section for a mesh side.  Input parameter 'fld' is the name of the
// field for which we are currently setting up a mask vector.
// 
// Presently allowed types for bcd (these correspond to FEML BC tag names):
//   D <==> Dirichlet/Essential.
//   N <==> Neumann/Natural.
//   M <==> Mixed/Robin BC.
//   H <==> "High-order" (computed, natural) pressure BC.        See KIO91.
//   A <==> "Axis" (selected, natural/essential) BC.             See BS04.
//   O <==> "Outflow": computed essential BC for pressure only.  See DKC14.
// ---------------------------------------------------------------------------
{
  const int_t N = _feml.attribute ("BCS", "NUMBER");
  char        groupc, fieldc, buf[StrMax];
  int_t       i, j, id, nbcs;

  for (i = 0; i < N; i++) {		 // -- Read through <BC> section.
    while (_feml.stream().peek() == '#') // -- Skip comments.
      _feml.stream().ignore (StrMax, '\n');

    _feml.stream() >> id >> groupc >> nbcs;

    if (groupc != grp) {
      _feml.stream().ignore (StrMax, '\n');
      for (j = 0; j < nbcs; j++)
	_feml.stream().getline (buf, StrMax);
    } else {			// -- Deal only with set of BCs for "grp".
      for (j = 0; j < nbcs; j++) {
	_feml.stream() >> buf >> fieldc;
	if (buf[1] == bcd && fieldc == fld)
	  return true;
	else
	  _feml.stream().ignore (StrMax, '\n');
      }
    }
  }

  return false;
}


void Mesh::extent (Point& lo,
		   Point& hi) const
// ---------------------------------------------------------------------------
// Traverse element vertices and extract maximum and minumum x, y locations.
// Could do this with the vertex list but some vertices might be unused.
// ---------------------------------------------------------------------------
{
  int_t  i, j;
  Elmt*  E;
  real_t x, xmax = -FLT_MAX, xmin = FLT_MAX;
  real_t y, ymax = -FLT_MAX, ymin = FLT_MAX;

  const int_t Ne = nEl();  
  for (i = 0; i < Ne; i++) {
    E = _elmtTable[i];
    const int_t Nn = E -> nNodes();
    for (j = 0; j < Nn; j++) {
      x = E -> node[j] -> loc.x;
      y = E -> node[j] -> loc.y;
      xmin = min (xmin, x);
      xmax = max (xmax, x);
      ymin = min (ymin, y);
      ymax = max (ymax, y);
    }
  }
  
  lo.x = xmin;
  lo.y = ymin;
  lo.z = 0.0;

  hi.x = xmax;
  hi.y = ymax;
  hi.z = 0.0;
}


//////////////////////////////////////////////////////////////////////////////
// Routines to deal with circular arc mesh sides.
//////////////////////////////////////////////////////////////////////////////


CircularArc::CircularArc (const int_t  id,
			  Mesh::Side*  S ,
			  const real_t R )
// ---------------------------------------------------------------------------
// Constructor for CircularArc.  R is the radius of arc, and its sign
// specifies the convexity of the element edge.
//
// R +ve ==> arc increases area enclosed by element (cf straight line),
// R -ve ==> arc decreases area enclosed by element.
// ---------------------------------------------------------------------------
{
  const char routine[] = "CircularArc::CircularArc";
  char       err[StrMax];

  curveSide = S;
  convexity = (R < 0.0) ? -1 : 1;
  radius    = fabs (R);

  Point P1  = curveSide -> startNode -> loc;
  Point P2  = curveSide -> endNode   -> loc;
  Point unitNormal, link, midpoint, centroid = {0.0, 0.0, 0.0};
  real_t  dx, dy, l, sign = 0.0;

  midpoint.x   = 0.5 * (P2.x + P1.x);
  midpoint.y   = 0.5 * (P2.y + P1.y);
  dx           =        P2.x - P1.x;
  dy           =        P2.y - P1.y;
  l            = hypot (dx, dy);
  unitNormal.x = -dy / l;
  unitNormal.y =  dx / l;

  if (2.0 * radius < l) {
    sprintf (err, "curve %1d:\narc, radius %f, can't span nodes %1d & %1d",
	     id, radius, 
	     curveSide -> startNode -> ID + 1, curveSide -> endNode -> ID + 1);
    message (routine, err, ERROR);
  } else
    semiangle = asin (0.5*l / radius);

  centroid = curveSide -> thisElmt -> centroid ();
  
  link.x = centroid.x - midpoint.x;
  link.y = centroid.y - midpoint.y;

  // -- Sign +1 if centre lies in direction of centroid from side midpoint.

  sign = link.x * unitNormal.x + link.y * unitNormal.y;
  sign = convexity * sign / fabs (sign);

  centre.x = midpoint.x + sign * cos (semiangle) * radius * unitNormal.x;
  centre.y = midpoint.y + sign * cos (semiangle) * radius * unitNormal.y;
}


void CircularArc::compute (const int_t   np     ,
			   const real_t* spacing,
			   Point*        knot   ) const
// ---------------------------------------------------------------------------
// Distribute np knots along arc according to spacing on -1, 1.
// ---------------------------------------------------------------------------
{
  const int_t nm = np - 1;
  Point       P1 = curveSide -> startNode -> loc;
  Point       P2 = curveSide -> endNode   -> loc;
  real_t      theta1, theta2, dtheta, phi;
  int_t       i;

  theta1 = atan2 (P1.y - centre.y, P1.x - centre.x);
  theta2 = atan2 (P2.y - centre.y, P2.x - centre.x);
  dtheta = theta2 - theta1;

  if (fabs (dtheta) > 2.0*semiangle + EPSSP)
    dtheta += (dtheta < 0.0) ? TWOPI : -TWOPI;

  knot[ 0].x = P1.x;  knot[ 0].y  = P1.y;
  knot[nm].x = P2.x;  knot[nm].y = P2.y;

  for (i = 1; i < nm; i++) {
    phi = theta1 + dtheta * 0.5 * (spacing[i] + 1.0);
    knot[i].x = centre.x + radius * cos (phi);
    knot[i].y = centre.y + radius * sin (phi);
  }
}


void CircularArc::printNek () const
// ---------------------------------------------------------------------------
// Print out information in NEKTON format.
// ---------------------------------------------------------------------------
{
  cout << setw (2)  << curveSide -> ID + 1
       << setw (5)  << curveSide -> thisElmt -> ID + 1
       << setw (14) << 1.0*convexity*radius
       << setw (14) << 0.0
       << setw (14) << 0.0
       << setw (14) << 0.0
       << setw (14) << 0.0
       << " C"
       << endl; 
}


//////////////////////////////////////////////////////////////////////////////
// Routines to deal with spline curve mesh sides.
//////////////////////////////////////////////////////////////////////////////

static Point             ga, gp;
static spline2D*         gs;
static vector<spline2D*> gcurve;


Spline::Spline (const int_t id      ,
		Mesh::Side* S       ,
		const char* filename)
// ---------------------------------------------------------------------------
// Constructor, calls getGeom, computes starting and ending arc parameters.
// ---------------------------------------------------------------------------
{
  strcpy ((_name = new char [strlen(filename) + 1]), filename);
  curveSide = S;
  _geom     = Spline::getGeom (filename);

  const char  routine[] = "Spline::Spline";
  const int_t verbose = Femlib::ivalue ("VERBOSE");
  Point       p1 = S -> startNode -> loc;
  Point       p2 = S -> endNode   -> loc;

  VERBOSE cerr << routine << ": --" << endl;

  // -- Set global variables used in locating arc parameters.

  gs   = _geom;

  gp   = p1;
  ga.x = p1.x - (p2.y - p1.y);
  ga.y = p1.y + (p2.x - p1.x);
  
  _startarc = this -> arcCoord ();

  VERBOSE cerr << "    location of p1: " << _startarc << endl;

  gp   = p2;
  ga.x = p2.x - (p2.y - p1.y);
  ga.y = p2.y + (p2.x - p1.x);
  
  _endarc = this -> arcCoord ();

  VERBOSE cerr << "    location of p2: " << _endarc << endl;
}


void Spline::compute (const int_t   np     ,
		      const real_t* spacing,
		      Point*        knot   ) const
// ---------------------------------------------------------------------------
// Compute knot points for spline curve.
// ---------------------------------------------------------------------------
{
  const int_t  nm = np - 1;
  const real_t darc = _endarc - _startarc;
  Point        P1 = curveSide -> startNode -> loc;
  Point        P2 = curveSide -> endNode   -> loc;
  real_t       s;
  int_t        i;

  knot[ 0].x = P1.x;  knot[ 0].y = P1.y;
  knot[nm].x = P2.x;  knot[nm].y = P2.y;

  for (i = 1; i < nm; i++) {
    s = _startarc + darc * 0.5 * (spacing[i] + 1.0);
    knot[i].x = Veclib::splint (_geom->x.size(), s, &_geom->arclen[0],
				&_geom->x[0], &_geom->sx[0]);
    knot[i].y = Veclib::splint (_geom->x.size(), s, &_geom->arclen[0],
				&_geom->y[0], &_geom->sy[0]);

  }
}


void Spline::printNek () const
// ---------------------------------------------------------------------------
// Print spline info in NEKTON format -- currently unknown.
// ---------------------------------------------------------------------------
{
  cout << "Unknown format for NEKTON splined edge -- fix me." << endl;
}


spline2D* Spline::getGeom (const char* fname)
// ---------------------------------------------------------------------------
// Return a natural spline curve of given name from internal store,
// otherwise open a file of the same name and compute coefficients
// based on the knots it contains.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "Spline::getGeom";
  const int_t verbose = Femlib::ivalue("VERBOSE");
  bool        found = false;
  spline2D*   c;
  char        err[StrMax];

  VERBOSE cerr << routine << ": spline filename: " << fname << endl;
 
  for (vector<spline2D*>::const_iterator k = gcurve.begin();
       !found && k != gcurve.end(); k++) {
    c = *k; found = (strcmp (c->name, fname) == 0);
  }

  if (!found) {

    VERBOSE cerr << routine << ": adding new spline curve" << endl;

    const int_t    SECTOR_MAX = 16;
    ifstream       file (fname);   
    real_t         x, y;
    int_t          i, j, N;
    vector<real_t> tmp;

    if (!file) {
      sprintf (err, "file: %s: not found", fname);
      message (routine, err, ERROR);
    }

    c = new spline2D;
    strcpy ((c->name = new char [strlen(fname) + 1]), fname);
    //    c->name = fname;
    c->pos  = 0;

    while (file.peek() == '#') file.ignore (StrMax, '\n');

    while (file >> x >> y) {
      c -> x.insert(c->x.end(), x);
      c -> y.insert(c->y.end(), y);
    }

    N = c->x.size();
    
    VERBOSE cerr << routine << ": " << N << " data points" << endl;

    c->sx.resize     (N);
    c->sy.resize     (N);
    c->arclen.resize (N);

    tmp.resize (N);
    
    Veclib::ramp   (N, 0, 1, &tmp[0], 1);
    Veclib::spline (N, FLT_MAX, FLT_MAX, &tmp[0], &c->x[0], &c->sx[0]);
    Veclib::spline (N, FLT_MAX, FLT_MAX, &tmp[0], &c->y[0], &c->sy[0]);

    for (c -> arclen[0] = 0.0, j = 0; j < N - 1; j++) {

      c -> arclen[j+1] = c -> arclen[j];

      for (i = 0; i < SECTOR_MAX; i++) {
	real_t p0x, p0y, p1x, p1y;

	if (i == 0) {
	  p0x = c->x[j];
	  p0y = c->y[j];
	  p1x = Veclib::splint (N, j+(i+1.0)/SECTOR_MAX,
			       &tmp[0],&c->x[0],&c->sx[0]);
	  p1y = Veclib::splint (N, j+(i+1.0)/SECTOR_MAX,
				&tmp[0], &c->y[0], &c->sy[0]);

	} else if (i == (SECTOR_MAX - 1)) {
	  p0x = Veclib::splint (N, j+(i*1.0)/SECTOR_MAX,
				&tmp[0], &c->x[0], &c->sx[0]);
	  p0y = Veclib::splint (N, j+(i*1.0)/SECTOR_MAX,
				&tmp[0], &c->y[0], &c->sy[0]);

	  p1x = c->x[j+1];
	  p1y = c->y[j+1];

	} else {
	  p0x = Veclib::splint (N, j+(i*1.0)/SECTOR_MAX,
				&tmp[0], &c->x[0], &c->sx[0]);
	  p0y = Veclib::splint (N, j+(i*1.0)/SECTOR_MAX,
				&tmp[0], &c->y[0], &c->sy[0]);
	  p1x = Veclib::splint (N, j+(i+1.0)/SECTOR_MAX,
				&tmp[0], &c->x[0], &c->sx[0]);
	  p1y = Veclib::splint (N, j+(i+1.0)/SECTOR_MAX,
				&tmp[0], &c->y[0], &c->sy[0]);

	}
	c -> arclen [j+1] += hypot (p1x-p0x, p1y-p0y);
      }
    }
    
    Veclib::spline (N, FLT_MAX, FLT_MAX, &c->arclen[0], &c->x[0], &c->sx[0]);
    Veclib::spline (N, FLT_MAX, FLT_MAX, &c->arclen[0], &c->y[0], &c->sy[0]);
    
    VERBOSE cerr << "arclength: " << c -> arclen[j] << endl;

    gcurve.insert (gcurve.end(), c);

    file.close();
  }

  return c;
}


int_t Spline::closest (const Point& p)
// ---------------------------------------------------------------------------
// Return the index of the knot point that lies closest to point gp. 
// Adjust for ends of interval.
// ---------------------------------------------------------------------------
{
  const int_t    M = gs->arclen.size();
  const int_t    N = M - gs->pos;
  const real_t   *x = &gs->x[0] + gs->pos, *y = &gs->y[0] + gs->pos;
  vector<real_t> len(N);
  int_t          i;

  for (i = 0; i < N; i++) len[i] = hypot (p.x - x[i], p.y - y[i]);

  i = Veclib::imin (N, &len[0], 1) + gs->pos;

  i = min (i, M - 2);

  if (i && i == gs -> pos) {
    gs -> pos = 0;
    i = closest (p);
  }
  return gs -> pos = i;

  return i;
}


static real_t getAngl (const real_t& s)
// ---------------------------------------------------------------------------
// This function computes an approximation (using the small angle
// formula) to the angle between the vector from (file-scope) control
// point "ga" through point "gp" and the vector from control point "ga"
// and the point on the splined curve that lies at arclength "s". It
// is this function that is minimised in order to estimate the
// position on the curve of point gp.
// ---------------------------------------------------------------------------
{
  real_t xs=Veclib::splint(gs->x.size(),s,&gs->arclen[0],&gs->x[0],&gs->sx[0]);
  real_t ys=Veclib::splint(gs->y.size(),s,&gs->arclen[0],&gs->y[0],&gs->sy[0]);

  real_t dxp = gp.x - ga.x;
  real_t dyp = gp.y - ga.y;

  real_t dxc = xs   - ga.x;
  real_t dyc = ys   - ga.y;

  return 1.0 - (dxp*dxc + dyp*dyc)/(hypot(dxp, dyp) * hypot(dxc, dyc));
}


real_t Spline::arcCoord ()
// ---------------------------------------------------------------------------
// Return the arclength that minimises the distance between gp and a
// point on the spline gs.
// ---------------------------------------------------------------------------
{
  const int_t i   = closest (gp);
  const int_t ip  = i + 1;
  real_t      TOL = 1.0e-6;
  real_t      s0, s1, s2;
  real_t      f0, f1, f2;

  s0 = gs->arclen[i];
  s1 = gs->arclen[ip];
#if 0  
  Recipes::mnbrak (s0, s1, s2, f0, f1, f2, ::getAngl);
  if (fabs (f1) > TOL) {
    Recipes::brent (s0, s1, s2, ::getAngl, TOL, f1);
    s1 = f1;
  }
#else
  Femlib::bracket (s0, s1, s2, f0, f1, f2, ::getAngl);
  if (fabs (f1) > TOL) {
    s1 = Femlib::brent (s0, s2, ::getAngl, TOL);
  }
#endif
  return s1;
}
