#ifndef MESH_H
#define MESH_H
///////////////////////////////////////////////////////////////////////////////
// mesh: header file for Mesh and related classes.
//
// $Id: mesh.h,v 8.1 2015/04/20 11:14:18 hmb Exp $
///////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>

#include "cfemdef.h"
#include "feml.h"

class Curve;


class Mesh
// ===========================================================================
// Mesh class is used as a storage for inter-element connectivities,
// BC tags, element-corner vertex locations and curved boundaries.
// It plays no direct part in the eventual solution of the discrete
// equations.  Presently, meshes are 2D only.
// 
// Mesh provides seven externally-visible facilities: 
// (1) Return number of elements in mesh,
// (2) Generation of element-boundary boundary-to-global mapping vector.
// (3) Generation of element-boundary Essential BC mask vector.
// (4) Generation of element mesh knot points.
// (5) Print up NEKTON-style .rea information (backwards compatibility).
// (6) Return x--y extent of mesh, based on traverse of element vertices.
// (7) Return true if an element touches the z axis (cylindrical coords).
//
// (2--4) can be independently computed for arbitrary polynomial orders.
// ===========================================================================
{
public:
  class Node;
  class Elmt;
  class Side;

  enum IDstatus {UNSET = -1};

  Mesh (FEML*, const bool = true);

  int_t nEl         () const { return _elmtTable.size(); }
  int_t buildMap    (const int_t, int_t*);
  void  buildMask   (const int_t, const char, int_t*);
  void  meshElmt    (const int_t, const int_t, const real_t*, const real_t*,
		     real_t*, real_t*) const;
  void  printNek    () const;
  void  extent      (Point&, Point&) const;
  void  assemble    (const bool = false);

  static void showGlobalID (Mesh&);
  static void showAssembly (Mesh&);

  class Node {
  public:
    int_t ID;
    int_t gID;
    Point loc;
    Node* periodic;
  };

  class Elmt {
    friend class Mesh;
  private:
    vector<Node*> node;
    vector<Side*> side;
  public:
    int_t ID;
    int_t nNodes   () const  { return node.size(); }
    Node* ccwNode  (int_t i) { return node[(i         +1) % nNodes()]; }
    Node* cwNode   (int_t i) { return node[(i+nNodes()-1) % nNodes()]; }
    Point centroid () const;
  };
  
  class Side {
    friend class Mesh;
  private:
    Side (int_t id, Node* n1, Node* n2) :
      ID(id), startNode(n1), endNode(n2), mateElmt(0) {
      gID.resize (0); mateSide = 0;
      axial = (n1 -> loc.y == 0.0) && (n2 -> loc.y == 0.0);
    }
    Side () {}
  public:
    int_t         ID;
    bool          axial;
    vector<int_t> gID;
    Node*         startNode;
    Node*         endNode;
    Elmt*         thisElmt;
    Elmt*         mateElmt;	// -- Doubles as a flag for union:
    union { Side* mateSide; char group; };
    void connect (const int_t, int_t&);
  };

private:
  FEML&          _feml;
  vector<Node*>  _nodeTable;
  vector<Elmt*>  _elmtTable;
  vector<Curve*> _curveTable;

  void surfaces      ();
  void curves        ();

  void checkAssembly ();
  void chooseNode    (Node*, Node*);
  void fixPeriodic   ();

  void describeGrp (char, char*)          const;
  void describeBC  (char, char, char*)    const;
  bool matchBC     (const char, const char, const char);

  void meshSide (const int_t, const int_t, const int_t,
		 const real_t*, Point*)                          const;
};


class Curve
// ===========================================================================
// Base class for curved edges.
// ===========================================================================
{
public:
  virtual ~Curve () { }

  virtual void compute  (const int_t, const real_t*, Point*) const = 0;
  virtual void printNek ()                                   const = 0;

  bool ismatch (const int_t e, const int_t s) const {
    return (e == curveSide -> thisElmt -> ID && s == curveSide -> ID);
  }

protected:
  Mesh::Side* curveSide;
};


class CircularArc : public Curve
// ===========================================================================
// Prototype curved edge class.   To make new kinds, retain the same
// interface structure and put new instances in Mesh::curves.
// ===========================================================================
{
public:
  CircularArc (const int_t, Mesh::Side*, const real_t);
  virtual void compute  (const int_t, const real_t*, Point*) const;
  virtual void printNek () const;

private:
  Point  centre;
  real_t radius;
  real_t semiangle;
  int_t  convexity;
};


class spline2D
// ===========================================================================
// Definition of spline2D struct used by Spline class.
// ===========================================================================
{
public:
  char*          name  ;	// name of file containing knot points
  vector<real_t> x     ;	// knot x-coordinates
  vector<real_t> y     ;	// knot y-coordinates
  vector<real_t> sx    ;	// spline x-coefficients
  vector<real_t> sy    ;	// spline y-coefficients
  vector<real_t> arclen;        // arclength along curve at knot locations
  int_t          pos   ;	// index of last confirmed position
};


class Spline : public Curve
// ===========================================================================
// A 2D splined curve defined in a named file.
// ===========================================================================
{
public:
  Spline (const int_t, Mesh::Side*, const char*);
  virtual void compute  (const int_t, const real_t*, Point*) const;
  virtual void printNek () const;

private:
  char*     _name;
  spline2D* _geom;
  real_t    _startarc;
  real_t    _endarc;

  spline2D* getGeom  (const char*);
  real_t    getAngl  (const real_t&);
  int_t     closest  (const Point&);
  real_t    arcCoord ();
};


#endif
