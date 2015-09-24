///////////////////////////////////////////////////////////////////////////////
// repmesh.cpp: read semtex session file: NODES and ELEMENTS, and
// generate planar reflection about either x or y axis, or generate
// rotation about named point.  Combine the result, deleting
// non-unique NODES.  Print on cout.
//
// Hand-editing will still be required to make valid SURFACES and
// CURVES sections.
//
// repmesh [-x || -y] [-r x0 y0 ang nrep] [-h] session
//
// Either X or Y reflections change the sense of rotation around
// elements.  Angular rotations are taken CCW, measured in degrees.
// Nrep is the integer number of times to repeat the rotation.
//
// Warning: This code is not robust to errors or irregularities in the
// input session file.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: repmesh.cpp,v 8.3 2015/09/08 21:43:00 hmb Exp $";

#include <sem.h>

class Node {
public:
  Node (const int_t& i, const Point& p) : id (i), loc (p) { }
  const int_t& ID  () const { return id;  }
  const Point& pos () const { return loc; }

  int_t id; Point loc; 
};

enum   Action { UNDEFINED, X, Y, R };
static char       prog[] = "repeatmesh";
static const real_t D2R    = 1.0 / 57.29577951308232087721;

static void  getArgs     (int,char**, Action&, Point&, real_t&, int_t&, FEML*&);
static int_t getVertices (FEML*, vector<Node*>&, const Action,
			  const Point&, const real_t, const int_t);
static int_t getElements (FEML*, vector<Node*>&, vector<Node*>*&, 
			  const int_t, const int_t);
static void  printUp     (ostream&, const vector<Node*>&,const vector<Node*>*&,
			  const int_t, const int_t, const int_t, const Action);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver for other routines.
// ---------------------------------------------------------------------------
{
  Action         mirror = UNDEFINED;
  Point          axis;
  real_t         angle = 0.;
  int_t          nvert, nel;	// -- Numbers read from input.
  int_t          nrep = 1;	// -- Number of times to repeat operation.
  vector<Node*>  vertices;
  vector<Node*>* elements;
  FEML*          F;

  Femlib::initialize (&argc, &argv);

  getArgs (argc, argv, mirror, axis, angle, nrep, F);

  nvert = getVertices (F, vertices, mirror, axis, angle, nrep);
  nel   = getElements (F, vertices, elements, nvert, nrep);

  printUp  (cout, const_cast<const vector<Node*>&>(vertices),
	    const_cast<const vector<Node*>*&>(elements),
	    nvert, nel, nrep, mirror);

  Femlib::finalize();
  
  return EXIT_SUCCESS;
}


static void getArgs (int     argc  ,
		     char**  argv  ,
		     Action& mirror,
		     Point&  axis  ,
		     real_t& angle ,
		     int_t&  nrep  ,
		     FEML*&  F     )
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last (optional) argument is the name of an input file.
// ---------------------------------------------------------------------------
{
  char buf[StrMax];
  char usage[] = "Usage: %s [-x || -y] [-r x0 y0 ang nrep] [-h] [file]\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'x':
      mirror = X;
      break;
    case 'y':
      mirror = Y;
      break;
    case 'r':
      mirror = R;
      --argc;
      axis.x = atof (*++argv);
      --argc;
      axis.y = atof (*++argv);
      --argc;
      angle  = atof (*++argv) * D2R;
      --argc;
      nrep   = atoi (*++argv);
      break;
    default:
      sprintf (buf, usage, prog);
      cerr << buf;
      exit (EXIT_FAILURE);
      break;
    }

  if (mirror == UNDEFINED) {
    sprintf (buf, usage, prog);
    cout << buf;
    exit (EXIT_FAILURE);
  }
        
  if (argc == 1) F = new FEML (*argv);
  else { sprintf (buf, usage, prog); cout << buf; exit (EXIT_FAILURE); }
}


static int_t getVertices (FEML*          F   ,
			  vector<Node*>& V   ,
			  const Action   A   ,
			  const Point&   O   ,
			  const real_t   angl,
			  const int_t    nrep)
// ---------------------------------------------------------------------------
// Read in declared vertices and create corresponding Nodes.
//
// Vertices for reflected Nodes are also generated in the upper locations
// of V, so there is a rule for getting from the original Nodes to the
// reflected Nodes and vice versa.  If the reflected Nodes land in the
// same locations (within rounding) as the originals, the Node pointers are
// the same as the unreflected ones.  Later on we have to construct
// a list of the unique Nodes to print out.
//
// Note that we check that input vertex IDs are given in order, since
// later the element node IDs can be used to index into V without clashes.
// ---------------------------------------------------------------------------
{
  char   routine[] = "getVertices", err[StrMax], buf[StrMax], *c;
  int_t  i, j, k, m, id, num, refId;
  bool   found;
  Node*  N;
  Point  P, Pd;
  real_t size, x, y, z, cc = cos (angl), ss = sin (angl);
  
  num = F -> attribute ("NODES", "NUMBER");
  refId = num; 
  
  V.resize (num + nrep*num);

  for (i = 0; i < num; i++) {    // -- Input Nodes (assumed unique).
    F -> stream() >> id >> x >> y >> z;
    if (id != i + 1)
      message (routine, "input list of NODES is out of order: check", ERROR);
    if (fabs(x) < EPSSP) x = 0.0;
    if (fabs(y) < EPSSP) y = 0.0;
    P.x = x; P.y = y; P.z = 0.0;
    V[i] = N = new Node (id, P);
  }

  // -- Generate reflected/rotated Nodes (only if geometrically unique).

  for (k = num, j = 0; j < nrep; j++) {

    cc = cos (angl * (j + 1));
    ss = sin (angl * (j + 1));

    for (i = 0; i < num; i++, k++) {

      x = V[i] -> pos().x;
      y = V[i] -> pos().y;

      switch (A) {
      case X: Pd.x = -x; Pd.y =  y; break;
      case Y: Pd.x =  x; Pd.y = -y; break;
      case R:
	Pd.x = x*cc - y*ss + O.x*(1.0-cc) + O.y*ss;
	Pd.y = x*ss + y*cc + O.y*(1.0-cc) - O.x*ss; break;
      }

      // -- Search over all nodes already installed.

      for (found = false, m = 0; !found && m < k; m++) {
	N = V[m];
	if (hypot (Pd.x-(N->pos().x), Pd.y-(N->pos().y)) < EPSm5) found = true;
      }
      if (!found) N = new Node (++refId, Pd);

      V[k] = N;
    }
  }

  return num;
}


static int_t getElements (FEML*           F    ,
			  vector<Node*>&  V    ,
			  vector<Node*>*& E    ,
			  const int_t     nvert,
			  const int_t     nrep )
// ---------------------------------------------------------------------------
// Read in declared elements, match up NODE IDs with Node*'s.
// Generate reflected elements.
//
// Note that we know how many elements there'll be (twice the number declared
// in input) but not the number of unique vertices, because some may retain
// their identity under reflection.
// ---------------------------------------------------------------------------
{
  char           routine[] = "getElements", buf[StrMax], *c;
  register int_t i, j, k;
  int_t          id, nel;

  nel = F -> attribute ("ELEMENTS", "NUMBER");

  E = new vector<Node*> [nel * (nrep+1)];
  
  for (i = 0; i < nel; i++) {
    F -> stream() >> id >> buf;
    E[i].resize (4);
    for (j = 0; j < 4; j++) {
      F -> stream() >> id;
      E[i][j] = V[id - 1];
    }
    F -> stream() . getline (buf, StrMax);
  }

  for (k = 1; k <= nrep; k++)
    for (i = 0; i < nel; i++) {
      E[i+ k*nel].resize (4);
      for (j = 0; j < 4; j++)
	E[i + k*nel][j] = V[ (E[i][j]->ID() - 1) + k*nvert];
    }

  return nel;
}


static void printUp (ostream&              S  ,
		     const vector<Node*>&  V  ,
		     const vector<Node*>*& E  ,
		     const int_t           nvert,
		     const int_t           nel,
		     const int_t           nrep,
		     const Action          A  )
// ---------------------------------------------------------------------------
// The output format is semtex NODE/ELEMENT information (same as it
// received on input).
// ---------------------------------------------------------------------------
{
  list<Node*>           U;
  list<Node*>::iterator p;
  register int_t        i, j;
  Node*                 n;
  int_t                 J;
  bool                  found;

  // -- Create & print up list of unique Nodes.

  for (i = 0; i < nvert * (nrep + 1); i++) {
    for (found = false, p = U.begin(); !found && p != U.end(); p++)
      found = *p == V[i];
    if (!found) U.push_back (V[i]);
  }

  S << "<NODES NUMBER=" << U.size() << ">" << endl;
  for (p = U.begin(); p != U.end(); p++) {
    n = (*p);
    S << n->ID() << '\t' << n->pos().x << '\t' << n->pos().y << "\t0" << endl;
  }

  S << "</NODES>" << endl << endl;
  
  // -- Print the original elements.

  S << "<ELEMENTS NUMBER=" << nel * (nrep + 1) << ">" << endl;
  for (i = 0; i < nel; i++) {
    S << i + 1 << "\t<Q>"; 
    J = E[i].size();
    for (j = 0; j < J; j++) S << setw(5) << E[i][j] -> ID();
    S << "\t</Q>" << endl;
  }
  
  // -- Print up reflected elements.

  for (i = 0; i < nrep * nel; i++) {
    S << nel + i + 1 << "\t<Q>";
    J = E[i + nel].size();
    for (j = 0; j < J; j++)
      S << setw(5)
	<< ((A == R) ? E[i + nel][j] -> ID() : E[i + nel][J - j - 1] -> ID());
    S << "\t</Q>" << endl;
  }
  S << "</ELEMENTS>" << endl << endl;

#if 0
  // -- Print up boundary Nodes (this will need hand editing).

  S << endl << "2  Boundaries";

  S << endl;
  S << "1  1  ";
  j = 0;
  for (i = 0; i < N; i++) j += (V[i] -> interior()) ? 0 : 1;
  S << setw (5) << j;
  for (i = 0; i < N; i++)
    if (!V[i] -> interior()) S << setw (5) << V[i] -> ID();

  S << endl;
  S << "2  1  ";
  j = 0;
  for (i = N; i < N + N; i++) j += (V[i] -> interior()) ? 0 : 1;
  S << setw (5) << j;
  for (i = N + N - 1; i >= N; i--)
    if (!V[i] -> interior()) S << setw (5) << V[i] -> ID();  

  S << endl;
  S << endl << "0  Curves" << endl;
  
  S << "}" << endl;
#endif
}

