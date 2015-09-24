//////////////////////////////////////////////////////////////////////////////
// enumerate.C:  utility to generate mesh numbering from session file.
//
// Copyright (c) 1995 <--> $Date: 2015/04/20 11:14:18 $, Hugh Blackburn
//
// Usage: enumerate [options] file
//   options:
//   -h       ... display this message
//   -v       ... set verbose output
//   -n N     ... override element order to be N
//   -O [0-3] ... set level of bandwidth optimization
//
// Special action may need to be taken to generate numbering schemes
// for cylindrical coordinate flow problems.  See the discussion in
// header for field.C, and for routine Mesh::buildMask in mesh.C.
//
// Divergence problems sometimes arise when the highest numbered
// zero-mode pressure node occurs on the axis in 3D cylindrical
// simulations.  The code attempts to fix this problem, and lets you
// know if it can't.
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
// 02110-1301 USA
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: enumerate.cpp,v 8.1 2015/04/20 11:14:18 hmb Exp $";

#include <cstdlib>
#include <cstring>
#include <climits>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

#include <cfemdef.h>
#include <femlib.h>
#include <utility.h>
#include <mesh.h>
#include <veclib.h>

class Nsys {
friend void printup (const char*, vector<Nsys*>&, const int_t);

public:
  Nsys (char, vector<int_t>&, vector<int_t>&, const int_t, 
	const int_t, const int_t, const int_t,
	const int_t, const int_t, const int_t);

  bool match    (vector<int_t>&);
  void addField (char);

  int_t         nel;
  int_t         np_max;
  int_t         next_max;
  int_t         nint_max;
  int_t         ntotal;
  int_t         nbndry;
  int_t         nglobal;
  int_t         nsolve;
  int_t         nbandw;
  int_t         optlev;
  vector<char>  fields;
  vector<int_t> bndmap;
  vector<int_t> bndmsk;
  vector<int_t> axelmt;
  vector<int_t> axside;

  int_t sortGid         (int_t*, int_t*);
  void  renumber        (const int_t, const int_t = 0);
  int_t buildAdjncy     (vector<int_t>*)                              const;
  void  fillAdjncy      (vector<int_t>*, int_t*, int_t*, const int_t) const;
  void  connectivSC     (vector<int_t>*, const int_t*,
			 const int_t*, const int_t)                   const;
  int_t globalBandwidth ()                                            const;
  bool  highAxis        ()                                            const;
  int_t bandwidthSC     (const int_t*, const int_t*, const int_t)     const;
  void  rebuild         (FEML*, const int_t);

};

static char        prog[] = "enumerate";
static int_t       verb   = 0;
static const int_t FldMax = 16;

static void getargs   (int, char**, char*&, int_t&, int_t&, int_t&);
static char axial     (FEML*);
static void getfields (FEML*, char*, const bool);
static void checkVBCs (FEML*, const char*);
static void checkABCs (FEML*, const char);
       void printup   (const char*, vector<Nsys*>&, const int_t);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Determine, from BCs section of FEML file, list of fields for which
// numbering schemes are to be constructed.
//
// Generate a BC mask and initial numbering scheme for first field, using
// Mesh class routines.  Optimize numbering scheme according to selected level.
//
// For each succeeding field, first generate a BC mask and, if it matches
// a mask previously generated, add the field's name to the previous field's
// name vector but take no further action.  Otherwise, generate and optimize
// a new numbering system.
//
// Print up the masks and numbering schemes on cout.
// ---------------------------------------------------------------------------
{
  char   *session = 0, field[StrMax], axistag;
  FEML*  file;
  int_t  np = 0, opt = 1;
  bool   cyl3D = false;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, session, verb, np, opt);

  file = new FEML (session);

  if (verb)  Femlib::ivalue ("VERBOSE", verb);
  if   (np)  Femlib::ivalue ("N_P", np);
  else  np = Femlib::ivalue ("N_P");

  cyl3D = static_cast<bool>(Femlib::ivalue ("CYLINDRICAL"));

  getfields (file, field, (axistag = axial (file)) && cyl3D);
  if (axistag) checkABCs (file, axistag);
  if (cyl3D)   checkVBCs (file, field);
  
  Mesh          M (file);
  vector<Nsys*> S (strlen (field));

  int_t         i, j, k = 0;
  bool          found;
  const int_t   NEL      = M.nEl();
  const int_t   NP_MAX   = np;
  const int_t   NEXT_MAX = 4 * (NP_MAX - 1);
  const int_t   NINT_MAX = sqr (NP_MAX - 2);
  const int_t   NBNDRY   = NEL * NEXT_MAX;
  const int_t   NTOTAL   = NEL * sqr (NP_MAX);

  vector<int_t> btog (NBNDRY);
  vector<int_t> mask (NBNDRY);

  M.buildMask (np, field[0], &mask[0]);
  M.buildMap  (np, &btog[0]);

  S[k++] = new Nsys (field[0], btog, mask, opt,
		     NEL, NTOTAL, NBNDRY, NP_MAX, NEXT_MAX, NINT_MAX);

  for (i = 1; i < strlen (field); i++) {
    M.buildMask (np, field[i], &mask[0]);
    found = false;
    for (j = 0; !found && j < k; j++)
      if (found = S[j] -> match (mask)) 
	S[j] -> addField (field[i]);
    if (!found) {
      M.buildMap (np, &btog[0]);
      S[k++] = new Nsys (field[i], btog, mask, opt,
			 NEL, NTOTAL, NBNDRY, NP_MAX, NEXT_MAX, NINT_MAX);
    }
  }

  // -- Potentially have to fix numbering for cylindrical mode-zero pressure.

  if (axistag) {
    for (i = 0; i < k; i++)
      if (strchr (&S[i] -> fields[0], 'p')) {
	Nsys* pressure = S[i];
	if (!Veclib::any (pressure -> nbndry, &pressure -> bndmsk[0], 1)) {
	  pressure -> rebuild (file, clamp (static_cast<int>(opt), 2, 3));
	  break;
	}
      }
  }

  printup (field, S, k);

  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int    argc   , 
		     char** argv   ,
		     char*& session,
		     int_t& verb   ,
		     int_t& np     ,
		     int_t& opt    )
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "usage: enumerate [options] session\n"
                 "options:\n"
                 "  -h       ... display this message\n"
                 "  -v       ... set verbose output\n"
		 "  -n N     ... override number of element knots to be N\n"
		 "  -O [0-3] ... bandwidth optimization level [Default: 1]\n";
  char err[StrMax];
  char c;

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      cerr << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      for (verb = 1; *++argv[0] == 'v'; verb++);
      break;
    case 'n':
      if (*++argv[0])
	np = atoi (*argv);
      else {
	--argc;
	np = atoi (*++argv);
      }
      break;
    case 'O':
      if (*++argv[0])
	opt = atoi (*argv);
      else {
	--argc;
	opt = atoi (*++argv);
      }
      break;
    default:
      sprintf (err, "getargs: illegal option: %c\n", c);
      message (prog, err, ERROR);
      break;
    }

  if   (argc == 1) session = *argv;
  else             message (prog, "must provide session file", ERROR);
}


static char axial (FEML* file)
// ---------------------------------------------------------------------------
// Return the character tag corresponding to "axis" group, if that  exists.
// ---------------------------------------------------------------------------
{
  if (!file->seek ("GROUPS")) return '\0';

  int_t       i;
  char        nextc, buf[StrMax];
  const int_t N (file->attribute ("GROUPS", "NUMBER"));
  
  for (i = 0; i < N; i++) {
    while ((nextc = file->stream().peek()) == '#') // -- Skip comments.
      file->stream().ignore (StrMax, '\n');
    file->stream() >> buf >> nextc >> buf;
    if (strstr (buf, "axis")) return nextc;
  }

  return '\0';
}


static void getfields (FEML*      file     ,
		       char*      field    ,
		       const bool axisModes)
// ---------------------------------------------------------------------------
// Set up the list of fields according to the names found in the
// 'FIELDS' section of the FEML file, e.g.
//
// <FIELDS>
// # optional comment lines...
//   c u v w p
// </FIELDS>
//
// Input value "axisModes" flags that this is a problem with an axis boundary
// condition, cylindrical coordinates, and three space dimensions.
// In that case, there will be of extra numbering schemes set up for the
// non-zero Fourier modes for variables c, u, w & p (if they exist).
// The extra schemes get names C, U, W, & P.
// ---------------------------------------------------------------------------
{
  int_t i = 0;
  char  t[StrMax];

  // -- Set up string for the fields listed.

  if (file->seek ("FIELDS")) {

    file->stream().ignore (StrMax, '\n');
    
    while (file->stream().peek() == '#') // -- Skip comments.
      file->stream().ignore (StrMax, '\n');

    do {
      file->stream() >> field[i++];
    } while (field[i - 1] != '<' && i < StrMax);

    if (field[--i] == '<') {
      field[i] = '\0';
      file->stream() >> t;
      if (!(strstr (t,   "/FIELDS")))
	   message (prog, "FIELDS section not closed", ERROR);
    } else message (prog, "FIELDS section not closed", ERROR);
  } else   message (prog, "FIELDS section not found",  ERROR);

  // -- Add extra names for higher Fourier modes if required.

  if (axisModes) {
    if (strchr (field, 'c')) {
      field[i++] = 'C';
      field[i]   = '\0';
    }
    if (strchr (field, 'u')) {
      field[i++] = 'U';
      field[i]   = '\0';
    }
    if (strchr (field, 'w')) {
      field[i++] = 'W';
      field[i]   = '\0';
    }
    if (strchr (field, 'p')) {
      field[i++] = 'P';
      field[i]   = '\0';
    }
  }
}


static void checkVBCs (FEML*       file ,
		       const char* field)
// ---------------------------------------------------------------------------
// For cylindrical 3D fluids problems, the declared boundary condition types
// for velocity fields v & w must be the same for all groups, to allow
// for coupling of these fields (which uncouples the viscous substep).
//
// Check it out by running through each group's BCs and checking for 
// tag agreement on v & w BCs.
// ---------------------------------------------------------------------------
{
  if (!file->seek ("BCS")) return;
  if (!strchr (field, 'u')) return;
  if (!strchr (field, 'v') || !strchr (field, 'w')) return;

  int_t       i, j, id, nbcs;
  char          vtag, wtag, groupc, fieldc, tagc, tag[StrMax], err[StrMax];
  const int_t N (file->attribute ("BCS", "NUMBER"));

  for (i = 0; i < N; i++) {

    while ((groupc = file->stream().peek()) == '#') // -- Skip comments.
      file->stream().ignore (StrMax, '\n');

    file->stream() >> id >> groupc >> nbcs;
    vtag = wtag = '\0';

    for (j = 0; j < nbcs; j++) {

      file->stream() >> tag;
      if (strchr (tag, '<') && strchr (tag, '>') && (strlen (tag) == 3))
	tagc = tag[1];
      else {
	sprintf (err, "unrecognized BC tag format: %s", tag);
	message (prog, err, ERROR);
      }

      file->stream() >> fieldc;
      if      (fieldc == 'v') vtag = tagc;
      else if (fieldc == 'w') wtag = tagc;
      file->stream().ignore (StrMax, '\n');
    }
    
    if (!(vtag && wtag)) {
      sprintf (err, "group %c: BCs for fields 'v' & 'w' needed", groupc);
      message (prog, err, ERROR);
    }
    if (vtag != wtag) {
      sprintf (err, "group %c, fields 'v' & 'w': BC type mismatch", groupc);
      message (prog, err, ERROR);
    }
  }
}


static void checkABCs (FEML*      file ,
		       const char atag)
// ---------------------------------------------------------------------------
// Run through and ensure that for "axis" group, all BCs are of type <A>.
// ---------------------------------------------------------------------------
{
  if (!file->seek ("BCS")) return;

  int_t       i, j, id, nbcs;
  char        groupc, fieldc, tagc, tag[StrMax], err[StrMax];
  const int_t N (file->attribute ("BCS", "NUMBER"));

  for (i = 0; i < N; i++) {

    while ((groupc = file->stream().peek()) == '#') // -- Skip comments.
      file->stream().ignore (StrMax, '\n');

    file->stream() >> id >> groupc >> nbcs;

    for (j = 0; j < nbcs; j++) {

      file->stream() >> tag;
      if (strchr (tag, '<') && strchr (tag, '>') && (strlen (tag) == 3))
	tagc = tag[1];
      else {
	sprintf (err, "unrecognized BC tag format: %s", tag);
	message (prog, err, ERROR);
      }

      file->stream() >> fieldc;

      if (groupc == atag && tagc != 'A') {
	sprintf (err, "group '%c': field '%c' needs axis BC", groupc, fieldc);
	message (prog, err, ERROR);
      }
      file->stream().ignore (StrMax, '\n');
    }
  }
}


void printup (const char*    F   ,
	      vector<Nsys*>& S   ,
	      const int_t    nSys)
// ---------------------------------------------------------------------------
// Print up summary info followed by map & mask for each system.
// ---------------------------------------------------------------------------
{
  register int_t i, j, k, side, soff;
  const    int_t nedge = S[0] -> nbndry / (4 * S[0] -> nel);
  
  cout << "# FIELDS         :  " << F << endl;

  cout << "# ----------------";
  for (j = 0; j < nSys; j++) cout << "  -------------";
  cout << endl;

  cout << "# " << nSys << " NUMBER SETS  :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << &S[j] -> fields[0];
  }
  cout << endl;

  cout << "# NEL            :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << S[j] -> nel;
  }
  cout << endl;
  
  cout << "# NP_MAX         :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << S[j] -> np_max;
  }
  cout << endl;
  
  cout << "# NEXT_MAX       :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << S[j] -> next_max;
  }
  cout << endl;
  
  cout << "# NINT_MAX       :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << S[j] -> nint_max;
  }
  cout << endl;
  
  cout << "# NTOTAL         :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << S[j] -> ntotal;
  }
  cout << endl;
  
  cout << "# NBOUNDARY      :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << S[j] -> nbndry;
  }
  cout << endl;

  cout << "# NGLOBAL        :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << S[j] -> nglobal;
  }
  cout << endl;

  cout << "# NSOLVE         :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << S[j] -> nsolve;
  }
  cout << endl;

  cout << "# OPTIMIZATION   :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << S[j] -> optlev;
  }
  cout << endl;

  cout << "# BANDWIDTH      :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << S[j] -> nbandw;
  }
  cout << endl;

  cout << "# ----------------";
  for (j = 0; j < nSys; j++) cout << "  -------------";
  cout << endl;

  cout << "# elmt  side offst";
  for (j = 0; j < nSys; j++) cout << "     bmap  mask";
  cout << endl;

  for (i = 0, k = 1; k <= S[0] -> nel; k++)
    for (side = 1; side <= 4; side++)
      for (soff = 0; soff < nedge; soff++, i++) {
	cout << setw (6) << k << setw (6) << side << setw (6) << soff;
	for (j = 0; j < nSys; j++)
	  cout 
	    << setw (9) << S[j] -> bndmap[i]
	      << setw (6) << S[j] -> bndmsk[i];
	cout << endl;
      }
}


Nsys::Nsys (char           name    ,
	    vector<int_t>& map     ,
	    vector<int_t>& mask    ,
	    const  int_t   opt     ,
	    const  int_t   NEL     ,
	    const  int_t   NTOTAL  ,
	    const  int_t   NBNDRY  ,
	    const  int_t   NP_MAX  ,
	    const  int_t   NEXT_MAX,
	    const  int_t   NINT_MAX) :

	    nel           (NEL     ),
	    np_max        (NP_MAX  ),
	    next_max      (NEXT_MAX),
	    nint_max      (NINT_MAX),
	    ntotal        (NTOTAL  ),
	    nbndry        (NBNDRY  ),
	    optlev        (opt     )
// ---------------------------------------------------------------------------
// Constructor also carries out bandwidth optimization task.
// ---------------------------------------------------------------------------
{
  fields.resize (FldMax);
  memset (&fields[0], '\0', FldMax);
  fields[0] = name;
  bndmap    = map;
  bndmsk    = mask;
  nbndry    = bndmap.size();
  nglobal   = bndmap [Veclib::imax (nbndry, &bndmap[0], 1)] + 1;
  nsolve    = sortGid (&bndmap[0], &bndmsk[0]);

  renumber (optlev);
  nbandw = globalBandwidth ();
}


bool Nsys::match (vector<int_t>& test)
// ---------------------------------------------------------------------------
// Return true if test and bndmsk match.
// ---------------------------------------------------------------------------
{
  if (test.size() != bndmsk.size()) return false;

  return
    static_cast<bool>(Veclib::same (bndmsk.size(), &bndmsk[0],1, &test[0],1));
}


void Nsys::addField (char name)
// ---------------------------------------------------------------------------
// Add a new field which matches the present one.
// ---------------------------------------------------------------------------
{
  int_t k = 0;

  while (fields[k]) k++;

  if   (k == FldMax) message (prog, "too many fields", ERROR);
  else               fields  [k] = name;
}


static int cmp1 (const void *a,
		 const void *b)
// ---------------------------------------------------------------------------
// Used by qsort.  Compare first element (global node number) of two arrays.
// ---------------------------------------------------------------------------
{ return static_cast<int>
    (static_cast<const int_t *>(a)[0]-static_cast<const int_t *>(b)[0]); }


static int cmp2 (const void *a,
		 const void *b)
// ---------------------------------------------------------------------------
// Used by qsort.  Compare second element (solve mask) of two arrays.
// ---------------------------------------------------------------------------
{ return static_cast<int>
    (static_cast<const int_t *>(a)[1]-static_cast<const int_t *>(b)[1]); }



int_t Nsys::sortGid (int_t* bmap,
		     int_t* bmsk)
// ---------------------------------------------------------------------------
// Global node numbers get sorted to place essential-BC nodes last:
// this simplifies the later partition of global matrices.
//
// The non-essential type node numbers can be further sorted to
// optimize global matrix bandwidths, but this is not done here.
//
// A globally-numbered table (reOrder) is constructed, each entry of
// which is two ints: the global node number and the essential BC mask
// value (0/1).  This is then partitioned into "unknown" and "known"
// node numbers, with the "known" numbers last in the table.  Each
// partition is then sorted into ascending node number order, and the
// information is used to create a new element boundary-to-global
// numbering scheme in btog.
//
// Return number of global nodes at which solution is not set by
// essential BCs.
// ---------------------------------------------------------------------------
{
  vector<int_t> work (nbndry + 3 * nglobal);
  int_t         *bsave, *tmp, *reOrder;
  int_t         unknowns;

  bsave   = &work[0];
  tmp     = bsave + nbndry;
  reOrder = tmp  + nglobal;

  Veclib::copy  (nbndry,  bmap, 1, bsave, 1);
  Veclib::scatr (nbndry , bmsk, bsave, tmp);
  Veclib::ramp  (nglobal, 0,    1, reOrder,     2);
  Veclib::copy  (nglobal, tmp,  1, reOrder + 1, 2);
  
  unknowns = nglobal - Veclib::count (nglobal, tmp, 1);

  if (unknowns < nglobal) {

    // -- Partition into "unknown" nodes & "essential BC" nodes.

    qsort (reOrder,              nglobal,            2*sizeof (int_t), cmp2);

    // -- Sort each partition into ascending node number order.

    qsort (reOrder,              unknowns,           2*sizeof (int_t), cmp1);
    qsort (reOrder + 2*unknowns, nglobal - unknowns, 2*sizeof (int_t), cmp1);

    // -- Reload new gids.

    Veclib::copy  (nglobal, reOrder, 2,  tmp, 1);
    Veclib::ramp  (nglobal, 0,    1, reOrder, 1);
    Veclib::scatr (nglobal, reOrder, tmp, reOrder + nglobal);
    Veclib::gathr (nbndry , reOrder + nglobal, bsave, bmap);
  }

  return unknowns;
}


void Nsys::renumber (const int_t optlevel,
		     const int_t penalty )
// ---------------------------------------------------------------------------
// From the initial ordering specified in bndmap, use RCM to generate
// a reduced-bandwidth numbering scheme.  Reload into bndmap.
//
// Different optimization levels are allowed:
//
// 0: Do nothing (no renumbering).
// 1: Use FNROOT (trial root = 1) to find a pseudo-peripheral root node,
//    pass result to RCM for Reverse Cuthill McKee reordering.  Default level.
// 2: Use FNROOT to generate pseudo-peripheral nodes, but with trial roots
//    in steps of 10, up to n_solve.  Choose root to minimize global bandwidth.
// 3: Do not use FNROOT.  Try all unknown node numbers as trial roots for RCM.
//    Choose root to minimize global bandwidth.
//
// If penalty is set, then we apply a penalty to numbering schemes
// that have the highest-numbered pressure node on the axis.  By
// default penalty is off.
//
// Reference:
//    A. George and J. W-H. Liu
//    Computer Solution of Large Sparse Positive Definite Systems
//    Prentice-Hall (1981)
// ---------------------------------------------------------------------------
{
  if (!optlevel || !nsolve) return;

  if (verb)
    cout << "Bandwidth optimization (" << optlevel
	 << "), Field '" << &fields[0] << "'";

  register int_t i;
  int_t          root, nlvl;

  // -- Build node adjacency tables.
  
  vector<int_t>* adjncyList = new vector<int_t> [nsolve];
  const int_t    tabSize    = buildAdjncy (adjncyList);

  // -- Allocate memory.

  vector<int_t> work(tabSize + 1 + 4 * nsolve + 1 + nglobal + nbndry);
  int_t         *adjncy, *xadj, *perm, *mask, *xls, *invperm, *bsave;

  adjncy  = &work[0];
  xadj    = adjncy  + tabSize + 1;
  perm    = xadj    + nsolve  + 1;
  mask    = perm    + nsolve;
  xls     = mask    + nsolve;
  invperm = xls     + nsolve;
  bsave   = invperm + nglobal;
  
  Veclib::copy (nbndry, &bndmap[0], 1, bsave, 1);
  for (i = nsolve; i < nglobal; i++) invperm[i] = i;

  fillAdjncy (adjncyList, adjncy, xadj, tabSize);
  delete   [] adjncyList;

  switch (optlevel) {
  case 1: {
    root = 1;
    Veclib::fill   (nsolve, 1, mask, 1);
    Femlib::fnroot (root, xadj, adjncy, mask, nlvl, xls, perm);
    Femlib::rcm    (root, xadj, adjncy, mask, perm, nlvl, xls);
    break;
  }
  case 2: {
    int_t rtest, BWtest, BWmin = INT_MAX, best;

    if (verb) cout << ":" << endl;

    for (root = 1; root <= nsolve; root += 10) {
      rtest = root;
      Veclib::fill   (nsolve, 1, mask, 1);
      Femlib::fnroot (rtest, xadj, adjncy, mask, nlvl, xls, perm);
      Femlib::rcm    (rtest, xadj, adjncy, mask, perm, nlvl, xls);

      Veclib::sadd (nsolve, -1, perm, 1, perm, 1);
      for (i = 0; i < nsolve; i++) invperm[perm[i]] = i;
      Veclib::gathr (nbndry, invperm, bsave, &bndmap[0]);

      BWtest = globalBandwidth();
      if (penalty && highAxis()) BWtest += nglobal;
      if (BWtest < BWmin) {
	BWmin = BWtest;
	best  = rtest;
	if (verb) cout << "root = " << root << ", BW = " << BWmin << endl;
      }
    }

    Veclib::fill (nsolve, 1, mask, 1);
    Femlib::rcm  (best, xadj, adjncy, mask, perm, nlvl, xls );

    break;
  }
  case 3: {
    int_t BWtest, BWmin = INT_MAX, best;

    if (verb) cout << ":" << endl;

    for (root = 1; root <= nsolve; root++) {
      Veclib::fill (nsolve, 1, mask, 1);
      Femlib::rcm  (root, xadj, adjncy, mask, perm, nlvl, xls);

      Veclib::sadd (nsolve, -1, perm, 1, perm, 1);
      for (i = 0; i < nsolve; i++) invperm[perm[i]] = i;
      Veclib::gathr (nbndry, invperm, bsave, &bndmap[0]);

      BWtest = globalBandwidth();
      if (penalty && highAxis()) BWtest += nglobal;
      if (BWtest < BWmin) {
	BWmin = BWtest;
	best  = root;
	if (verb) cout << "root = " << root << ", BW = " << BWmin << endl;
      }
    }

    Veclib::fill (nsolve, 1, mask, 1);
    Femlib::rcm  (best, xadj, adjncy, mask, perm, nlvl, xls );

    break;
  }
  default:
    break;
  }

  Veclib::sadd (nsolve, -1, perm, 1, perm, 1);
  for (i = 0; i < nsolve; i++) invperm[perm[i]] = i;
  Veclib::gathr (nbndry, invperm, bsave, &bndmap[0]);

  if (penalty && highAxis())
    message (prog, "Highest numbered pressure node remains on axis", WARNING);

  if (verb) cout << endl;
}


int_t Nsys::buildAdjncy (vector<int_t>* adjncyList) const
// ---------------------------------------------------------------------------
// Traverse elements and build up a vector of linked lists that
// describe the global nodes adjacent to each global node.
//
// Return the total amount of storage required when information is packed
// into an int_t vector, as required by genrcm.
// ---------------------------------------------------------------------------
{
  register int_t k, ntab;
  const int_t    next = nbndry / nel;

  for (k = 0, ntab = 0; k < nel; k++) {
    connectivSC (adjncyList, &bndmap[0] + ntab, &bndmsk[0] + ntab, next);
    ntab += next;
  }

  for (k = 0, ntab = 0; k < nsolve; k++)
    ntab += adjncyList[k].size();

  return ntab;
}


void Nsys::fillAdjncy (vector<int_t>* adjncyList,
		       int_t*         adjncy    ,
		       int_t*         xadj      ,
		       const int_t    tabSize   ) const
// ---------------------------------------------------------------------------
// Load the information contained in adjncyList into the two vectors
// adjncy & xadj required by genrcm.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Nsys::fillAdjncy";
  register int_t        i, k;
  vector<int_t>::iterator p;

  for (i = 0, k = 1; i < nsolve; i++) {
    xadj[i] = k;
    for (p = adjncyList[i].begin(); p != adjncyList[i].end(); p++, k++)
      adjncy[k - 1] = *p + 1;
  }
  
  if (k != tabSize + 1)
    message (routine, "after traversing list, k != tabSize + 1", ERROR);

  adjncy[tabSize] = 0;
  xadj  [ nsolve] = k;
}


void Nsys::connectivSC (vector<int_t>* adjList,
			const int_t*   bmap   ,
			const int_t*   mask   ,
			const int_t    next   ) const
// ---------------------------------------------------------------------------
// AdjList is an array of linked lists, each of which describes the
// global nodes that have connectivity with the the current node,
// i.e. which make a contribution to the weighted-residual integral
// for this node.  This routine fills in the contribution from the
// current element.
//
// For general finite elements, all nodes of an element are
// interconnected, while for statically-condensed elements, only the
// boundary nodes are considered (since internal nodes are not
// global).
//
// Essential-BC nodes are ignored, since we're only interested in
// mimimizing bandwidths of global matrices.
// ---------------------------------------------------------------------------
{
  register int_t          i, j, found, gidCurr, gidMate;
  vector<int_t>::iterator a;
  
  for (i = 0; i < next; i++) {
    if (! mask[i]) {
      gidCurr = bmap[i];
      
      for (j = 0; j < next; j++) {
	if (i != j && ! mask[j]) {
	  for (a = adjList[gidCurr].begin(), gidMate = bmap[j], found = 0;
	       !found && a != adjList[gidCurr].end(); a++)
	    found = *a == gidMate;
	  if (!found) adjList[gidCurr].push_back (gidMate);
	}
      }
    }
  }
}


int_t Nsys::globalBandwidth () const
// --------------------------------------------------------------------------
// Precompute the bandwidth of assembled global matrix (including
// diagonal).
// --------------------------------------------------------------------------
{
  register int_t k, noff, nband = 0;
  const int_t    next = nbndry / nel;

  for (k = 0, noff = 0; k < nel; k++) {
    nband = max (bandwidthSC (&bndmap[0]+noff, &bndmsk[0]+noff, next), nband);
    noff += next;
  }

  ++nband; // -- Diagonal.

  return nband;
}



int_t Nsys::bandwidthSC (const int_t* bmap,
			 const int_t* mask,
			 const int_t  next) const
// ---------------------------------------------------------------------------
// Find the global equation bandwidth of this element, excluding
// diagonal.
// ---------------------------------------------------------------------------
{
  register int_t i;
  register int_t Min  = INT_MAX;
  register int_t Max  = INT_MIN;

  for (i = 0; i < next; i++) {
    if (!mask[i]) {
      Min = min (bmap[i], Min);
      Max = max (bmap[i], Max);
    }
  }

  return Max - Min;
}


void Nsys::rebuild (FEML*       file  ,
		    const int_t optlev)
// ---------------------------------------------------------------------------
// The pressure numbering scheme gets rebuilt if the highest-numbered
// pressure node lies on the axis.  Before getting here, we are sure
// that this is the appropriate Nsys for the zero-mode pressure and
// has no essential BCs, so pressure system is singular.
// ---------------------------------------------------------------------------
{
  // -- Build a table of elements and sides that touch the axis.
  
  const char  axisBC = axial (file);
  const int_t nsurf (file->attribute ("SURFACES", "NUMBER"));
  char        tagc, tag[StrMax], err[StrMax];
  int_t       i, j, id, el, si, naxis = 0;

  for (i = 0; i < nsurf; i++) {
    while ((tagc = file->stream().peek()) == '#') // -- Skip comments.
      file->stream().ignore (StrMax, '\n');
    
    file->stream() >> id >> el >> si >> tag;
    if (strchr (tag, '<') && strchr (tag, '>') && strlen (tag) == 3) {
      if (tag[1] == 'B') {
	file->stream() >> tagc;
	if (tagc == axisBC) naxis++;
      }
    } else {
      sprintf (err, "unrecognized surface tag format: %s", tag);
      message (prog, err, ERROR);
    }

    file->stream().ignore (StrMax, '\n');
  }

  if (!naxis) return;

  axelmt.resize (naxis);
  axside.resize (naxis);

  file->attribute ("SURFACES", "NUMBER");

  for (j = 0, i = 0; i < nsurf; i++) {
    while ((tagc = file->stream().peek()) == '#') // -- Skip comments.
      file->stream().ignore (StrMax, '\n');
    
    file->stream() >> id >> el >> si >> tag;
    if (tag[1] == 'B') {
      file->stream() >> tagc;
      if (tagc == axisBC) {
	axelmt[j] = --el;
	axside[j] = --si;
	j++;
      }
    }

    file->stream().ignore (StrMax, '\n');
  }

  if (j != naxis) message (prog, "mismatch of axis surfaces", ERROR);

  // -- Now we are assured that there is at least one axial BC node.
  //    Before throwing out the numbering scheme, we must check if
  //    the highest-numbered node lies on the axis.

  if (highAxis()) {
    message(prog,"highest pressure node lies on axis. Renumbering...",REMARK);
    renumber (optlev, 1); 
  }
}


bool Nsys::highAxis() const
// ---------------------------------------------------------------------------
// Return true if the highest numbered pressure node lies on the axis,
// otherwise false.  Internal tables axelmt and axside must have been
// built in advance.
// ---------------------------------------------------------------------------
{
  const int_t pmax  = nglobal - 1;
  const int_t naxis = axelmt.size();
  int_t       i, j, loff, soff, elmtID, sideID;

  for (i = 0; i < nbndry; i++)
    if (bndmap[i] == pmax) {
      loff   = i % next_max;
      elmtID = (i - loff) / next_max;
      sideID = (loff - loff % np_max) / np_max;
      soff   = loff - sideID * np_max;
      for (j = 0; j < naxis; j++)
	if (axelmt[j] == elmtID && axside[j] == sideID) return true;
      if (soff == 0) {		// -- Check also end of CCW side.
	sideID = (sideID + 3) % 4;
	for (j = 0; j < naxis; j++)
	  if (axelmt[j] == elmtID && axside[j] == sideID) return true;
      }
    }
  
  return false;
}
