///////////////////////////////////////////////////////////////////////////////
// rectmesh.C: create a session file for a rectangular mesh of elements.
//
// Copyright (c) 2000 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn
//
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
//
// Usage:
// -----
// rectmesh [-b <num>] [-e <num>] [-v <num>] [file]
//   -b <num> ... output in <num> blocks, contiguous in x. [Default: 1]
//   -e <num> ... offset first element number by <num>.
//   -v <num> ... offset first vertex number by <num>.
//
// Files:
// -----
// Input consists of a list of x, followed by y, locations of element
// boundaries, one per line.  A single blank line separates x from y
// locations.  Output consists of a (2D) session file with an element
// order of 7, and "wall" group boundaries around the domain border.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: rectmesh.cpp,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <sem.h>

static char prog[] = "rectmesh";
static void getargs (int, char**, int_t&, int_t&, int_t&, istream*&);
static void header  ();


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char                    line[STR_MAX];
  istream*                input;
  real_t                  x, y;
  stack <real_t>          X, Y;
  vector<vector<Point*> > vertex;
  int_t                   Nx = 0, Ny = 0, Nb = 1, NelB;
  int_t                   eOffset = 0, vOffset = 0;
  int_t                   i, j, k, b;
  string                  s;

  getargs (argc, argv, Nb, eOffset, vOffset, input);

  // -- Read x, then y locations onto two stacks.

  while (input -> getline(line, STR_MAX).gcount() > 1) {
    istringstream ss (s = line);
    ss >> x;
    X.push (x);
    Nx++;
  }

  if ((Nx - 1) % Nb)
    message (prog,
	     "Nx-1 (No. of elements in x) must be an integer multiple of Nb",
	     ERROR);

  while (input -> getline(line, STR_MAX)) {
    istringstream ss (s = line);
    ss >> y;
    Y.push (y);
    Ny++;
  }
  
  NelB = (Nx-1)/Nb*(Ny-1);

  // -- Insert into vertex matrix.

  vertex.resize (Ny);
  for (i = 0; i < Ny; i++) vertex[i].resize(Nx);
  for (i = 0; i < Ny; i++)
    for (j = 0; j < Nx; j++) {
      vertex[i][j] = new Point;
      vertex[i][j] -> z = 0.0;
    }

  j = Nx;
  while (j--) { vertex[0][j] -> x = X.top(); X.pop(); }
  i = Ny;
  while (i--) { vertex[i][0] -> y = Y.top(); Y.pop(); }

  for (i = 0; i < Ny; i++)
    for (j = 0; j < Nx; j++) {
      vertex[i][j] -> x = vertex[0][j] -> x;
      vertex[i][j] -> y = vertex[i][0] -> y;
    }

  header();

  // -- Print up vertex list.

  cout << "<NODES NUMBER=" << Nx*Ny << ">" << endl;

  for (k = vOffset, i = 0; i < Ny; i++)
    for (j = 0; j < Nx; j++)
      cout << setw(5)  << ++k << "\t"
	   << setw(15) << vertex[i][j] -> x
	   << setw(15) << vertex[i][j] -> y
	   << setw(15) << vertex[i][j] -> z
	   << endl;
  
  cout << "</NODES>" << endl;

  // -- Print up elements.

  cout << endl << "<ELEMENTS NUMBER=" << NelB*Nb << ">" << endl;

#if 1
  for (k = eOffset+1, b = 0; b < Nb; b++)
    for (i = 0; i < (Ny - 1); i++)
      for (j = b*((Nx-1)/Nb); j < (b+1)*((Nx-1)/Nb); j++, k++)
	cout << setw(5) << k << "\t" << "<Q>"
	     << setw(5) << vOffset + j +  i      * Nx + 1
	     << setw(5) << vOffset + j +  i      * Nx + 2
	     << setw(5) << vOffset + j + (i + 1) * Nx + 2
	     << setw(5) << vOffset + j + (i + 1) * Nx + 1
	     << "    </Q>" << endl;
#else
  for (k = eOffset+1, i = 0; i < (Ny - 1); i++)
    for (j = 0; j < (Nx - 1); j++, k++)
      cout << setw(5) << k << "\t" << "<Q>"
	   << setw(5) << vOffset + j +  i      * Nx + 1
	   << setw(5) << vOffset + j +  i      * Nx + 2
	   << setw(5) << vOffset + j + (i + 1) * Nx + 2
	   << setw(5) << vOffset + j + (i + 1) * Nx + 1
	   << "    </Q>" << endl;
#endif
    
  cout << "</ELEMENTS>" << endl;

  // -- Print up surfaces.

  cout << endl << "<SURFACES NUMBER=" << 2*((Nx-1)+(Ny-1)) << ">" << endl;

#if 1
  for (k = 1, b = 0; b < Nb; b++)
    for (j = 0; j < (Nx-1)/Nb; j++, k++)
      cout << setw(5) << k << setw(5) 
	   << eOffset + b*NelB + j + 1
	   << "    1"
	   << "    <B> w </B>" << endl;
  for (i = 0; i < (Ny - 1); i++, k++)
    cout << setw(5) << k 
	 << setw(5) << eOffset + (Nb - 1)*NelB + (i + 1)*(Nx - 1)/Nb
	 << "    2"
	 << "    <B> w </B>" << endl;
  for (b = Nb; b > 0; b--)
    for (j = (Nx-1)/Nb; j > 0; j--, k++)
      cout << setw(5) << k 
	   << setw(5) <<  eOffset + b*NelB - (Nx - 1)/Nb + j
	   << "    3"
	   << "    <B> w </B>" << endl;
  for (i = Ny - 1; i > 0; i--, k++)
    cout << setw(5) << k 
	 << setw(5) << eOffset +  (i - 1) * (Nx - 1)/Nb + 1   
	 << "    4"
	 << "    <B> w </B>" << endl;
#else
  for (k = 1, j = 0; j < (Nx - 1); j++, k++)
    cout << setw(5) << k << setw(5) 
	 <<  eOffset + j + 1
	 << "    1"
	 << "    <B> w </B>" << endl;
  for (i = 0; i < (Ny - 1); i++, k++)
    cout << setw(5) << k 
	 << setw(5) << eOffset + (i + 1) * (Nx - 1)
	 << "    2"
	 << "    <B> w </B>" << endl;
  for (j = Nx - 1; j > 0; j--, k++)
    cout << setw(5) << k 
	 << setw(5) <<  eOffset + j + (Nx - 1) * (Ny - 2)
	 << "    3"
	 << "    <B> w </B>" << endl;
  for (i = Ny - 1; i > 0; i--, k++)
    cout << setw(5) << k 
	 << setw(5) << eOffset + (i - 1) * (Nx - 1) + 1   
	 << "    4"
	 << "    <B> w </B>" << endl;
#endif
    
  cout << "</SURFACES>" << endl;

  return EXIT_SUCCESS;
}


static void getargs (int       argc ,
		     char**    argv ,
		     int_t&    Nb   ,
		     int_t&    eOff ,
                     int_t&    vOff ,
		     istream*& input)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: rectmesh [options] [file]\n"
    "  options:\n"
    "  -h       ... print this message\n"
    "  -b <num> ... output in <num> blocks, contiguous in x [Default: 1]\n"
    "  -e <num> ... offset first element number by <num>\n"
    "  -v <num> ... offset first vertex number by <num>\n";

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'b':
      if (*++argv[0]) Nb = atoi (*argv);
      else {Nb = atoi (*++argv); argc--;}
      break;
    case 'e':
      if (*++argv[0]) eOff = atoi (*argv);
      else {eOff = atoi (*++argv); argc--;}
      break;
    case 'v':
      if (*++argv[0]) vOff = atoi (*argv);
      else {vOff = atoi (*++argv); argc--;}
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> bad()) message (prog, "unable to open geometry file", ERROR);
  } else input = &cin;
}


static void header ()
// ---------------------------------------------------------------------------
// Output header information.  Mesh is valid (eg for meshpr) without this.
// ---------------------------------------------------------------------------
{
  cout << "<FIELDS>\n\tu\tv\tp\n</FIELDS>" << endl << endl;

  cout << "<GROUPS NUMBER=1>\n\t1\tw\twall\n</GROUPS>" << endl << endl;

  cout << "<BCS NUMBER=1>\n\t1\tw\t3" << endl;
  cout << "\t\t\t<D>\tu = 0.0\t</D>"  << endl;
  cout << "\t\t\t<D>\tv = 0.0\t</D>"  << endl;
  cout << "\t\t\t<H>\tp = 0.0\t</H>"  << endl;
  cout << "</BCS>"                    << endl << endl;

  cout << "<TOKENS>"                << endl;
  cout << "\tKINVIS    = 2e-6"      << endl;
  cout << "\tD_T       = 0.005"     << endl;
  cout << "\tN_STEP    = 100"       << endl;
  cout << "\tN_TIME    = 2"         << endl;
  cout << "\tN_P       = 7"         << endl;
  cout << "\tN_Z       = 1"         << endl;
  cout << "\tBETA      = 1.0"       << endl;
  cout << "\tIO_CFL    = 50"        << endl;
  cout << "\tIO_FLD    = 1000"      << endl;
  cout << "\tCHKPOINT  = 1"         << endl;
  cout << "</TOKENS>"               << endl << endl;
}
