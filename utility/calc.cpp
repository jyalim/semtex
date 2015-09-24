///////////////////////////////////////////////////////////////////////////////
// calc: a basic calculator using the function parser.
//
// Copyright (c) 1994 <--> $Date: 2015/04/20 11:14:18 $, Hugh Blackburn
//
// Usage: calc [-h] [file]
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

static char RCS[] = "$Id: calc.cpp,v 8.1 2015/04/20 11:14:18 hmb Exp $"; 

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#include <cfemdef.h>
#include <utility.h>
#include <femlib.h>

static char prog[] = "calc";
static void getargs (int, char**, istream*&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char     buf[StrMax];
  istream* input;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, input);

  while (input -> getline (buf, FILENAME_MAX))
    if (strstr (buf, "="))
      Femlib::value (buf);
    else
      cout << setprecision(17) << Femlib::value (buf) << endl;
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int       argc ,
		     char**    argv ,
		     istream*& input)
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char buf[StrMax];
  char usage[] = "Usage: %s [-h] [file]\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cerr << "-- Calculator operators, functions and procedures:" << endl;
      yy_help ();
      cerr << endl << "-- Preset internal variables:"  << endl;
      yy_show ();
      exit (EXIT_SUCCESS);
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> bad()) {
      cerr << usage;
      sprintf (buf, "unable to open file: %s", *argv);
      message (prog, buf, ERROR);
    }
  } else input = &cin;
}

