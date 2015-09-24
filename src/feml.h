#ifndef FEML_H
#define FEML_H
///////////////////////////////////////////////////////////////////////////////
// feml: Finite Element markup Language (FEML) header file.
//
// Copyright (c) 1997 <--> $Date: 2015/09/08 21:43:00 $, Hugh Blackburn
//
// $Id: feml.h,v 8.2 2015/09/08 21:43:00 hmb Exp $
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

using namespace std;

#include <cfemdef.h>

#define FEML_KEYWORD_MAX 32


class FEML
// ===========================================================================
// Routines provide facilities to position input stream at location of
// a given keyword.
//
// This is the list of currently-defined (reserved) section keywords:
//   TOKENS
//   FIELDS
//   GROUPS
//   BCS
//   NODES
//   ELEMENTS
//   SURFACES
//   CURVES
//   USER
//   HISTORY
//   CUT
//   FORCE
// Keywords are stored upper case, input is case-insensitive.
// The FEML class does not require that any of the above sections are actually
// used in an input file, it just treats them as reserved section tag-names.
// ===========================================================================
{

public:
  FEML  (const char*);
  ~FEML () { _feml_file.close(); }

  bool valueFromSection (real_t *, const char *, const char *);
  bool valueFromSection (int_t  *, const char *, const char *);
  bool valueFromSection (char   *, const char *, const char *);

  bool        seek      (const char*);
  int_t       attribute (const char*, const char*);
  istream&    stream    ()       { return _feml_file; }
  const char* root      () const { return _feml_root; }
  int_t       sections  (vector <const char*>&);
  bool        echo      (ostream&, const char*);

private:
  char*     _feml_root;		        // Name of FEML file, suffix removed.
  ifstream  _feml_file;		        // Input stream.
  streampos _keyPosn[FEML_KEYWORD_MAX]; // Locations corresponding to keywords.
  char*     _keyWord[FEML_KEYWORD_MAX]; // Keywords used.
  int_t     _nKey;                      // Number of keywords used.

  bool      tokens ();
  void      check_ASCII();
};

#endif

