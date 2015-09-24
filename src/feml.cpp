///////////////////////////////////////////////////////////////////////////////
// feml.C:  Finite Element Markup Language (FEML) routines.
//
// Copyright (c) 1994 <--> $Date: 2015/09/08 21:43:00 $, Hugh Blackburn
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
// --
//
// After initialization, FEML files are prescanned to find locations of
// keywords.  These locations are stored, and reset by the seek function.
//
// Keywords are set in SGML section-block notation, with open and close tags.
//
// <keyword [options]>
//   .
//   .
//   .
// </keyword>
//
// 1. Keywords case-insensitive on input, converted to uppercase internally.
// 2. There are reserved keywords; see feml (and below).
// 3. There is a maximum number of keywords (KEYWORDS_MAX) set in feml.
// 4. No restriction is placed on options, or on the contents of a section.
// 5. After seeking, input stream is repositioned just following keyword.
// 6. Because of the repositioning requirement, feml cannot work with a
//    stream like standard input: input must come from a file.
//
// TO DO: move to XML.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: feml.cpp,v 8.2 2015/09/08 21:43:00 hmb Exp $";

#include <cstdlib>		// -- C standard headers.
#include <cstdio>
#include <cstring>
#include <cctype>

#include <vector>		// -- C++ STL.

using namespace std;

#include <cfemdef.h>		// -- Semtex headers.
#include <utility.h>
#include <feml.h>
#include <femlib.h>


FEML::FEML (const char* session) :
  _nKey (0)
// ---------------------------------------------------------------------------
// Attempt to open session file, prescan to locate sections.  Set
// _feml_root.  Load tokens.
//
// The array _keyWord is initially filled with all the reserved words,
// but these are not necessarily present in the session file. If a
// section is actually present in the file, that is flagged by a
// non-NULL corresponding stream pointer in array _keyPosn.
// ---------------------------------------------------------------------------
{
  const char routine[] = "FEML::FEML";
  char       c, err[STR_MAX], key[STR_MAX], yek[STR_MAX];
  char*      u;
  int_t      i, N;
  bool       OK, found;

  const char* reserved[] = {
    "TOKENS",
    "FIELDS",
    "GROUPS",
    "BCS",
    "NODES",
    "ELEMENTS",
    "SURFACES",
    "CURVES",
    "USER",
    "HISTORY",
    "CUT",
    "FORCE",
    0
  };

  _feml_file.open (session);

  if (!_feml_file) {
    sprintf (err, "couldn't open session file %s", session);
    message (routine, err, ERROR);
  }

  check_ASCII();
  
  strcpy ((_feml_root = new char [strlen (session) + 1]), session);

  for (i = 0; reserved[i] && i < FEML_KEYWORD_MAX; i++) {
    _keyWord[i] = strcpy ((new char [strlen (reserved[i]) + 1]), reserved[i]);
    _keyPosn[i] = NULL;
  }
  _keyWord[i] = NULL;

  if (i == FEML_KEYWORD_MAX) {
    sprintf (err, "Number of reserved keywords exceeds table size (%1d)", i);
    message (routine, err, ERROR);
  }
  
  while (_feml_file >> c) {

    if (c == '<') {

      // -- Next word may be a keyword.
      //    First massage it to get correct form.

      _feml_file >> key;

      N = strlen (key);
      if (key[N - 1] == '>') {
	c = '>';
	_feml_file.putback (c);
	key[N - 1] = '\0';
	N--;
      }

      u = key; while (*u = toupper (*u)) u++;

      // -- Check if key is a keyword and if so:
      //    1. install file position immediately following keyword in table;
      //    2. move on to find closing tag.

      for (found = false, i = 0; !found && _keyWord[i]; i++) {

	if (strcmp (key, _keyWord[i]) == 0) {

	  // -- Locate closing '>'.

	  _keyPosn[i] = _feml_file.tellg ();
	  while (!found && _feml_file >> c) found = c == '>';

	  if (!found) {
	    sprintf (err, "closing '>' not found for keyword %s", key);
	    message (routine, err, ERROR);
	  }

	  // -- Locate closing "</key>".

	  OK = false;

	  while ((!OK) && (_feml_file >> c)) {
	    if (c == '<') {
	      _feml_file >> c;

	      if (c == '/') {
		_feml_file >> yek;

		N = strlen (yek);
		if (yek[N - 1] == '>') {
		  c = '>';
		  _feml_file.putback (c);
		  yek[N - 1] = '\0';
		  N--;
		}

		u = yek; while (*u = toupper (*u)) u++;

		if (OK = strcmp (key, yek) == 0) {
		  while (_feml_file >> c) if (c == '>') break;

		  if (c != '>') {
		    sprintf (err, "closing '>' not found for /%s", key);
		    message (routine, err, ERROR);
		  }
		}
	      }
	    }
	  }

	  if (!OK) {
	    sprintf (err, "couldn't locate </%s> to match <%s>", key, key);
	    message (routine, err, ERROR);
	  } else ++_nKey;
	}
      }
    }
  }

  if (!found) message (routine, "no keywords located", ERROR);

  _feml_file.clear ();		// -- Reset EOF error condition.
  _feml_file.seekg (0);		// -- And rewind.

  tokens ();			// -- Initialize Femlib parser.
}


void FEML::check_ASCII()
// ---------------------------------------------------------------------------
// Check for non-ASCII characters in session file.
// ---------------------------------------------------------------------------
{
  const char routine[] = "FEML::check_ASCII";
  int col, line = 0;
  char the_line[STR_MAX], err[STR_MAX];
  
  while (_feml_file.getline (the_line, STR_MAX)) {
    line ++;
    if (the_line[0] == '#') continue;
    for (col = 0; col < strlen(the_line); col++)
      if (!isascii(the_line[col])) {
        sprintf (err, "Non-ASCII character on line %i, col %i", line, col+1);
        message(routine, err, ERROR);
      }
  }
  _feml_file.clear ();          // -- Reset EOF error condition.
  _feml_file.seekg (0);         // -- And rewind.
}


bool FEML::seek (const char* keyword)
// ---------------------------------------------------------------------------
// Look for keyword in stored table.
// If present, stream is positioned after keyword and true is returned.
// If not, stream is rewound and false is returned.
// ---------------------------------------------------------------------------
{
  int_t i;
  bool  found = false;

  for (i = 0; !found && _keyWord[i]; i++)
    found = (strstr   (keyword, _keyWord[i]) != 0 &&
	     static_cast<int_t>(_keyPosn[i]) != 0);

  if (!found) {
    _feml_file.clear ();
    _feml_file.seekg (0);
    return false;
  } else {
    _feml_file.clear ();
    _feml_file.seekg (_keyPosn[i - 1]);
  }

  return true;
}


int_t FEML::attribute (const char* tag ,
		       const char* attr)
// ---------------------------------------------------------------------------
// Tag attributes are given as options in form <tag attr=int [attr=int ...]>
// Return integer value following '='.  No whitespace allowed in attributes.
// On return, FEML's stream is set to start of next line.
// ---------------------------------------------------------------------------
{
  const char routine[] = "FEML::attribute";
  char       buf[STR_MAX], err[STR_MAX];
  char*      v;
  int_t      n = 0;

  if (!seek (tag)) {
    sprintf (err, "couldn't locate tag %s in feml file", tag);
    message (routine, err, ERROR);
  }

  while (_feml_file >> buf)
    if (strstr (buf, attr)) {
      if (buf[strlen (buf) - 1] == '>') buf[strlen (buf) - 1] = '\0';
      v = buf;
      while (*v && *v++ != '=');
      if (!(*v)) {
	sprintf (err, "attribute syntax error in %s", buf);
	message (routine, err, ERROR);
      }

      n = atoi (v);
      break;
    } else if (strchr (buf, '>')) {
      sprintf (err, "%s not found in tag %s", attr, tag);
      message (routine, err, ERROR);
    }

  _feml_file.ignore (STR_MAX, '\n');

  return n;
}


bool FEML::tokens ()
// ---------------------------------------------------------------------------
// Install token table.  Return false if no TOKEN section is found.
// NUMBER attribute ignored if present.  Fix any inconsistent values.
// Parser must have been initialized before entry.
// Lines with '#' at the start are ignored.
// ---------------------------------------------------------------------------
{
  const char     routine[] = "FEML::tokens";
  char           buf[STR_MAX];
  register char* u;

  if (seek ("TOKENS")) {
    _feml_file.ignore (STR_MAX, '\n');

    while (_feml_file.getline (buf, STR_MAX)) {
      if (strstr (buf, "=") && buf[0] != '#') Femlib::value (buf);
      u = buf; while (*u = toupper (*u)) u++;
      if (strstr (buf, "TOKENS")) break;
    }

    if (Femlib::ivalue ("IO_FLD") > Femlib::ivalue ("N_STEP"))
        Femlib::ivalue ("IO_FLD",   Femlib::ivalue ("N_STEP"));

    if (Femlib::ivalue ("N_TIME") > 3) {
      message (routine, "N_TIME too large, reset to 3", WARNING);
      Femlib::ivalue ("N_TIME", 3);
    }

    return true;
  }

  return false;
}


int_t FEML::sections (vector <const char*>& present)
// ---------------------------------------------------------------------------
// Return a vector with pointers to the names of the sections that
// were present in the session file.
// ---------------------------------------------------------------------------
{
  int_t i = 0, j = 0;

  present.resize (_nKey);

  while (_keyWord[i]) { if (_keyPosn[i]) present[j++] = _keyWord[i]; i++; }

  return _nKey;
}

bool FEML::echo (ostream&    stream,
		 const char* key   )
// ---------------------------------------------------------------------------
// Print out the text of appropriate FEML file section on stream.
// Do nothing if named section is not present.
// ---------------------------------------------------------------------------
{
  const char routine[] = "FEML::echo";
  int_t loc, N;
  char  c, yek[128], err[128], *u;
  bool  found = false;

  for (loc = 0; !found && _keyWord[loc]; loc++)
    found = (strstr       (key, _keyWord[loc]) != 0 &&
	     static_cast<int_t>(_keyPosn[loc]) != 0);

  if (!found) { _feml_file.clear (); _feml_file.seekg (0); return false; }

  --loc;
  _feml_file.clear ();
  _feml_file.seekg (_keyPosn[loc]);

  stream << '<' << _keyWord[loc];

  _feml_file >> noskipws;

  _feml_file >> c; stream << c;

  // -- We already know the structure of the file is valid and the
  //    present section will be appropriately closed.  We now scan
  //    forward and print up until that point is reached.

  found = false;

  while ((!found) && (_feml_file >> c)) {
    stream << c;
    if (c == '<') {
      _feml_file >> c; stream << c;

      if (c == '/') {
	_feml_file >> yek;

	N = strlen (yek);
	if (yek[N - 1] == '>') {
	  c = '>';
	  _feml_file.putback (c);
	  yek[N - 1] = '\0';
	  N--;
	}

	u = yek; while (*u = toupper (*u)) u++;

	stream << yek;

	if (found = strcmp (key, yek) == 0) {
	  while (_feml_file >> c) { stream << c; if (c == '>') break; }

	  if (c != '>') {
	    sprintf (err, "closing '>' not found for /%s", key);
	    message (routine, err, ERROR);
	  }
	}
      }
    }
  }
  stream << endl;

  if (!found) {
    sprintf (err, "couldn't locate </%s> to match <%s>", key, key);
    message (routine, err, ERROR);
  }

  _feml_file >> skipws;

        return true;
}

// added by Thomas Albrecht for evaluate force section


bool FEML::valueFromSection (real_t     *value  ,
			     const char *section,
                             const char *token  )
// ---------------------------------------------------------------------------
// Search for 'token' in 'section' of input file. If found, copy to 'value',
// otherwise, 'value' is unchanged.
// Returns TRUE on sucess, FALSE if section or token not found.
// ---------------------------------------------------------------------------
{
  char routine[]   = "FEML::valueFromSection";
  char endsection[StrMax], s[StrMax];

  sprintf(endsection, "</%s>", section);

  if (seek (section)) {
    stream().ignore (StrMax, '\n');

    while (!stream().eof())
    {
      stream().getline(s, StrMax);
      if (s[0] == '#')
      {
        stream().ignore (StrMax, '\n');
        continue;
      }
      if (strcmp (s, endsection) == 0) break;
      char *tok;
      if ((tok = strtok (s, "\t= ")) == NULL) continue;

      if (strcmp (tok, token) == 0) {
        tok = strtok (NULL, "=\n");
        *value = Femlib::value(tok);
        return true;
      }
    }
    return false; // token not found
  }  else {
    char msg[StrMax];
    sprintf(msg, "%s section not found", section);
    message (routine, msg,  ERROR);
  }

  return false;
}


bool FEML::valueFromSection (int_t      *value  ,
			     const char *section,
			     const char *token  )
// ---------------------------------------------------------------------------
// as above, integer version
// ---------------------------------------------------------------------------
{
  char routine[]   = "FEML::valueFromSection";
  char endsection[StrMax], s[StrMax];

  sprintf(endsection, "</%s>", section);

  if (seek (section)) {
    stream().ignore (StrMax, '\n');

    while (!stream().eof())
    {
      stream().getline(s, StrMax);
      if (s[0] == '#')
      {
        stream().ignore (StrMax, '\n');
        continue;
      }
      if (strcmp (s, endsection) == 0) break;
      char *tok;
      if ((tok = strtok (s, "\t= ")) == NULL) continue;

      if (strcmp (tok, token) == 0) {
        tok = strtok (NULL, "=\n");
        *value = Femlib::ivalue(tok);
        return true;
      }
    }
    return false; // token not found
  }  else {
    char msg[StrMax];
    sprintf(msg, "%s section not found", section);
    message (routine, msg,  ERROR);
  }

  return false;
}


bool FEML::valueFromSection (char       *buf    ,
			     const char *section,
			     const char *token  )
// ---------------------------------------------------------------------------
// as above, string version
// ---------------------------------------------------------------------------
{
  char routine[]   = "FEML::valueFromSection";
  char endsection[StrMax], s[StrMax];

  sprintf(endsection, "</%s>", section);

  if (seek (section)) {
    stream().ignore (StrMax, '\n');

    while (!stream().eof())
    {
      stream().getline(s, StrMax);
      if (s[0] == '#')
      {
        stream().ignore (StrMax, '\n');
        continue;
      }
      if (strcmp (s, endsection) == 0) break;
      char *tok;
      if ((tok = strtok (s, "\t= ")) == NULL) continue;

      if (strcmp (tok, token) == 0) {
        tok = strtok (NULL, "=\n");
        strcpy(buf, tok);
        return true;
      }
    }
    return false; // token not found
  }  else {
    char msg[StrMax];
    sprintf(msg, "%s section not found", section);
    message (routine, msg,  ERROR);
  }

  return false;
}
