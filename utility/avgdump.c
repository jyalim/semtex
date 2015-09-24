/*****************************************************************************
 * avgdump.c:  Compute (running) averages of field files.
 *
 * Copyright (c) 1999 <--> $Date: 2015/04/20 11:14:18 $, Hugh Blackburn
 *
 * SYNOPSIS
 * --------
 * Input is two field files.  The first ("old.file") is assumed to
 * contain a running average of previous dumps, while the second
 * ("new.file") is assumed to contain new data that will be added into
 * the running average.  Both files must be binary (but not
 * necessarily matching) format and storage (number of elements,
 * interpolation orders, scalar fields) must conform.  The step number
 * in "old.file" is taken as the number of averages that have
 * contributed to it, and "new.file"'s data is added in with
 * appropriate weight.  Averaging is initialised using the
 * command-line flag "-i", which adds the data in the two files with
 * equal weight.  Only the first dump in each file is dealt with.
 *
 * A new binary file is written to stdout.  The step number is set to
 * reflect the number of averages that have been done to the data.
 *
 * With a non-zero relaxation factor given on the command line (-r
 * <eps>) do recursive relaxation using
 *
 * data_next = eps * data_new + (1 - eps) * data_old
 *
 * The amounts to first-order recursive filtering with a time constant
 * -ln(1 - eps).  In this case the step number in the old file is not
 * used although it gets updated and stored.
 *
 * USAGE
 * -----
 * avgdump [options] old.file new.file
 * options:
 * -h       ... print this message
 * -i       ... initialise averaging
 * -r <eps> ... weight new file by eps and old file by (1-eps)
 *
 * --
 * This file is part of Semtex.
 * 
 * Semtex is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 * 
 * Semtex is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Semtex (see the file COPYING); if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 *****************************************************************************/

static char RCS[] = "$Id: avgdump.c,v 8.1 2015/04/20 11:14:18 hmb Exp $";

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include <cfemdef.h>
#include <cveclib.h>

static char  prog[]    = "avgdump";
static const char* hdr_fmt[] = {	 /* -- Header output formatting. */
  "%-25s "             "Session\n",
  "%-25s "             "Created\n",
  "%-5d%-5d%-5d%-10d " "Nr, Ns, Nz, Elements\n",
  "%-25d "             "Step\n",
  "%-25.6g "           "Time\n",
  "%-25.6g "           "Time step\n",
  "%-25.6g "           "Kinvis\n",
  "%-25.6g "           "Beta\n",
  "%-25s "             "Fields written\n",
  "%-25s "             "Format\n"
};
typedef struct {		 /* -- Data information structure. */
  char     session [STR_MAX];
  char     creation[STR_MAX];
  int      nr               ;
  int      ns               ;
  int      nz               ;
  int      nel              ;
  int      step             ;
  double   time             ;
  double   timestep         ;
  double   kinvis           ;
  double   beta             ;
  char     field [STR_MAX]  ;
  char     format[STR_MAX]  ;
  double** data             ;
} Dump;


static void getargs   (int, char**, FILE**, FILE**, int*, double*);
static void getheader (FILE*, Dump*);
static void getdata   (FILE*, Dump*);
static void runavg    (Dump*, Dump*, int, double);
static void printup   (FILE*, Dump*);


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Driver routine.
 * ------------------------------------------------------------------------- */
{
  FILE   *oldfile = 0, *newfile = 0;
  Dump   *olddump = 0, *newdump = 0;
  int    init     = 0;
  double relax    = -1.0;	/* -- If positive we do relaxation. */

  getargs (argc, argv, &oldfile, &newfile, &init, &relax);

  olddump = (Dump*) calloc (1, sizeof (Dump));
  newdump = (Dump*) calloc (1, sizeof (Dump));

  getheader (oldfile, olddump);
  getheader (newfile, newdump);

  if (olddump -> nr  != newdump -> nr  ||
      olddump -> ns  != newdump -> ns  ||
      olddump -> nz  != newdump -> nz  ||
      olddump -> nel != newdump -> nel)
    message (prog, "structure of files don't match",           ERROR);

  if (!strstr (olddump -> field, newdump -> field))
    message (prog, "average fields don't match dumped fields", ERROR);

  getdata (oldfile, olddump);
  getdata (newfile, newdump);
  runavg  (olddump, newdump, init, relax);
  printup (stdout,  olddump);

  return (EXIT_SUCCESS);
}


static void getargs (int     argc   ,
		     char**  argv   ,
		     FILE**  oldfile,
		     FILE**  newfile,
		     int*    init   ,
		     double* relax  )
/* ------------------------------------------------------------------------- *
 * Parse command-line arguments.
 * ------------------------------------------------------------------------- */
{
  const char usage[] =
    "usage: avgdump [options] old.file new.file\n"
    "options:\n"
    "  -h       ... display this message\n"
    "  -i       ... initialise averaging\n"
    "  -r <eps> ... weight new file by eps and old file by (1-eps)\n";
    
  char err[STR_MAX], c;

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fprintf (stderr, usage); exit (EXIT_SUCCESS); break;
    case 'i':
      *init = 1; break;
    case 'r': 
      if (*++argv[0]) *relax = atof (*argv);
      else {*relax = atof (*++argv); argc--;}
      break;
    default:
      sprintf (err, "illegal option: %c\n", c);
      message (prog, err, ERROR);
      break;
    }

  if (argc == 2) {
    *oldfile = efopen (argv[0], "r");
    *newfile = efopen (argv[1], "r");
  } else {
    fprintf (stderr, usage); exit (EXIT_FAILURE);
  }  
}


static void getheader (FILE* f,
		       Dump* h)
/* ------------------------------------------------------------------------- *
 * Find header information.
 * ------------------------------------------------------------------------- */
{
  char buf[STR_MAX];

  fgets  (h -> session,  STR_MAX, f);
  fgets  (h -> creation, STR_MAX, f);
  fscanf (f, "%d %d %d %d", &h -> nr, &h -> ns, &h -> nz, &h -> nel);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%d", &h -> step);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%lf", &h -> time);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%lf", &h -> timestep);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%lf", &h -> kinvis);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%lf", &h -> beta);
  fgets  (buf, STR_MAX, f);
  fscanf (f, "%s", h -> field);
  fgets  (buf, STR_MAX, f);
  fgets  (h -> format, STR_MAX, f);

  if (!strstr (h -> format,      "binary"))
    message (prog, "input field file not in binary format",     ERROR);
  else if (!strstr (h -> format, "endian"))
    message (prog, "input field file in unknown binary format", WARNING);
}


static void getdata (FILE* f,
		     Dump* h)
/* ------------------------------------------------------------------------- *
 * Find the number of fields, allocate storage & do binary read.
 * ------------------------------------------------------------------------- */
{
  char      localfmt[STR_MAX];
  int       i, swab;
  const int npts    = h -> nr * h -> ns * h -> nz * h -> nel;
  const int nfields = strlen (h -> field);
  const int ntot    = npts * nfields;

  h -> data = dmatrix (0, nfields - 1, 0, npts - 1);

  if (fread (h -> data[0], sizeof (double), ntot, f) != ntot) {
    sprintf (localfmt, "could not read fields: %s", h -> field);
    message (prog, localfmt, ERROR);
  }

  format (localfmt);
  swab = (strstr (h -> format, "big") && strstr (localfmt,    "little")) || 
         (strstr (localfmt,    "big") && strstr (h -> format, "little"));

  if (swab) dbrev (ntot, h -> data[0], 1, h -> data[0], 1);
}


static void runavg (Dump*  a   ,
		    Dump*  b   ,
		    int    init,
		    double eps )
/* ------------------------------------------------------------------------- *
 * Running average of a & b, leave in a.  Optionally initialise, and if
 * eps is positive, do recursive filtering as outlined in file header.
 * ------------------------------------------------------------------------- */
{
  int       i;
  double    fac;
  const int nfields = strlen (a -> field);
  const int npts    = a -> nr * a -> ns * a -> nz * a -> nel;

  if (eps > EPSSP) {
    if (init) a -> step = 2; else a -> step++;
    for (i = 0; i < nfields; i++) {
      dsmul  (npts, 1.0 - eps, a -> data[i], 1, a -> data[i], 1);
      dsvtvp (npts, eps, b -> data[i], 1, a -> data[i], 1, a -> data[i], 1);
    }
  } else {
    if (init) {
      fac = 0.5;
      a -> step = 2;
      for (i = 0; i < nfields; i++)
	dsvvpt (npts, fac, a -> data[i], 1, b -> data[i], 1, a -> data[i], 1);
    } else {
      fac = (double) a -> step;
      a -> step++;
      for (i = 0; i < nfields; i++) {
	dsvtvp (npts, fac, a -> data[i], 1, b -> data[i], 1, a -> data[i], 1);
	dsmul  (npts, 1.0 / (fac + 1.0),    a -> data[i], 1, a -> data[i], 1);
      }
    }
  }
}


static void printup (FILE* f,
		     Dump* h)
/* ------------------------------------------------------------------------- *
 * Output the modified data.
 * ------------------------------------------------------------------------- */
{
  int       i;
  const int ntot = h -> nr * h -> ns * h -> nz * h -> nel;
  const int nfields = strlen (h -> field);
  char      s1[STR_MAX], s2[STR_MAX];
  time_t    tp = time (0);

  fprintf  (f, h -> session);
  strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
  sprintf  (s1, hdr_fmt[1], s2);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[2], h -> nr, h -> ns, h -> nz, h -> nel);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[3], h -> step);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[4], h -> time);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[5], h -> timestep);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[6], h -> kinvis);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[7], h -> beta);
  fprintf  (f, s1);
  sprintf  (s1, hdr_fmt[8], h -> field);
  fprintf  (f, s1);
  sprintf  (s2, "binary "); format (s2 + strlen (s2));
  sprintf  (s1, hdr_fmt[9], s2);
  fprintf  (f, s1);

  for (i = 0; i < nfields; i++) 
    if (fwrite (h -> data[i], sizeof (double), ntot, f) != ntot)
      message (prog, "error occurred while writing", ERROR);
}
