/*****************************************************************************
 * convert:  Format conversion program for sem-compatible field files.
 *
 * Copyright (c) 1994 <--> $Date: 2015/04/20 11:14:18 $, Hugh Blackburn
 *
 * USAGE: convert [-h] [-format] [-v] [-n dump] [-o output] [-z] [input[.fld]
 *
 * Default action is to convert binary to ASCII files or vice versa.
 *
 * Optional argument format can be one of:
 * -a: force ASCII output;
 * -b: force binary output in current machine IEEE format;
 * -s: force binary output in byte-swapped    IEEE format.
 *
 * If both -s & -b are specified then whichever is last on the comand line
 * takes precedence. -z sets the Step and Time to zero.
 *
 * Each input is read into an internal buffer in machine's double binary
 * format prior to output.
 *
 * Output is either to stdout or an optional file argument.
 *
 * sample                      Session
 * Mon Apr 22 18:23:13 1991    Created
 * 9 9 1 8                     Nr, Ns, Nz, Nel
 * 50                          Step
 * 0.05                        Time
 * 0.001                       Time step
 * 0.025                       Kinvis
 * 1                           Beta-z
 * uvp                         Fields written
 * ASCII                       Format
 *
 * Other recognized formats are: ascii (for backwards compatibility)
 *                               binary/BINARY assumed to be machine's default.
 *                               binary IEEE little-endian
 *                               binary IEEE big-endian
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

static char RCS[] = "$Id: convert.c,v 8.1 2015/04/20 11:14:18 hmb Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef enum { UNKNOWN, ASCII, IEEE_BIG, IEEE_LITTLE } FORMAT;

static char  prog[]  = "convert";

static int   verbose = 0;
static char  usage[] = "Usage: convert [-format] [-h] [-v] [-o output] "
                       "[input[.fld]]\n"
                       "format can be one of:\n"
                       "  -a ... force ASCII output\n"
                       "  -b ... force IEEE-binary output\n"
                       "  -s ... force IEEE-binary output (byte-swapped)\n"
                       "other options are:\n"
                       "  -h        ... print this message\n"
                       "  -n dump   ... select dump number\n"
                       "  -v        ... be verbose\n"
                       "  -o output ... output to named file\n"
                       "  -z        ... zero Time and Step in output\n";
    

static void   getargs      (int, char**, FILE**, FILE**, FORMAT*, int*, int*);
static void   error        (const char*);
static void   get_data     (FILE*, const int, const FORMAT, const FORMAT,
			    const int, const int, double**);
static void   put_data     (FILE*, const FORMAT, const FORMAT,
			    const int, const int, double**);
static void   dswap        (const int, double*);
static int    count_fields (const char*);
static FORMAT architecture (void);
static FORMAT classify     (const char*);
static int    iformat      (void);


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Driver routine.
 * ------------------------------------------------------------------------- */
{
  char     buf[BUFSIZ];
  double** data;
  int      nfields, npts, n, nr, ns, nz, nel;
  int      ndump = 0, nread = 0, selected = 1, zero = 0;
  FILE*    fp_in    = stdin;
  FILE*    fp_out   = stdout;
  FORMAT   inputF   = UNKNOWN,
           outputF  = UNKNOWN,
           machineF = architecture();

  getargs (argc, argv, &fp_in, &fp_out, &outputF, &ndump, &zero);

  while (fgets (buf, BUFSIZ, fp_in)) {

    if (ndump) selected = ndump == ++nread;

    /* -- Find size information. */

    n = 3;
    while (--n) {
      if (selected) fputs (buf, fp_out);
      fgets (buf, BUFSIZ, fp_in);
    }

    if (sscanf(buf, "%d%d%d%d", &nr, &ns, &nz, &nel) != 4)
      error ("unable to read the file size");               
    npts = nr * ns * nz * nel;
    
    /* -- Zero time and step if required. */

    if (selected) fputs (buf, fp_out);
    fgets (buf, BUFSIZ, fp_in);
    if (zero) sprintf (buf, "%-25d Step\n", 0);
    if (selected) fputs (buf, fp_out);
    fgets (buf, BUFSIZ, fp_in);
    if (zero) sprintf (buf, "%-25.6g Time\n", 0.0);
    
    n = 5;
    while (--n) {
      if (selected) fputs (buf, fp_out);
      fgets (buf, BUFSIZ, fp_in);
    }

    nfields = count_fields (buf);
    if (selected) fputs (buf, fp_out);

    /* -- Set input and ouput formats. */

    fgets (buf, BUFSIZ, fp_in);
    inputF = classify (buf);

    if (outputF == UNKNOWN) outputF = (inputF == ASCII) ? machineF : ASCII;

    /* -- Allocate storage. */

    data = (double**) malloc (nfields * sizeof (double*));
    for (n = 0; n < nfields; n++) 
      data[n] = (double*) malloc (npts*sizeof (double));

    /* -- Do input & output of field data. */

    if (verbose)
      fprintf (stderr, "%s: converting %1d fields, %1d points ",
	       prog, nfields, npts);
    
    get_data (fp_in, selected, inputF, machineF, npts, nfields, data);
    if (selected) put_data (fp_out, inputF, outputF,  npts, nfields, data);

    /* -- Deallocate storage. */

    for (n = 0; n < nfields; n++)
      free (data[n]);
    free (data);

    if (ndump && selected) /* -- No need to do any more work: quit. */ break;
  }

  return EXIT_SUCCESS;
}


static int count_fields (const char *s)
/* ------------------------------------------------------------------------- *
 * Count the number of field names in a string.
 * ------------------------------------------------------------------------- */
{
  int n = 0;

  while (isalnum (s[n])) n++;
  
  return n;
}


static void error (const char *s)
/* ------------------------------------------------------------------------- *
 * Print an error message and die.
 * ------------------------------------------------------------------------- */
{
  fprintf (stderr, "%s: %s\n", prog, s);
  exit (EXIT_FAILURE);
}


static void getargs (int     argc  ,
		     char**  argv  ,
		     FILE**  fp_in ,
		     FILE**  fp_out,
		     FORMAT* outf  ,
		     int*    ndump ,
		     int*    zero  )
/* ------------------------------------------------------------------------- *
 * Parse command line arguments.
 * ------------------------------------------------------------------------- */
{
  char c;
  int  i;
  char fname[FILENAME_MAX];

  while (--argc && (*++argv)[0] == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fputs (usage, stderr);
      exit  (EXIT_SUCCESS);
      break;

    case 'a':
      *outf = ASCII;
      break;
      
    case 'b':
      i = iformat ();
      switch (i) {
      case 1:
	*outf = IEEE_LITTLE;
	break;
      case 0:
	*outf = IEEE_BIG;
	break;
      default:
	fprintf (stderr, "%s: non-IEEE internal storage -- fix me", prog);
	exit (EXIT_FAILURE);
	break;
      }
      break;

    case 'n':
      if (*++argv[0])
	*ndump = atoi (*argv);
      else {
	*ndump = atoi (*++argv);
	argc--;
      }
      break;
      
    case 's':
      i = iformat ();
      switch (i) {
      case 1:
	*outf = IEEE_BIG;
	break;
      case 0:
	*outf = IEEE_LITTLE;
	break;
      default:
	fprintf (stderr, "%s: non-IEEE internal storage -- fix me", prog);
	exit (EXIT_FAILURE);
	break;
      }
      break;

    case 'v':
      verbose = 1;
      break;

    case 'o':
      if (*++argv[0])
	strcpy (fname, *argv);
      else {
	strcpy (fname, *++argv);
	argc--;
      }
      if ((*fp_out = fopen (fname,"w")) == (FILE*) NULL) {
	fprintf (stderr, "convert: unable to open output file: %s\n", fname);
	exit (EXIT_FAILURE);
      }
      *argv += strlen (*argv) - 1;
      break;

    case 'z':
      *zero = 1;		/* -- Cryptic, huh? */
      break;

    default:
      fprintf (stderr, "%s: unknown option -- %c\n", prog, c);
      break;
    }

  if (argc == 1)
    if ((*fp_in = fopen (*argv, "r")) == (FILE*) NULL) {
      sprintf (fname, "%s.fld", *argv);
      if ((*fp_in = fopen (fname, "r")) == (FILE*) NULL) {
	fprintf(stderr, "%s: unable to open input file -- %s or %s\n",
		prog, *argv, fname);
	exit (EXIT_FAILURE);
      }
    }

  return;
}


static void dswap (const int n,
		   double*   x)
/* ------------------------------------------------------------------------- *
 * Byte-reversal routine.
 * ------------------------------------------------------------------------- */
{
  register double *cx = x;
  register char   *c  = (char*) x;
  register int    i,j;

  for (i = 0; i < n; i++, cx++, c = (char*) cx)
     for (j = 0; j < 4; j++){
      char d = c[j];
      c[j]   = c[7-j];
      c[7-j] = d;
    }

  return;
}


static int iformat (void)
/* ------------------------------------------------------------------------- *
 * Return 1 if machine floating-point format is IEEE little-endian,
 * 0 if IEEE big-endian, -1 for unrecognized format.
 * ------------------------------------------------------------------------- */
{
  union { float f;  int i;    unsigned char c[4]; } v;
  union { double d; int i[2]; unsigned char c[8]; } u;
  int   reverse = (-1);
  u.d = 3;
  v.i = 3;
  if      (u.c[0] == 64 && u.c[1] == 8 && v.c[3] == 3) reverse = 0;
  else if (u.c[7] == 64 && u.c[6] == 8 && v.c[0] == 3) reverse = 1;

  return (reverse);
}


static FORMAT architecture (void)
/* ------------------------------------------------------------------------- *
 * What is this machine?  Die if unrecognized.
 * ------------------------------------------------------------------------- */
{
  switch (iformat ()) {
  case 1:
    return IEEE_LITTLE;
    break;
  case 0:
    return IEEE_BIG;
    break;
  case -1: default:
    fprintf (stderr, "%s: unrecognized machine architecture", prog);
    exit (EXIT_FAILURE);
    break;
  }

  return UNKNOWN;		/* Never happen. */
}


static FORMAT classify (const char* s)
/* ------------------------------------------------------------------------- *
 * Figure out what the input format is.
 * If ASCII, no problem.
 * Otherwise, we have to cope with backwards compatibility:
 *   if plain binary, set to machine's default binary,
 *   else set to declaration.
 * ------------------------------------------------------------------------- */
{
  const char* c = s;
  while (isspace (*c)) c++;

  switch (*c) {

  case 'a': case 'A':
    return ASCII;
    break;

  case 'b': case 'B':
    if      (!strstr (s, "IEEE"))   return architecture();
    else if ( strstr (s, "little")) return IEEE_LITTLE;
    else if ( strstr (s, "big"))    return IEEE_BIG;
    break;

  default:
    fprintf (stderr, "%s: unknown format specifier: %s\n", prog, s);
    exit (EXIT_FAILURE);
  }

  return ASCII;
}


static void get_data (FILE*        fp      ,
		      const int    selected,
		      const FORMAT format  ,
		      const FORMAT machine ,
		      const int    npts    ,
		      const int    nfields ,
		      double**     data    )
/* ------------------------------------------------------------------------- *
 * Read into data according to signalled format of input stream.
 * 
 * If this is a binary read, we can just fseek over data if this is
 * not the selected set.
 * ------------------------------------------------------------------------- */
{
  register int i, j;
  char         err[FILENAME_MAX];
  const int    swap = format != machine;

  switch (format) {

  case ASCII:
    
    if (verbose) fprintf (stderr, "(ASCII --> ");
    for (j = 0; j < npts; j++) {
      for (i = 0; i < nfields; i++)
	if (fscanf (fp, "%lf", data[i] + j) != 1) {
	  sprintf (err, "unable to read number: line %d, field %d\n", j+1,i+1);
	  error (err);
	}
      fgets (err, FILENAME_MAX, fp);
    }
    break;
      
  case IEEE_BIG:
    if (verbose) fprintf (stderr, "(IEEE-BIG_ENDIAN --> ");
    goto READ_BINARY;

  case IEEE_LITTLE:
    if (verbose) fprintf (stderr, "(IEEE-LITTLE_ENDIAN --> ");

  READ_BINARY:
    if (selected) {
      for (i = 0; i < nfields; i++)
	if (fread (data[i], sizeof (double), npts, fp) != npts) {
	  sprintf (err,"%s: an error has occured while reading (binary)",prog);
	  error   (err);
	}
      if (swap) for (i = 0; i < nfields; i++) dswap (npts, data[i]);
    } else {
      if (fseek (fp, sizeof(double)*npts*nfields, 1)) {
	sprintf (err,"%s: an error has occured while seeking (binary)",prog);
	error   (err);
      }
    }
    break;

  default: case UNKNOWN:
    sprintf (err, "%s: unknown input format", prog);
    error (err);
    break;

  }
}


static void put_data (FILE*        fp     ,
		      const FORMAT inputF ,
		      const FORMAT outputF,
		      const int    npts   ,
		      const int    nfields,
		      double**     data   )
/* ------------------------------------------------------------------------- *
 * Write from data to fp, according to input and desired output formats.
 * ------------------------------------------------------------------------- */
{
  char         err[FILENAME_MAX];
  register int i, j;
  const int    swap = outputF != architecture ();

  switch (outputF) {

  case ASCII:
    if (verbose) fprintf (stderr, "ASCII)\n");

    fprintf (fp, "ASCII                     Format\n");

    for (j = 0; j < npts; j++) {
      for (i = 0; i < nfields; i++)
	if (fprintf (fp, "%#16.10g ", data[i][j]) < 0) {
	  sprintf (err, "error occured while writing (ASCII)");
	  error (err);
	}
      fputs ("\n", fp);
    }

    break;
      
  case IEEE_BIG:
    if (verbose) fprintf (stderr, "IEEE-BIG_ENDIAN)\n"); 
    fprintf (fp, "binary IEEE big-endian    Format\n");
    goto WRITE_BINARY;

  case IEEE_LITTLE:
   if (verbose) fprintf (stderr, "IEEE-LITTLE_ENDIAN)\n"); 
   fprintf (fp, "binary IEEE little-endian Format\n");

  WRITE_BINARY:
    
    if (swap) for (i = 0; i < nfields; i++) dswap (npts, data[i]);
    for (i = 0; i < nfields; i++)
      if (fwrite (data[i], sizeof (double), npts, fp) != npts) {
	sprintf (err, "an error has occured while writing (binary)");
	error (err);
      }
    break;

  default: case UNKNOWN:
    sprintf (err, "%s: unknown input format", prog);
    error (err);
    break;

  }  
}
