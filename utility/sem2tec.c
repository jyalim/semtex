/*****************************************************************************
 * sem2tec: convert a semtex field file to AMTEC Tecplot format.
 *
 * Copyright (c) 1990 <--> $Date: 2015/04/20 11:14:19 $, 
 *   Ron Henderson, Hugh Blackburn
 *
 * Usage: sem2tec [-h] [-o output] [-m mesh] [-n #] [-d #] [-w] input[.fld]
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

static char RCS[] = "$Id: sem2tec.c,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "cfemdef.h"
#include "cveclib.h"
#include "cfemlib.h"

#define MAXFIELDS 50

static char usage[] = 
  "usage: sem2tec [options] session[.fld]\n"
  "options:\n"
  "-h       ... print this message\n"
  "-o file  ... write output to the named file instead of running preplot\n" 
  "-m file  ... read the mesh from the named file (instead of stdin)\n"
  "-d <num> ... extract dump <num> from file\n"
  "-n <num> ... evaluate the solution on an evenly-spaced mesh with N X N\n"
  "             points.  If N = 0, then no interpolation is done, i.e., the\n"
  "             output mesh will be on a standard GLL-spectral element mesh\n"
  "-w       ... extend the data by one additional plane in the z-direction\n";

static FILE    *fp_fld = 0,          /* default input files */
               *fp_msh = 0;
static char    *tecfile;             /* output file name */

static int     nr, ns, nz, nel, nfields;
static int     nzp = 0, preplot_it = 1, np = 1, dump = 1;
static char    type[MAXFIELDS];
static double  *data[MAXFIELDS], *x, *y, *z;

static void    write_tec   (FILE*);
static void    parse_args  (int, char**);
static void    read_mesh   (FILE*);
static void    read_data   (FILE*);
static void    interpolate (void);
static void    wrap        (void);
static double* do_interp   (const double*, const double*, 
			    const double*, const double*, double*);


int main (int    argc,
          char** argv)
/* ------------------------------------------------------------------------- *
 * Driver.
 * ------------------------------------------------------------------------- */
{
  char fname[STR_MAX];
  char buf  [STR_MAX];
  FILE *fp, *fp_tec;
  
  fp_msh = stdin;

  if ((fp = fopen (tmpnam (fname),"w+")) == (FILE*) NULL) {
    fprintf (stderr, "sem2tec: unable to open a temporary file\n");
    exit    (EXIT_FAILURE);
  }

  parse_args  (argc, argv);

  read_mesh   (fp_msh);
  while (dump--) read_data (fp_fld);
  interpolate ();
  wrap        ();
  write_tec   (fp);

  if (preplot_it) {
    sprintf (buf, "preplot %s %s > /dev/null", fname, tecfile);
    system  (buf);
    remove  (fname);
  } else {
    rewind (fp);
    fp_tec = fopen (tecfile, "w");
    while  (fgets(buf, STR_MAX, fp)) fputs(buf, fp_tec);
    fclose (fp_tec);
    fclose (fp);
  }

  return EXIT_SUCCESS;
}


static void parse_args (int    argc,
			char** argv)
/* ------------------------------------------------------------------------- *
 * Read command-line, open files.
 * ------------------------------------------------------------------------- */
{
  char c, fname[FILENAME_MAX];

  while (--argc && (*++argv)[0] == '-') 
    while (c = *++argv[0])
      switch (c) {
      case 'h':
	fputs (usage, stdout);
	exit  (EXIT_SUCCESS);
	break;
      case 'm':
	sprintf(fname, "%s", *++argv);
	if ((fp_msh = fopen(fname, "r")) == (FILE*) NULL) {
	  fprintf(stderr, "sem2tec: unable to open the mesh file -- %s\n", 
		  fname);
	  exit(EXIT_FAILURE);
	}
	argv[0] += strlen(*argv)-1; argc--;
	break;
      case 'o':
	tecfile = (char*) malloc(STR_MAX);
	strcpy(tecfile, *++argv);
	argv[0] += strlen(*argv)-1; argc--;
	preplot_it = 0;
	break;
      case 'd':
	if (*++argv[0])
	  dump = atoi(*argv);
	else {
	  dump = atoi(*++argv);
	  argc--;
	}
	(*argv)[1] = '\0';
	break;
      case 'n':
	if (*++argv[0])
	  np = atoi(*argv);
	else {
	  np = atoi(*++argv);
	  argc--;
	}
	(*argv)[1] = '\0';
	break;
      case 'w':
	nzp = 1;
	break;
      default:
	fprintf(stderr, "sem2tec: unknown option -- %c\n", c);
	break;
      }

  if (argc != 1) {
    fputs(usage, stderr);
    exit(EXIT_FAILURE);
  }

#if 0
  /* open the input file, must end in .fld */

  if (!strstr (*argv, ".fld")) sprintf (fname, "%s.fld", *argv);
  else                         sprintf (fname, "%s",     *argv);
  if ((fp_fld = fopen(fname, "r")) == (FILE*) NULL) {
    fprintf(stderr, "sem2tec: unable to open %s or %s\n", *argv, fname);
    exit(EXIT_FAILURE);
  }
#else
  /* open the input file */

  if ((fp_fld = fopen(*argv, "r")) == (FILE*) NULL) {
    sprintf(fname, "%s.fld", *argv);
    if ((fp_fld = fopen(fname, "r")) == (FILE*) NULL) {
      fprintf(stderr, "sem2tec: unable to open %s or %s\n", *argv, fname);
      exit(EXIT_FAILURE);
    }
  }
#endif

  /* get the name of the ouput file (if not supplied) */

  if (!tecfile) {
    int len = strlen(*argv);
    tecfile = (char*) calloc (STR_MAX, sizeof(int));
    if   (strcmp(*argv + len-4, ".fld") == 0) strncpy(tecfile, *argv, len-4);
    else                                      strcpy (tecfile, *argv);
    strcat (tecfile, ".plt");
  }
}


static void read_mesh (FILE *fp)
/* ------------------------------------------------------------------------- *
 * Read in mesh data, as generated by meshpr.
 * ------------------------------------------------------------------------- */
{
  int  n, i;
  char buf[STR_MAX];

  fgets (buf, STR_MAX, fp);
  if (sscanf (buf, "%d%d%d%d", &nr, &ns, &nz, &nel) != 4) {
    fputs ("error while reading mesh\n", stderr);
    exit  (EXIT_FAILURE);
  }

  n   = nr * ns * nel;
  x   = dvector (0, n - 1);
  y   = dvector (0, n - 1);
  z   = (nz > 1) ? dvector (0, nz) : 0;
  nzp = (z && nzp) ? nz + 1 : nz;

  for (i = 0; i < n; i++) {
    fgets (buf, STR_MAX, fp);
    if (sscanf (buf, "%lf%lf", x+i, y+i) != 2) {
      fputs ("error while reading mesh data\n", stderr);
      exit  (EXIT_FAILURE);
    }
  }
  
  if (z) 
    for (i = 0; i <= nz; i++) {
      fgets (buf, STR_MAX, fp);
      if (sscanf (buf, "%lf", z+i) != 1) {
	fputs ("error while reading z mesh data\n", stderr);
	exit  (EXIT_FAILURE);
      }
    }
}


static void read_data (FILE *fp)
/* ------------------------------------------------------------------------- *
 * Read in data files.  If NZ > 1, it is assumed that data are in the file
 * in plane-by-plane order.  For each plane, the ordering of data varies
 * depending on whether the file is in ASCII or binary format: for ASCII, the
 * fields are in column order, element-by-element (row-major), whereas for
 * binary formats the fields are written in the file sequentially.
 *
 * Automatic conversion between little- and big-endian binary formats.
 * ------------------------------------------------------------------------- */
{
  int  i, m, n, nplane;
  int  nr_chk, ns_chk, nz_chk, nel_chk;
  char buf[STR_MAX], *c;
  
  /* -- Read the header down to the field list, check size of input. */

  for (n = 0; n < 3; n++) {
    fgets (buf, STR_MAX, fp);
  }

  if (sscanf (buf, "%d%d%d%d", &nr_chk, &ns_chk, &nz_chk, &nel_chk) != 4) {
    fputs ("error reading size of field\n", stderr);
    exit  (EXIT_FAILURE);
  }

  if (nr != nr_chk || ns != ns_chk || nel != nel_chk) {
    fputs ("2D structure of mesh and field file do not agree\n", stderr);
    exit (EXIT_FAILURE);
  }

  for (n = 3; n < 9; n++) fgets(buf, STR_MAX, fp);

  /* -- Read the list of fields. */

  n       = 0;
  c       = buf;
  nfields = 0;
  while (isalpha(*c) && nfields < MAXFIELDS) type[nfields++] = (*c++);

  if (nfields > MAXFIELDS) {
    fprintf(stderr, "sem2tec: a maximum of %d fields may be converted.\n", 
	    MAXFIELDS);
    exit(EXIT_FAILURE);
  }

  /* -- Allocate memory. */

  nplane = nr * ns * nel;
  for (n = 0; n < nfields; n++)
    data[n] = (double*) malloc (nzp * nplane * sizeof (double));

  /* -- Check the format. */

  c = fgets(buf, STR_MAX, fp); 
  while (isspace(*c)) c++;

  switch (tolower(*c)) {                     /* ASCII or binary read? */

  case 'a':
    for (m = 0; m < nz; m++)
      for (i = 0; i < nplane; i++)
	for (n = 0; n < nfields; n++)
	  if (fscanf(fp, "%lf", data[n] + m * nplane + i) < 0) {
	    fputs("sem2tec: field file (ASCII) read error\n", stderr);
	    exit (EXIT_FAILURE);
	  }
    break;

  case 'b': {
    int swab, machine  = iformat();

    swab = (strstr (buf, "little") && machine == 0 ||
	    strstr (buf, "big"   ) && machine == 1  ) ? 1 : 0;

    for (n = 0; n < nfields; n++) {
      if (fread (data[n], sizeof(double), nz * nplane, fp) != nz * nplane) {
	fputs("sem2tec: field file (binary) read error\n", stderr);
	  exit (EXIT_FAILURE);
      }
      if (swab) dbrev (nz * nplane, data[n], 1, data[n], 1);
    }
    break;
  }

  default:
    fprintf (stderr, "sem2tec: unknown format flag: '%c'\n", *c);
    exit    (EXIT_FAILURE);
    break;
  }
}


static void interpolate (void)
/* ------------------------------------------------------------------------- *
 * Interpolate from the GLL mesh to an evenly-spaced mesh.
 * ------------------------------------------------------------------------- */
{
  register int k, m, nplane_new;
  const int    nplane_old = nr * ns * nel;
  const double *imr, *itmr, *ims, *itms;
  double       *mesh_x, *mesh_y;
  double       **newplane = (double**) malloc (nz * sizeof (double*));

  switch (np) {
  case 0:              /* interpolation turned off */
    return;
    break;
    
  case 1:              /* no size specified ... use (NR|NS) */
    np = MAX (nr, ns);
    break;

  default:             /* size specified on the command line */
    break;
  }
  
  nplane_new = np * np * nel;

  /* -- Compute interpolation matrices. */

  proj (&imr, &itmr, nr, GLJ, 0.0, 0.0, np, TRZ, 0.0, 0.0);
  proj (&ims, &itms, ns, GLJ, 0.0, 0.0, np, TRZ, 0.0, 0.0);

  /* -- Interpolate the mesh. */

  mesh_x = do_interp (imr, itmr, ims, itms, x);
  mesh_y = do_interp (imr, itmr, ims, itms, y);

  free (x); x = mesh_x;
  free (y); y = mesh_y;

  /* -- Interpolate data plane-by-plane. */

  for (k = 0; k < nfields; k++) {

    for (m = 0; m < nz; m++)
      newplane[m] = do_interp (imr, itmr, ims, itms, data[k] + m * nplane_old);

    free (data[k]);
    data[k] = (double*) malloc (nplane_new * nzp * sizeof (double));

    for (m = 0; m < nz; m++) {
      dcopy (nplane_new, newplane[m], 1, data[k] + m * nplane_new, 1);
      free  (newplane[m]);
    }
  }

  nr = ns = np;
}


static double* do_interp (const double* imr ,
			  const double* itmr,
			  const double* ims ,
			  const double* itms,
			  double*       data)
/* ------------------------------------------------------------------------- *
 * Wrapper for 2D tensor-product interpolation.
 * ------------------------------------------------------------------------- */
{
  register int k;
  const int    nrns = nr * ns;
  const int    ntot = np * np;
  double       *new  = dvector (0, ntot * nel - 1),
               *tmp  = dvector (0, np*nr),
               *p    = new;

  for (k = 0; k < nel; k++, data += nrns, p += ntot) {
    dmxm ((double*) ims, np, data, ns, tmp, nr);
    dmxm (tmp, np, (double*) itmr, nr, p  , np);
  }

  return new;
}


static void wrap (void)
/* ------------------------------------------------------------------------- *
 * Extend data in the (periodic) z-direction so that it wraps around.
 * ------------------------------------------------------------------------- */
{
  register int i;
  const int    nplane_new = nr * ns * nel;

  if (nzp == nz) return;

  for (i = 0; i < nfields; i++)
    dcopy (nplane_new, data[i], 1, data[i] + nz * nplane_new, 1);
}


static void write_tec (FILE *fp)
/* ------------------------------------------------------------------------- *
 * Write in ASCII to temporary file, tecplot ASCII format.
 * ------------------------------------------------------------------------- */
{
  register int i, j, k, m;
  const int    nrns = nr * ns, nplane = nr * ns * nel;

  fprintf (fp, "VARIABLES = \"X\" \"Y\" ");
  if (z) fprintf (fp, "\"Z\" ");
#if 0				/* Old code version. */
  for (i = 0; i < nfields; i++) fprintf (fp, "\"%c\" ", toupper(type[i]));
#else
  for (i = 0; i < nfields; i++) fprintf (fp, "\"%c\" ", type[i]);
#endif
  fprintf (fp, "\n");

  for (k = 0; k < nel; k++) {
    fprintf (fp, "ZONE T=\"Element %d\", I=%d, J=%d,", k+1, nr, ns);
    if (z) fprintf (fp, " K=%d,", nzp);
    fprintf (fp, " F=POINT\n");
    for (m = 0; m < nzp; m++) {
      for (i = 0; i < nrns; i++) {
	fprintf (fp, "%#14.7g %#14.7g ", x[k*nrns + i], y[k*nrns + i]);
	if (z) fprintf (fp, "%#14.7g ",  z[m]);
	for (j = 0; j < nfields; j++)
	  fprintf(fp, "%#14.7g ", data[j][m*nplane + k*nrns + i]);
	fprintf(fp, "\n");
      }
    }
  }

  for (i = 0; i < nfields; i++) free (data[i]);
  fflush (fp);
}
