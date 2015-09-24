/*****************************************************************************
 * sem2vtk: convert a SEM field file to VTK format.
 *
 * Usage:
 * sem2vtk [-h] [-c] [-o output] [-m mesh] [-n #] [-d #] [-w] input[.fld]
 *
 * Based on code sem2tec by Ron Henderson.
 * Modified by Stefan Schmidt, Chris Cantwell and Thomas Albrecht.
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

static char RCS[] = "$Id: sem2vtk.c,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <fcntl.h>
#include <cfemdef.h>
#include <cveclib.h>
#include <cfemlib.h>

#define MAXFIELDS 32

static char usage[] =
  "usage: sem2vtk [options] session[.fld]\n"
  "options:\n"
  "-h       ... print this message\n"
  "-o file  ... write output to the named file instead of running preplot\n"
  "-c       ... if nz > 1, generate cylindrical output\n"
  "-m file  ... read the mesh from the named file (instead of stdin)\n"
  "-d <num> ... extract dump <num> from file\n"
  "-n <num> ... evaluate the solution on an evenly-spaced mesh with N X N\n"
  "             points.  If N = 0, then no interpolation is done, i.e., the\n"
  "             output mesh will be on a standard GLL-spectral element mesh\n"
  "-w       ... extend the data by one additional plane in the z-direction\n";

static FILE    *fp_fld = 0,          /* default input files */
               *fp_msh = 0;
static char    *vtkfile;             /* output file name */

static int     nr, ns, nz, nel, nfields;
static int     nzp = 0, preplot_it = 1, np = 1, dump = 1, cylindrical=0;
static char    type[MAXFIELDS];
static double  *data[MAXFIELDS], *x, *y, *z;

static void    write_vtk   (FILE*);
static void    parse_args  (int, char**);
static void    read_mesh   (FILE*);
static void    read_data   (FILE*);
static void    interpolate (void);
static void    wrap        (void);
static double* do_interp   (const double*, const double*,
			    const double*, const double*, double*);
#include <errno.h>


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
  strcpy(fname, "tmp.XXXXXX");

  if ((fp = fdopen(mkstemp(fname), "w+")) == (FILE *) NULL) {
    fprintf (stderr, "sem2vtk: unable to open a temporary file\n");
    exit    (EXIT_FAILURE);
  }

  parse_args  (argc, argv);

  read_mesh   (fp_msh);
  while (dump--) read_data (fp_fld);
  interpolate ();
  wrap        ();
  write_vtk   (fp);

  if (preplot_it) {
    sprintf (buf, "preplot %s %s > /dev/null", fname, vtkfile);
    system  (buf);
    remove  (fname);
  } else {
    rewind (fp);
    fp_tec = fopen (vtkfile, "w");
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
	  fprintf(stderr, "sem2vtk: unable to open the mesh file -- %s\n",
		  fname);
	  exit(EXIT_FAILURE);
	}
	argv[0] += strlen(*argv)-1; argc--;
	break;
      case 'o':
	vtkfile = (char*) malloc(STR_MAX);
	strcpy(vtkfile, *++argv);
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
      case 'c':
        cylindrical = 1;
        break;
      default:
	fprintf(stderr, "sem2vtk: unknown option -- %c\n", c);
	break;
      }

  if (argc != 1) {
    fputs(usage, stderr);
    exit(EXIT_FAILURE);
  }

  /* open the input file */

  if ((fp_fld = fopen(*argv, "r")) == (FILE*) NULL) {
    sprintf(fname, "%s.fld", *argv);
    if ((fp_fld = fopen(fname, "r")) == (FILE*) NULL) {
      fprintf(stderr, "sem2vtk: unable to open %s or %s\n", *argv, fname);
      exit(EXIT_FAILURE);
    }
  }

  /* get the name of the ouput file (if not supplied) */

  if (!vtkfile) {
    int len = strlen(*argv);
    vtkfile = (char*) calloc (STR_MAX, sizeof(int));
    if   (strcmp(*argv + len-4, ".fld") == 0) strncpy(vtkfile, *argv, len-4);
    else                                      strcpy (vtkfile, *argv);
    strcat (vtkfile, ".plt");
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
 * Read in data files.  If NZ > 1, it is assumed that data are in the
 * file in plane-by-plane order.  For each plane, the ordering of data
 * varies depending on whether the file is in ASCII or binary format:
 * for ASCII, the fields are in column order, element-by-element
 * (row-major), whereas for binary formats the fields are written in
 * the file sequentially.
 *
 * Automatic conversion between little- and big-endian binary formats.
 * ------------------------------------------------------------------------- */
{
  int  i, m, n, nplane;
  int  nr_chk, ns_chk, nz_chk, nel_chk;
  char buf[STR_MAX], *c;

  /* -- Read the header down to the field list, check size of input. */

  for (n = 0; n < 3; n++) fgets(buf, STR_MAX, fp);

  if (sscanf (buf, "%d%d%d%d", &nr_chk, &ns_chk, &nz_chk, &nel_chk) != 4) {
    fputs ("error while reading mesh\n", stderr);
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
  while (n++ < 25 && nfields < MAXFIELDS)
    if (isalnum(*c++)) type[nfields++] = *(c-1);

  if (nfields > MAXFIELDS) {
    fprintf(stderr, "sem2vtk: a maximum of %d fields may be converted.\n",
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
	    fputs("sem2vtk: field file (ASCII) read error\n", stderr);
	    exit (EXIT_FAILURE);
	  }
    break;

  case 'b': {
    int swab, machine  = iformat();

    swab = (strstr (buf, "little") && machine == 0 ||
	    strstr (buf, "big"   ) && machine == 1  ) ? 1 : 0;

    for (n = 0; n < nfields; n++) {
      if (fread (data[n], sizeof(double), nz * nplane, fp) != nz * nplane) {
	fputs("sem2vtk: field file (binary) read error\n", stderr);
	exit (EXIT_FAILURE);
      }
      if (swab) dbrev (nz * nplane, data[n], 1, data[n], 1);
    }
    break;
  }

  default:
    fprintf (stderr, "sem2vtk: unknown format flag: '%c'\n", *c);
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


static void write_vtk (FILE *fp)
/* ------------------------------------------------------------------------- *
 * Write VTK file in ASCII format.
 * ------------------------------------------------------------------------- */
{
  register int i, j, k, m;
  const int    nrns = nr * ns, nplane = nr * ns * nel;
  /* write an extra plane of cells to link up last plane with first */
  const int    nzpc = (cylindrical) ? nzp : nzp-1;
  unsigned int offset1 = 0, offset2 = 0;

  /* write the VTK header--  could customise the descriptive string */

  fprintf (fp, "# vtk DataFile Version 2.0\n");
  fprintf (fp, "Semtex output file\n");
  fprintf (fp, "ASCII\n");

  /* define grid */
  fprintf (fp, "DATASET UNSTRUCTURED_GRID\n");
  fprintf (fp, "POINTS %i float\n", nel*nzp*nrns);
  for (m = 0; m < nzp; m++) {	   /* go over z planes nzp=1 for 2D      */
    for (k = 0; k < nel; k++) {	   /* go over elements                   */
      for (i = 0; i < nrns; i++) { /* go over x-y points in each element */

	/* write out x,y coordinate, adjusting if cylindrical coords */

	if (cylindrical)
	  fprintf(fp, "%#14.7g %#14.7g ", x[k*nrns + i],
		  y[k*nrns + i] * sin(z[m%nzp]));
	else
	  fprintf(fp, "%#14.7g %#14.7g ", x[k*nrns + i],
		  y[k*nrns + i]);
	/* write out z coordinate if we're in 3D, again adjusting
	   for cylindrical coordinates if necessary */

	if (z) {
	  if (cylindrical)
	    fprintf(fp, "%#14.7g ", y[k*nrns+i]*cos(z[m%nzp]));
	  else
	    fprintf(fp, "%#14.7g ", z[m]);
	}
	else
	  fprintf(fp, "         0.0 ");
	/* done writing this point, new line */
	fprintf(fp,"\n");
      }
    }
  }
  fprintf(fp, "\n");

  /* define cells */
  if (z)
    fprintf (fp, "CELLS %i %i\n",
	     nzpc*nel*(nr-1)*(ns-1), nzpc*nel*(nr-1)*(ns-1)*9);
  else
    fprintf (fp, "CELLS %i %i\n", nel*(nr-1)*(ns-1), nel*(nr-1)*(ns-1)*5);
  for (m = 0; m < nzp; m++) {
    if (!cylindrical && nzp > 1 && m == nzp - 1) break;
    for (k = 0; k < nel; k++) {
      /* compute plane/element offsets */
      offset1 = m*nplane + k*nrns;
      /* if cylindrical, last plane is the first plane to link up */
      offset2 = (cylindrical) ? (m+1)%nzp*nplane : (m+1)*nplane;
      offset2 += k*nrns;
      for (j = 0; j < nr - 1; j++) {
	for (i = 0; i < ns - 1; i++) {
	  if (z) fprintf(fp, "8 ");
	  else fprintf(fp, "4 ");
	  fprintf(fp, "%i ", offset1 + j*ns + i);
	  fprintf(fp, "%i ", offset1 + j*ns + i + 1);
	  fprintf(fp, "%i ", offset1 + (j+1)*ns + i + 1);
	  fprintf(fp, "%i ", offset1 + (j+1)*ns + i);
	  if (z) {
	    fprintf(fp, "%i ", offset2 + j*ns + i);
	    fprintf(fp, "%i ", offset2 + j*ns + i + 1);
	    fprintf(fp, "%i ", offset2 + (j+1)*ns + i + 1);
	    fprintf(fp, "%i ", offset2 + (j+1)*ns + i);
	  }
	  fprintf(fp, "\n");
	}
      }
    }
  }
  fprintf(fp, "\n");

  /* define cell types */
  fprintf (fp, "CELL_TYPES %i\n", nzp*nel*(nr-1)*(ns-1));
  for (i = 0; i < nzp*nel*(nr-1)*(ns-1); i++) {
    if (z) fprintf( fp, "12\n");
    else fprintf( fp, "9\n");
  }
  fprintf(fp, "\n");

  /* now write out the data for each point */
  fprintf(fp, "POINT_DATA %i\n", nel*nzp*nrns);

  /* print out vector velocity values */
  j = 0;
  if (nfields > 1) {
    fprintf(fp, "VECTORS Velocity float\n");
    for (m = 0; m < nzp; m++) {
      for (k = 0; k < nel; k++) {
	/* plane/element offset */
	offset1 = (cylindrical) ? (m%nzp)*nplane : m*nplane;
	offset1 += k*nrns;
	for (i = 0; i < nrns; i++) {
	  /* u, v values are in 1st and 2nd data array */
	  fprintf(fp, "%#14.7g ", data[0][offset1 + i]);
	  fprintf(fp, "%#14.7g ", data[1][offset1 + i]);
	  if (type[2] == 'w') fprintf(fp, "%14.7g ", data[2][offset1 + i]);
	  else fprintf(fp, "0.0 ");
	  fprintf(fp, "\n");
	}
      }
    }
    if (type[2] == 'w') j = 3; else j = 2;
  }

  /* print out scalar values */

  for (; j < nfields; j++) {
    fprintf(fp, "SCALARS %c float\n", type[j]);
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (m = 0; m < nzp; m++) {
      for (k = 0; k < nel; k++) {
	/* plane/element offset */
	offset1 = (cylindrical) ? (m%nzp)*nplane : m*nplane;
	offset1 += k*nrns;
	for (i = 0; i < nrns; i++) {
	  /* pressure is in 3rd data array */
	  fprintf(fp, "%#14.7g\n", data[j][offset1 + i]);
	}
      }
    }
    fprintf(fp, "\n");
  }

  for (i = 0; i < nfields; i++)
    free (data[i]);
  fflush (fp);
}
