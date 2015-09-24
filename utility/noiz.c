/*****************************************************************************
 * noiz.c: add a random Gaussian perturbation to a velocity field.
 *         Optionally filter out a named mode.
 *
 * Copyright (c) 1996 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn
 *
 * USAGE
 * -----
 * noiz [-h] [-f] [-o output] [-p perturb] [-m mode] [-s seed] [input[.fld]
 *
 * SYNOPSIS
 * --------
 * Noiz reads a field file and adds a gaussian-distributed random variable
 * of specified standard deviation to each velocity datum.  Fields may be in
 * ASCII or binary format, output is in same format.  Optionally, noise
 * is added just to a prescribed Fourier mode (mode numbers begin at zero).
 *
 * NOTES
 * -----
 * Default value of perturbation is 0.0.
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

static char RCS[] = "$Id: noiz.c,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include <cfemlib.h>
#include <cfemdef.h>
#include <cveclib.h>

#define IA    16807
#define IM    2147483647
#define AM    (1.0/IM)
#define IQ    127773
#define IR    2836
#define NTAB  32
#define NDIV  (1+(IM-1)/NTAB)
#define RNMX  (1.0-EPSDP)
#define UNSET -1
static long   seed = 0;

static void   getargs (int, char**, FILE**, FILE**, double*, int*, int*);
static void   a_to_a  (int, int, int, int,
		       FILE*, FILE*, char*, double, int, int);
static void   b_to_b  (int, int, int, int,
		       FILE*, FILE*, char*, double, int, int, int);
static double gasdev  (long*);   
static void   perturb (double*, const int, const int, const int, const double);
static void   filter  (double*, const int, const int, const int);


static char prog[] = "noiz";
static const char *hdr_fmt[] = { 
  "%-25s "    "Session\n",
  "%-25s "    "Created\n",
  "%-25s "    "Nr, Ns, Nz, Elements\n",
  "%-25d "    "Step\n",
  "%-25.6g "  "Time\n",
  "%-25.6g "  "Time step\n",
  "%-25.6g "  "Kinvis\n",
  "%-25.6g "  "Beta\n",
  "%-25s "    "Fields written\n",
  "%-25s "    "Format\n"
};


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Wrapper.
 * ------------------------------------------------------------------------- */
{
  char   buf[STR_MAX], fields[STR_MAX], fmt[STR_MAX], *c;
  int    i, n, np, nz, nel, mode = UNSET, swab = 0, filt = 0;
  FILE   *fp_in  = stdin,
         *fp_out = stdout;
  double pert    = 0.0;

  getargs (argc, argv, &fp_in, &fp_out, &pert, &mode, &filt);

  format (fmt);

  while (fgets (buf, STR_MAX, fp_in)) { 

    n = 3;
    while (--n) {
      fputs (buf, fp_out);
      fgets (buf, STR_MAX, fp_in);
    }

    if (sscanf (buf, "%d%*s%d%d", &np, &nz, &nel) != 3)
      message (prog, "unable to read the file size", ERROR);

    if (mode != UNSET && 2 * mode > nz) {
      sprintf (fields, "too many modes (%1d) for input (nz = %1d)", mode, nz);
      message (prog, fields, ERROR);
    }

    n = 6;
    while (--n) {
      fputs (buf, fp_out);
      fgets (buf, STR_MAX, fp_in);
    }
    fputs(buf, fp_out);   

    fgets(fields, STR_MAX, fp_in);
    
    for (n = 0, i = 0; i < 25; i++) if (isalnum(fields[i])) n++;

    fputs (fields, fp_out);
    fgets (buf, STR_MAX, fp_in);
    for (i = 0; i < strlen (buf); i++) buf[i] = tolower (buf[i]);

    switch (buf[0]) {
    case 'a':
      fprintf (fp_out, hdr_fmt[9], "ASCII");
      a_to_a  (np, nz, nel, n, fp_in, fp_out, fields, pert, mode, filt);
      break;

    case 'b':
      if (!strstr (buf, "endian"))
	message (prog, "input field file in unknown binary format", WARNING);
      else {
	swab = (   (strstr (buf, "big") && strstr (fmt, "little"))
		|| (strstr (fmt, "big") && strstr (buf, "little")) );
      }
      sprintf (buf, "binary ");
      strcat  (buf, fmt);
      fprintf (fp_out, hdr_fmt[9], buf);
      b_to_b  (np, nz, nel, n, fp_in, fp_out, fields, pert, mode, swab, filt);
      break;

    default:
      sprintf (buf, "unknown format flag -- %c", buf[0]);
      message (prog, buf, ERROR);
      break;
    }
  } 
  
  return EXIT_SUCCESS;
}


static void a_to_a (int    np     ,
		    int    nz     ,
		    int    nel    ,
		    int    nfields, 
		    FILE*  in     , 
		    FILE*  out    , 
		    char*  fields , 
		    double pert   ,
		    int    mode   ,
		    int    filt   )
/* ------------------------------------------------------------------------- *
 * ASCII input.
 * ------------------------------------------------------------------------- */
{
  register int i, j, kr, ki;
  const int    nplane = np * np * nel,
               npts   = nz * nplane,
               ntot   = nfields * npts;
  double**     data;
  char         buf[STR_MAX];
  double       datum;

  data = dmatrix (0, nfields - 1, 0, npts - 1);

  for (j = 0; j < npts; j++)
    for (i = 0; i < nfields; i++)
      if (fscanf (in, "%lf", &datum) != 1) {
	sprintf (buf, "unable to read a number -- line %d, field %d\n",
		j+1, i+1);
	message (prog, buf, ERROR);
      }

  for (i = 0; i < nfields; i++) {
    switch (fields[i]) {
    case 'u': case 'v': case 'w': 
      if (pert > 0.0) perturb (data[i], mode, nz, nplane, pert);
    default:
      if (filt)       filter  (data[i], mode, nz, nplane);
      break;
    }
  }

  for (j = 0; j < npts; j++) {
    for (i = 0; i < nfields; i++)
      if (fprintf (out, "%#16.10g ", datum) < 1)
	message (prog, "an error has occured while writing", ERROR);
    fprintf (out, "\n");
  }

  freeDmatrix (data, 0, 0);
}


static void b_to_b (int    np     ,
		    int    nz     ,
		    int    nel    ,
		    int    nfields, 
		    FILE*  in     , 
		    FILE*  out    , 
		    char*  fields , 
		    double pert   ,
		    int    mode   ,
		    int    swab   ,
		    int    filt   )
/* ------------------------------------------------------------------------- *
 * Binary input.
 * ------------------------------------------------------------------------- */
{
  register int i, j, kr, ki;
  const int    nplane = np * np * nel,
               npts   = nz * nplane,
               ntot   = nfields * npts;
  double**     data;

  data = dmatrix (0, nfields - 1, 0, npts - 1);

  for (i = 0; i < nfields; i++) {

    if (fread (data[i], sizeof (double), npts, in) != npts)
      message (prog, "an error has occured while reading", ERROR);
    if (swab) dbrev (npts, data[i], 1, data[i], 1);

    switch (fields[i]) {
    case 'u': case 'v': case 'w': 
      if (pert > 0.0) perturb (data[i], mode, nz, nplane, pert);
    default:
      if (filt)       filter  (data[i], mode, nz, nplane);
      break;
    }
  }

  if (fwrite (data[0], sizeof (double), ntot, out) != ntot)
    message (prog, "an error has occured while writing", ERROR);
  
  freeDmatrix (data, 0, 0);
}


static void getargs (int     argc  ,
		     char**  argv  ,
		     FILE**  fp_in ,
		     FILE**  fp_out,
		     double* pert  ,     
		     int*    mode  ,
		     int*    filter)
/* ------------------------------------------------------------------------- *
 * Parse command line arguments.
 * ------------------------------------------------------------------------- */
{
  char c;
  int  i;
  char fname[FILENAME_MAX];
  char usage[] = "usage: noiz [options] [input[.fld]]\n"
    "options:\n"
    "-h         ... print this help message\n"
    "-f         ... filter instead of perturb\n"
    "-o output  ... write to named file\n"
    "-p perturb ... standard deviation of perturbation  [Default: 0.0]\n"
    "-m mode    ... add noise only to this Fourier mode [Default: all modes]\n"
    "-s seed    ... set random number seed              [Default: 0]\n";

  while (--argc && (*++argv)[0] == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fputs (usage, stderr);
      exit  (EXIT_SUCCESS);
      break;
    case 'f':
      *filter = 1;
      break;
    case 'o':
      if (*++argv[0]) strcpy(fname, *argv);
      else {strcpy(fname, *++argv); argc--;}
      if ((*fp_out = fopen(fname,"w")) == (FILE*) NULL) {
	fprintf(stderr, "%s: unable to open the output file -- %s\n", 
		prog, fname);
	exit (EXIT_FAILURE);
      }
      *argv += strlen (*argv)-1;
      break;
    case 'p':
      if (*++argv[0]) *pert = atof (*argv);
      else {*pert = atof (*++argv); argc--;}
      *argv += strlen (*argv)-1;
      break;
    case 'm':
      if (*++argv[0]) *mode = atoi (*argv);
      else {*mode = atoi (*++argv); argc--;}
      break;
    case 's':
      if (*++argv[0]) seed = atoi (*argv);
      else {seed = atoi (*++argv); argc--;}
      break;
    default:
      fprintf (stderr, "%s: unknown option -- %c\n", prog, c);
      break;
    }
  
  if (argc == 1)
    if ((*fp_in = fopen(*argv, "r")) == (FILE*) NULL) {
      sprintf(fname, "%s.fld", *argv);
      if ((*fp_in = fopen(fname, "r")) == (FILE*) NULL) {
	fprintf(stderr, "%s: unable to open input file -- %s or %s\n",
		prog, *argv, fname);
	exit (EXIT_FAILURE);
      }
    }

  return;
}


static double ran1 (long *idum)
/* ------------------------------------------------------------------------- *
 * Generate IUD random variates on (0, 1).  Numerical Recipes.
 * ------------------------------------------------------------------------- */
{
  register int j;
  long         k;
  static long  iy = 0;
  static long  iv[NTAB];
  double       temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum = 1;
    else *idum = -(*idum);
    for (j=NTAB+7; j>=0; j--) {
      k = (*idum)/IQ;
      *idum = IA * (*idum - k*IQ) - IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j]= *idum;
    }
    iy = iv[0];
  }
  k = (*idum) / IQ;
  *idum = IA * (*idum - k*IQ) - IR*k;
  if (*idum < 0) *idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if   ((temp = AM * iy) > RNMX) return RNMX;
  else                           return temp;
}
 

static double gasdev (long *idum)
/* ------------------------------------------------------------------------- *
 * Generate normally distributed deviate with zero mean & unit variance.
 * Numerical Recipes.
 * ------------------------------------------------------------------------- */
{
  static int    iset = 0;
  static double gset;
  double        fac, r, v1, v2;
    
  if  (iset == 0) {
    do {
      v1 = 2.0 * ran1 (idum) - 1.0;
      v2 = 2.0 * ran1 (idum) - 1.0;
      r  = v1 * v1 + v2 * v2;
    } while (r >= 1.0);
    fac = sqrt (-2.0 * log (r) / r);
    gset = v1 * fac;
    iset = 1;
    return v2 * fac;
  } else {
    iset = 0;
    return gset;
  }
}


static void perturb (double*      data  ,
		     const int    mode  ,
		     const int    nz    , 
		     const int    nplane,
		     const double pert  )
/* ------------------------------------------------------------------------- *
 * Add perturbation to data field.
 * ------------------------------------------------------------------------- */
{
  register int j;
  const int    npts = nz * nplane;
  const int    kr   = (2 * mode)     * nplane;
  const int    ki   = (2 * mode + 1) * nplane;
  double       eps;

  if (mode != UNSET) {	/* -- Perturb only specified Fourier mode. */

    dDFTr (data, nz, nplane, +1);
    
    eps = pert * nz;		/* -- Account for scaling of modes. */

    if      (mode == 0)
      for (j = 0; j < nplane; j++) data[j] += eps * gasdev (&seed);

    else if (mode == (nz >> 1))
      for (j = 0; j < nplane; j++) data[nplane + j] += eps * gasdev (&seed);

    else {
      for (j = 0; j < nplane; j++) data[kr + j] += eps * gasdev (&seed);
      for (j = 0; j < nplane; j++) data[ki + j] += eps * gasdev (&seed);
    }
    
    dDFTr (data, nz, nplane, -1);

  } else {			/* -- Perturb all modes. */
    
    eps = pert;
    for (j = 0; j < npts; j++) data[j] += eps * gasdev (&seed);
    
  }
}


static void filter (double*   data  ,
		    const int mode  ,
		    const int nz    , 
		    const int nplane)
/* ------------------------------------------------------------------------- *
 * Zero data in named mode.
 * ------------------------------------------------------------------------- */
{
  if (mode != UNSET) {

    dDFTr (data, nz, nplane, +1);

    if      (mode == 0)         dzero (nplane, data, 1);
    else if (mode == (nz >> 1)) dzero (nplane, data + nplane, 1);
    else                        dzero (2 * nplane, data + 2*mode*nplane, 1);
    
    dDFTr (data, nz, nplane, -1);
  }
}
