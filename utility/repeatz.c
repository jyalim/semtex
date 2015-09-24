/*****************************************************************************
 * repeatz.c: read a field file, and output the number of repetitions
 * in the z direction indicated on the command line.  Field must be
 * binary format.  By default, output the original file, no repeats.
 * Adust Nz and Beta in header as appropriate.
 *
 * NB: Nz must be adjusted so it has prime factors of 2,3,5 to ensure
 * the resulting field can be Fourier transformed in z.  This 'feature'
 * can be defeated with command-line flag -f.
 *
 * Copyright (c) 2002 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn.
 *
 * USAGE
 * -----
 * repeatz [-h] [-n <rep>] [-f] [input[.fld]
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

static char RCS[] = "$Id: repeatz.c,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include <cfemlib.h>
#include <cfemdef.h>
#include <cveclib.h>

static void getargs (int, char**, FILE**, int*, int*);
static int  roundup (const int);
static void pack    (const double*, const int, double*,
		     const int, const int, const int);

static char prog[] = "repeatz";
static const char *hdr_fmt[] = { 
  "%-25s "              "Session\n",
  "%-25s "              "Created\n",
  "%-5d%-5d%-5d%-10d "  "Nr, Ns, Nz, Elements\n",
  "%-25d "              "Step\n",
  "%-25.6g "            "Time\n",
  "%-25.6g "            "Time step\n",
  "%-25.6g "            "Kinvis\n",
  "%-25.6g "            "Beta\n",
  "%-25s "              "Fields written\n",
  "%-25s "              "Format\n"
};


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Wrapper.
 * ------------------------------------------------------------------------- */
{
  char   buf[STR_MAX], fmt[STR_MAX];
  int    i, j, k, n, np, nzin, nzout, nel, nrep = 1, force = 0;
  int    nfields, nplane, nptin, nptout, ntot, swab;
  FILE   *fp_in = stdin, *fp_out = stdout;
  double *datain, *dataout, beta;

  getargs (argc, argv, &fp_in, &nrep, &force);
  format  (fmt);

  while (fgets (buf, STR_MAX, fp_in)) { 

    fputs (buf, fp_out); fgets (buf, STR_MAX, fp_in);
    fputs (buf, fp_out); fgets (buf, STR_MAX, fp_in);
    
    if (sscanf (buf, "%d%*s%d%d", &np, &nzin, &nel) != 3)
      message (prog, "unable to read the file size", ERROR);

    if (!force) {
      if ((nzout = roundup (nzin)) != nzin)
	message (prog, "input nz does not have 2, 3, 5 factors", ERROR);
      nzout = roundup (nzin * nrep);
    } else
      nzout = nzin * nrep;

    fprintf (fp_out, hdr_fmt[2], np, np, nzout, nel);

    n = 4;
    while (n--) { fgets (buf, STR_MAX, fp_in); fputs (buf, fp_out); }

    fgets (buf, STR_MAX, fp_in);
    sscanf (buf, "%lf", &beta);
    beta /= nrep;
    fprintf (fp_out, hdr_fmt[7], beta);   

    fgets (buf, STR_MAX, fp_in); fputs (buf, fp_out);
    for (nfields = 0, i = 0; i < 25; i++) if (isalnum(buf[i])) nfields++;

    fgets (buf, STR_MAX, fp_in);
    if (!strstr(buf, "binary"))
      message (prog, "input file not binary format", ERROR);
    swab = (strstr (buf, "big") && strstr (fmt, "little")) || 
           (strstr (fmt, "big") && strstr (buf, "little"));
    strcat (strcpy (buf, "binary "), fmt);
    fprintf (fp_out, hdr_fmt[9], buf);
    
    /* -- Set sizes, allocate storage. */

    nplane  = np * np * nel;
    ntot    = nplane + (nplane & 1);
    nptin   = nzin  * ntot;
    nptout  = nzout * ntot;
    datain  = dvector (0, nptin  - 1);
    dataout = dvector (0, nptout - 1);
    
    /* -- Read and write all data fields. */

    for (i = 0; i < nfields; i++) {
      dzero (nptin,  datain,  1);
      dzero (nptout, dataout, 1);

      for (j = 0; j < nzin; j++) {
	if (fread (datain+j*ntot, sizeof (double), nplane, fp_in) != nplane)
	  message (prog, "an error occured while reading", ERROR);
	if (swab) dbrev (ntot, datain+j*ntot, 1, datain+j*ntot, 1);
      }

      if (force) { /* -- We can just copy in physical space. */
	for (k = 0; k < nrep; k++) { 
	  for (j = 0; j < nzin; j++) {
	    if (fwrite (datain+j*ntot, sizeof (double), nplane, fp_out) != nplane)
	      message (prog, "an error occured while writing", ERROR);
	  }
	}
      } else {	   /* -- Have to go to Fourier space for padding. */
	dDFTr (datain,  nzin,  ntot, FORWARD);
	pack  (datain,  nzin,  dataout, nzout, nrep, ntot);
	dDFTr (dataout, nzout, ntot, INVERSE);

	for (j = 0; j < nzout; j++)
	  if (fwrite (dataout+j*ntot, sizeof (double), nplane, fp_out) != nplane)
	    message (prog, "an error occured while writing", ERROR);
      }
    }
  } 

  freeDvector (datain,  0);
  freeDvector (dataout, 0);
  
  return EXIT_SUCCESS;
}


static void getargs (int    argc ,
		     char** argv ,
		     FILE** fp_in,
		     int*   nrep ,
		     int*   force)
/* ------------------------------------------------------------------------- *
 * Parse command line arguments.
 * ------------------------------------------------------------------------- */
{
  char c, fname[FILENAME_MAX];
  char usage[] = "repeatz [-h] [-n <rep>] [-f] [input[.fld]\n";

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fputs (usage, stderr);
      exit  (EXIT_SUCCESS);
      break;
    case 'f':
      *force = 1;
      break;
    case 'n':
      if (*++argv[0]) *nrep = atoi (*argv);
      else           {*nrep = atoi (*++argv); argc--;}
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
}


static int roundup (const int nzdes)
/* ------------------------------------------------------------------------- *
 * Roundup number of z planes to suit Fourier transformation.
 * ------------------------------------------------------------------------- */
{
  int n, nz = nzdes, ip, iq, ir, ipqr2;

  if (nz != 1) {
    do {
      n = nz;
      prf235 (&n, &ip, &iq, &ir, &ipqr2);
      nz += (n) ? 0 : 2;
    } while (n == 0);
  }

  return nz;
}


static void pack (const double* datain ,  
		  const int     nzin   ,  
		  double*       dataout, 
		  const int     nzout  , 
		  const int     nrep   ,
		  const int     psize  )
/* ------------------------------------------------------------------------- *
 * Place planes of datain in the appropriate planes of dataout
 * (Fourier space).
 * ------------------------------------------------------------------------- */
{
  const int nmode = nzin  >> 1;
  const int msize = psize << 1;
  int       i;

  dcopy (psize, datain, 1, dataout, 1);
  if (nzin > 1) dcopy (psize, datain + psize, 1, dataout + psize, 1);

  for (i = 1; i < nmode; i++)
    dcopy (msize, datain + i * msize, 1, dataout + i * nrep * msize, 1);
}
