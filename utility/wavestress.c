/*****************************************************************************
 * wavestress.c: from a 2D3C but complex modal data file, compute 2D3C
 * distributions of streamwise-averaged Reynolds stresses.  Input data
 * must be binary format and contain only fields u v w p.
 *
 * NB this version does not symmetrise the Fourier-direction stress.
 *
 * Copyright (c) 2011 <--> $Date: 2015/08/13 06:16:20 $, Hugh Blackburn
 *
 * USAGE
 * -----
 * wavestress [-h] [input[.fld]]
 *
 * INPUT FILE
 * ----------
 * Contains only fields uvwp and has N_Z = 2 (Real and Imaginary parts).
 *
 * OUTPUT FILE
 * -----------
 * Is a standard 2D/real (N_Z = 1) Reynolds stress file containing
 * uvwpABCDEF, with
 *
 * u = 0.5 * sqrt (u.Re^2 + u.Im^2 + v.Re^2 + v.Im^2)
 * v = 0.5 * sqrt (u.Re^2 + u.Im^2 + v.Re^2 + v.Im^2 + w.Re^2 + w.Im^2)
 * w = 0 FIXME (but it's OK because this area is not used by anything else)
 * p = sqrt (p.Re^2 + p.Im^2)
 * A = 2*(u.Re^2    + u.Im^2)
 * B = 2*(u.Re*v.Re + u.Im*v.Im)
 * C = 2*(v.Re^2    + v.Im^2)
 * D = 2*(u.Re*w.Re + u.Im*w.Im)
 * E = 2*(v.Re*w.Re + v.Im*w.Im)
 * F = 2*(w.Re^2    + w.Im^2)
 *
 *****************************************************************************/

static char RCS[] = "$Id: wavestress.c,v 8.1 2015/08/13 06:16:20 hmb Exp $";

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include <cfemlib.h>
#include <cfemdef.h>
#include <cveclib.h>

static void getargs (int, char**, FILE**);
static char prog[] = "wavestress";
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
  char   buf[STR_MAX], fields[STR_MAX], fmt[STR_MAX];
  int    i, j, n, np, nz, nel, swab = 0;
  int    nfields, nplane;
  FILE   *fp_in = stdin, *fp_out = stdout;
  double **idata, **odata, *vcmpt1, *vcmpt2;

  getargs (argc, argv, &fp_in);
  format  (fmt);

  while (fgets (buf, STR_MAX, fp_in)) {

    /* -- Process header. */

    fputs (buf, fp_out); fgets (buf, STR_MAX, fp_in);
    fputs (buf, fp_out); fgets (buf, STR_MAX, fp_in);

    if (sscanf (buf, "%d%*s%d%d", &np, &nz, &nel) != 3)
      message (prog, "unable to read the file size", ERROR);

    if (nz != 2) {
      sprintf (fields, "input must have nz = 2 (here %1d)", nz);
      message (prog, fields, ERROR);
    }
    fprintf (fp_out, hdr_fmt[2], np, np, 1, nel);

    n = 6;
    while (--n) {
      fgets (buf, STR_MAX, fp_in);
      fputs (buf, fp_out);
    }

    fgets(fields, STR_MAX, fp_in);
    memset(fields+25, '\0', STR_MAX-25);
    for (nfields = 0, i = 0; i < 25; i++) if (isalpha(fields[i])) nfields++;
    if (!((nfields == 4) && (strstr (fields, "uvwp"))))
	message (prog, "input must have only fields u v w p.", ERROR);
    fprintf (fp_out, hdr_fmt[8], "uvwpABCDEF");

    fgets (buf, STR_MAX, fp_in);
    for (i = 0; i < strlen (buf); i++) buf[i] = tolower (buf[i]);

    if (!strstr(buf, "binary"))
      message (prog, "input file not binary format", ERROR);
    if (!strstr (buf, "endian"))
      message (prog, "input field file in unknown binary format", WARNING);
    else
      swab = ((strstr (buf, "big") && strstr (fmt, "little")) ||
	      (strstr (fmt, "big") && strstr (buf, "little")) );
    sprintf (buf, "%s %s", "binary", fmt);
    fprintf (fp_out, hdr_fmt[9], buf);

    /* -- Set sizes, allocate storage, set to zero. */

    nplane = np * np * nel;

    idata = dmatrix (0, 3, 0, nplane * 2);
    odata = dmatrix (0, 9, 0, nplane); /* -- uvwpABCDEF = 10 */

    dzero (4*nplane,  idata[0], 1);
    dzero (10*nplane, odata[0], 1);

    /* -- Read in all data fields. */

    for (i = 0; i < nfields; i++) {
      if (fread (idata[i], sizeof (double), nplane * 2, fp_in) != nplane * 2)
	message (prog, "an error occured while reading", ERROR);
      if (swab) dbrev (nplane*2, idata[i], 1, idata[i], 1);
    }

    /* -- Compute A. */

    vcmpt1 = idata[0];		/* -- Real part of u. */
    vcmpt2 = idata[0] + nplane; /* -- Imag part of u. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt1, 1, odata[4], 1, odata[4], 1);
    dvvtvp (nplane, vcmpt2, 1, vcmpt2, 1, odata[4], 1, odata[4], 1);
    dscal  (nplane, 2.0, odata[4], 1);

    /* -- Compute B . */

    vcmpt1 = idata[0];		/* -- Real part of u. */
    vcmpt2 = idata[1];		/* -- Real part of v. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt2, 1, odata[5], 1, odata[5], 1);

    vcmpt1 = idata[0] + nplane; /* -- Imag part of u. */
    vcmpt2 = idata[1] + nplane;	/* -- Imag part of v. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt2, 1, odata[5], 1, odata[5], 1);

    dscal  (nplane, 2.0, odata[5], 1);

    /* -- Compute C. */

    vcmpt1 = idata[1];		/* -- Real part of v. */
    vcmpt2 = idata[1] + nplane; /* -- Imag part of v. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt1, 1, odata[6], 1, odata[6], 1);
    dvvtvp (nplane, vcmpt2, 1, vcmpt2, 1, odata[6], 1, odata[6], 1);
    dscal  (nplane, 2.0, odata[6], 1);

    /* -- Compute D . */

    vcmpt1 = idata[0];          /* -- Real part of u. */
    vcmpt2 = idata[2];          /* -- Real part of w. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt2, 1, odata[7], 1, odata[7], 1);

    vcmpt1 = idata[0] + nplane; /* -- Imag part of u. */
    vcmpt2 = idata[2] + nplane; /* -- Imag part of w. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt2, 1, odata[7], 1, odata[7], 1);

    dscal  (nplane, 2.0, odata[7], 1);

    /* -- Compute E . */

    vcmpt1 = idata[1];          /* -- Real part of v. */
    vcmpt2 = idata[2];          /* -- Real part of w. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt2, 1, odata[8], 1, odata[8], 1);

    vcmpt1 = idata[1] + nplane; /* -- Imag part of v. */
    vcmpt2 = idata[2] + nplane; /* -- Imag part of w. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt2, 1, odata[8], 1, odata[8], 1);

    dscal  (nplane, 2.0, odata[8], 1);

    /* -- Compute F. */

    vcmpt1 = idata[2];          /* -- Real part of w. */
    vcmpt2 = idata[2] + nplane; /* -- Imag part of w. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt1, 1, odata[9], 1, odata[9], 1);
    dvvtvp (nplane, vcmpt2, 1, vcmpt2, 1, odata[9], 1, odata[9], 1);
    dscal  (nplane, 2.0, odata[9], 1);

    /* -- Compute p. */

    vcmpt1 = idata[3];		/* -- Real part of p. */
    vcmpt2 = idata[3] + nplane; /* -- Imag part of p. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt1, 1, odata[2], 1, odata[2], 1);
    dvvtvp (nplane, vcmpt2, 1, vcmpt2, 1, odata[2], 1, odata[2], 1);
    dvsqrt (nplane, odata[2], 1, odata[2], 1);

    /* -- Compute u & v. */

    vcmpt1 = idata[0];		/* -- Real part of u. */
    vcmpt2 = idata[0] + nplane; /* -- Imag part of u. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt1, 1, odata[0], 1, odata[0], 1);
    dvvtvp (nplane, vcmpt2, 1, vcmpt2, 1, odata[0], 1, odata[0], 1);
    dvvtvp (nplane, vcmpt1, 1, vcmpt1, 1, odata[1], 1, odata[1], 1);
    dvvtvp (nplane, vcmpt2, 1, vcmpt2, 1, odata[1], 1, odata[1], 1);
    vcmpt1 = idata[1];		/* -- Real part of v. */
    vcmpt2 = idata[1] + nplane; /* -- Imag part of v. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt1, 1, odata[0], 1, odata[0], 1);
    dvvtvp (nplane, vcmpt2, 1, vcmpt2, 1, odata[0], 1, odata[0], 1);
    dvvtvp (nplane, vcmpt1, 1, vcmpt1, 1, odata[1], 1, odata[1], 1);
    dvvtvp (nplane, vcmpt2, 1, vcmpt2, 1, odata[1], 1, odata[1], 1);
    vcmpt1 = idata[1];		/* -- Real part of w. */
    vcmpt2 = idata[1] + nplane; /* -- Imag part of w. */
    dvvtvp (nplane, vcmpt1, 1, vcmpt1, 1, odata[1], 1, odata[1], 1);
    dvvtvp (nplane, vcmpt2, 1, vcmpt2, 1, odata[1], 1, odata[1], 1);

    dvsqrt (nplane, odata[0], 1, odata[0], 1);
    dsmul  (nplane, 0.5, odata[0], 1, odata[0], 1);
    dvsqrt (nplane, odata[1], 1, odata[1], 1);
    dsmul  (nplane, 0.5, odata[1], 1, odata[1], 1);

    /* FIXME: compute w if needed in future. */

    /* -- Write out uvwpABCDEF in binary. */

    for (i = 0; i < 10; i++)
      if (fwrite (odata[i], sizeof (double), nplane, fp_out) != nplane)
	message (prog, "an error occured while writing", ERROR);

    freeDmatrix (idata, 0, 0);
    freeDmatrix (odata, 0, 0);
  }

  return EXIT_SUCCESS;
}


static void getargs (int    argc ,
		     char** argv ,
		     FILE** fp_in)
/* ------------------------------------------------------------------------- *
 * Parse command line arguments.
 * ------------------------------------------------------------------------- */
{
  char c, fname[FILENAME_MAX];
  char usage[] = "wavestress [-h] [input[.fld]]\n";

  while (--argc && (*++argv)[0] == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fputs (usage, stderr);
      exit  (EXIT_SUCCESS);
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
