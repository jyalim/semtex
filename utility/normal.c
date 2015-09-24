/*****************************************************************************
 * normal.c: use veclib's random number generation routines (ran2 +
 * gasdev) to produce Gaussian-distributed pseudo-random numbers with
 * nominated mean and standard deviation (defaults: 0 & 1).
 *
 * Copyright (c) 2008 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn
 *
 * USAGE
 * -----
 * normal [-h] [-n num] [-m mean] [-d sdev] [-s seed]
 *
 * A positive value of seed causes the random number seed to be
 * generated from wall-clock time (and the value is otherwise
 * irrelevant). A negative value is used directly for the seed.
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

static char RCS[] = "$Id: normal.c,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>

#include <cfemlib.h>
#include <cfemdef.h>
#include <cveclib.h>

static void getargs (int, char**, int_t*, int_t*, real_t*, real_t*);
static char prog[] = "normal";


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Wrapper.
 * ------------------------------------------------------------------------- */
{
  int_t  i, n = 1, seed = 0;
  real_t mean = 0.0, sdev = 1.0;

  getargs (argc, argv, &n, &seed, &mean, &sdev);

#if 1
  raninit (seed);
#else
  raninit (-((short)time(NULL)));
#endif
  for (i = 0; i < n; i++) printf ("%g\n", dnormal(mean, sdev));

  return EXIT_SUCCESS;
}


static void getargs (int     argc,
		     char**  argv,
		     int_t*  n   ,
		     int_t*  seed,
		     real_t* mean,
		     real_t* sdev)
/* ------------------------------------------------------------------------- *
 * Parse command line arguments.
 * ------------------------------------------------------------------------- */
{
  char c;
  int  i;
  char fname[FILENAME_MAX];
  char usage[] = "usage: random [options]\n"
    "options:\n"
    "-h         ... print this help message\n"
    "-n num     ... generate num random numbers [Default: 1]\n"
    "-s seed    ... set random number seed      [Default: 0]\n"
    "-m mean    ... set mean value              [Default: 0.0]\n"
    "-d sdev    ... set standard deviation      [Default: 1.0]\n";


  while (--argc && (*++argv)[0] == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fputs (usage, stderr);
      exit  (EXIT_SUCCESS);
      break;
    case 'n':
      if (*++argv[0]) *n = atoi (*argv);
      else {*n = atoi (*++argv); argc--;}
      break;
    case 's':
      if (*++argv[0]) *seed = atoi (*argv);
      else {*seed = atoi (*++argv); argc--;}
      break;
    case 'm':
      if (*++argv[0]) *mean = atof (*argv);
      else {*mean = atof (*++argv); argc--;}
      break;
    case 'd':
      if (*++argv[0]) *sdev = atof (*argv);
      else {*sdev = atof (*++argv); argc--;}
      break;
    default:
      fprintf (stderr, "%s: unknown option -- %c\n", prog, c);
      break;
    }

  return;
}
