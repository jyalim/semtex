/******************************************************************************
 * slit.c: reproduce specified columns of input on output.
 *
 * Copyright (c) 1990 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn
 *
 * Usage: slit [-c <colstr>] [file], where <colstr> is a
 * comma-separated list of column numbers.
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

static char RCS[] = "$Id: slit.c,v 8.1 2015/04/20 11:14:19 hmb Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define MAXCOL	256
#define	MAXSTR	256
#define NUMDIG	8

int  parse    (char *strin, char *strout, int *pos, char sep);
void slit     (FILE *fp, int n, int *col);
int  get_line (FILE *fp, char coltext[][MAXSTR], int *nwords);


main(int argc, char *argv[])
/* ========================================================================= *
 * This does adminstration for the routines which do all the work.
 * ========================================================================= */
{
  int	i, ncol=1, colnum[MAXCOL];
  char	num[NUMDIG];
  FILE *fp;
  
  colnum[0] = 0;
  while (argc > 1 && argv[1][0] == '-') {
    if (argv[1][1] == 'c') {
      argc--; argv++;
      i = 0;  ncol = 0;
      while (parse(argv[1], num, &i, ',') != 0 && ncol < MAXCOL) {
	colnum[ncol] = atoi(num) - 1;
	if (colnum[ncol] < 0) {
	  (void)fprintf(stderr, "slit: bad column\n");
	  exit(1);
	}
	ncol++;
      }
    }
    else {
      (void)fprintf(stderr, "%s: unknown arg %s\n",
		    argv[0], argv[1]);
      exit(1);
    }
    argc--; argv++;
  }
  
  if (argc == 1) {	/* no file given, slit stdin */
    slit(stdin, ncol, colnum);
  }
  else {
    for (i = 1; i < argc; i++) {
      if ((fp=fopen(argv[i],"r")) == NULL) {
	(void)fprintf(stderr,"slit: can't open %s\n",argv[i]);
	(void)fprintf(stderr,"Usage: slit [-c <chanstr>]");
	(void)fprintf(stderr," [<filename>]\n");
	exit(1);
      }
      else {
	slit(fp, ncol, colnum);
	fclose(fp);
      }
    }
  }
  exit(0);
}


int parse(char *strin, char *strout, int *pos, char sep)
/* =========================================================================
 * Parse strin into strout until a non-solid character or separator
 * character occurs in strin or end of strin is reached.
 *
 * Pos is the index of the first character in strin to be examined.
 * Update pos to point at the next non-separator character in strin.
 *
 * Parse returns the number of characters parsed from strin.
 * ========================================================================= */
{
  int	i = 0;
  
  while (    (strin[*pos] != sep ) 
	 && !isspace(strin[*pos] ) 
	 &&  (strin[*pos] != '\0') ) {
    strout[i] = strin[*pos];
    (*pos)++;
    i++;
  }

  if (strin[*pos] == sep) (*pos)++;
  strout[i] = '\0';

  return i;
}


void slit(FILE *fp, int n, int *col)
/* =========================================================================
 * Slit into columns, print up.
 * ========================================================================= */
{
  int	i, nwords;
  char	coltext[MAXCOL][MAXSTR];
  
  while (get_line(fp, coltext, &nwords) != EOF) {
    if (nwords > 0) {
      (void)printf("%s", coltext[col[0]]);
      for (i=1; i<n; i++) 
	(void)printf(" %s", coltext[col[i]]);
    }
    (void)printf("\n");
  }
}


int get_line(FILE *fp, char coltext[][MAXSTR], int *nwords)
/* =========================================================================
 * The parsing of each input line is done here.
 * ========================================================================= */
{
  int	c, i, rec, inword;
  
  rec    = -1;
  inword =  0;
  *nwords = 0;

  while ((c=getc(fp)) != EOF && c != '\n') {
    if (!isspace(c)) {
      if (!inword) {
	inword = 1;
	(*nwords)++;
	rec++;
	i = 0;
	coltext[rec][i] = c;
	i++;
      } else {
	coltext[rec][i] = c;
	i++;
      }
    } else
      inword = 0;

    if (rec >= 0) coltext[rec][i] = '\0';
  }

  return c;
}

