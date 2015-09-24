/*****************************************************************************
 * TOPRE.C: make a NEKTON preprocessor input file. (See NEK5000 site.)
 *
 * Copyright (c) 1994 <--> $Date: 2015/04/20 11:14:19 $, Hugh Blackburn
 *
 * Generate, from an indexed list of 2D points and a list of index
 * quadruplets, the element-description part of a text file to be used as
 * input for the NEKTON preprocessor, PRE.  The quadruplets give the indices
 * of quadrilateral element corner vertices in CCW traverse.
 *
 * Example input file for two elements follows:
 * # Any number of comment lines beginning with "#";
 * 6 VERTICES
 * 1 c  0.0   0.0
 * 2 c  1.0   0.0
 * 3 c  1.0   1.0
 * 4 c  0.0   1.0
 * 5 c  1.0   2.0
 * 6 p  2.0  90.0
 * [ blank line]
 * 2 ELEMENTS
 * 1  1 2 3 4
 * 2  4 3 5 6
 * [ EOF ]
 *
 * Vertex locations can be specified in Cartesian ("c") notation, or in
 * polar ("p") notation (radius, theta; theta in degrees).  Points are
 * stored internally in Cartesian form.
 *
 * No internal checking of consistency is done, but, as a part of every run,
 * an ASCII file called "sm.dat" is generated, so that sm can be used for a
 * quick visual check of the mesh (using Ron's "meshplot" macro set).
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

static char RCS[] = "$Id: topre.c,v 8.1 2015/04/20 11:14:19 hmb Exp $";



#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include <cveclib>

static char prog[] = "topre";

#define  readLine(fp)  fgets(buf, STR_MAX, (fp))
#define  strip(s)      (s)[strlen((s))-1] = '\0'
#define  skipComments  if (*buf == '#') continue; else if (*buf == '\n') break
#define  upperCase(s)  { char *z=(s); while (*z=toupper(*z)) z++; }

#define  d2r      0.01745329251994329576
#define  r2d     57.29577951308232087721
#define  rad(x)  (x * d2r)
#define  deg(x)  (x * r2d)

static   zomplex  P2C    (zomplex);
static   zomplex  C2P    (zomplex);
static   void     header (void);





int main()
/* ========================================================================= *
 * Driver.                                                                   *
 * ========================================================================= */
{
  int        i, n, nvert, nel;
  zomplex   *vertex, Z;
  int      **elmt;
  FILE      *sm;


  while (readLine(stdin)) if (*buf != '#') break;

  /* INPUT VERTEX INFORMATION */

  if (strstr(buf, "VERTICES"))
    sscanf(buf, "%d", &nvert);
  else
    message("topre: expected number of vertices, got: ", buf, ERROR);

  vertex = zvector(1, nvert);

  n = 0;
  while (readLine(stdin)) {
    skipComments;
    sscanf(buf, "%d %*s %lf %lf", &i, &Z.Re, &Z.Im);
    if (i > nvert)
      message(prog, "vertex number exceeds previous declaration", ERROR);
    if (strstr(buf, "p")) {
      Z.Im = rad(Z.Im);
      Z    = P2C(Z);
    } else if (!strstr(buf, "c"))
      message("topre: bad vertex specification", buf, ERROR);
    vertex[i] = Z;
    n++;
  }
  if (n != nvert) {
    sprintf (buf, "declared %1d vertices, read %1d", nvert, n);
    message(prog, buf, ERROR);
  }


  /* INPUT ELEMENT INFORMATION */
    
  readLine(stdin); strip(buf); upperCase(buf);
  if (strstr(buf, "ELEMENTS"))
    sscanf(buf, "%d", &nel);
  else
    message("topre: expected number of elements, got: ", buf, ERROR);

  elmt = imatrix(1, nel, 0, 3);
  n = 0;
  while (readLine(stdin) && n < nel) {
    skipComments;
    n++;
    sscanf(buf, "%d", &i);
    if (i > nel)
      message("topre: an element index exceeds delared max: ", buf, ERROR);
    sscanf(buf, "%*s %d %d %d %d", elmt[i], elmt[i]+1, elmt[i]+2, elmt[i]+3);
  }

  if (n != nel)
    message(prog, "number of elements declared exceeds number read", ERROR);

  /* GENERATE SM DATA FILE */

  sm = fopen("sm.dat", "w");
  fprintf(sm , "2 2 1 %1d NR NS NZ NEL\n", nel);
  for (i=1; i<=nel; i++) {
    fprintf(sm, "%15.6f%15.6f\n", vertex[elmt[i][0]].Re,vertex[elmt[i][0]].Im);
    fprintf(sm, "%15.6f%15.6f\n", vertex[elmt[i][1]].Re,vertex[elmt[i][1]].Im);
    fprintf(sm, "%15.6f%15.6f\n", vertex[elmt[i][3]].Re,vertex[elmt[i][3]].Im);
    fprintf(sm, "%15.6f%15.6f\n", vertex[elmt[i][2]].Re,vertex[elmt[i][2]].Im);
  }
  fclose(sm);

  /* GENERATE NEKTON FORMAT OUTPUT */

  header();

  printf(" **MESH DATA** 1st line is X of corner 1,2,3,4. 2nd line is Y.\n");
  printf("%6d%6d%6d           NEL,NDIM,NELV\n", nel, 2, nel);
  for (i=1; i<=nel; i++) {
    printf("            ELEMENT%5d   [  1A]    GROUP     0\n", i);
    printf("%14.6f%14.6f%14.6f%14.6f\n",
	   vertex[elmt[i][0]].Re, vertex[elmt[i][1]].Re,
	   vertex[elmt[i][2]].Re, vertex[elmt[i][3]].Re);  
    printf("%14.6f%14.6f%14.6f%14.6f\n",
	   vertex[elmt[i][0]].Im, vertex[elmt[i][1]].Im,
	   vertex[elmt[i][2]].Im, vertex[elmt[i][3]].Im);
  } 

  return 0;
}





static zomplex P2C(zomplex z)
/* ========================================================================= *
 * Polar-->Cartesian conversion.                                             *
 * ========================================================================= */
{
  zomplex c;

  c.Re = z.Re * cos(z.Im);
  c.Im = z.Re * sin(z.Im);

  return c;
}





static zomplex C2P(zomplex z)
/* ========================================================================= *
 * Cartesian-->polar conversion.                                             *
 * ========================================================================= */
{
  zomplex c;

  c.Re = sqrt(z.Re*z.Re + z.Im*z.Im);
  c.Im = atan2(z.Im, z.Re);

  return c;
}





static void header(void)
/* ========================================================================= *
 * Print PRE style header.                                                   *
 * ========================================================================= */
{
  char *hdr[] = {
    " ****** PARAMETERS *****\n",
    "    2.500000     NEKTON VERSION \n",
    "            2 DIMENSIONAL RUN\n",
    "            9 PARAMETERS FOLLOW\n",
    "   1.00000         DENSITY   \n",
    "   1.00000         KINVIS    \n",
    "   1.00000         NSTEPS    \n",
    "  0.010000E+00     DT        \n",
    "  0.000000E+00     EQTYPE    \n",
    "  -1.00000         INTYPE    \n",
    "   5.00000         NORDER    \n",
    "  0.100000E-01     TOLREL    \n",
    "   1.00000         TOLABS    \n",
    "      4  Lines of passive scalar data follows2 CONDUCT; 2RHOCP\n",
    "   1.00000       1.00000       1.00000       1.00000       1.00000    \n",
    "   1.00000       1.00000       1.00000       1.00000    \n",
    "   1.00000       1.00000       1.00000       1.00000       1.00000    \n",
    "   1.00000       1.00000       1.00000       1.00000    \n",
    "           11  LOGICAL SWITCHES FOLLOW\n",
    "  T     IFFLOW\n",
    "  F     IFHEAT\n",
    "  T     IFTRAN\n",
    "  F F F F F F F F F F F IFNAV & IFADVC (Advection in P.S. fields)\n",
    "  F F T T T T T T T T T IFTMSH (IF mesh for this field is T mesh)\n",
    "  F     IFAXIS\n",
    "  F     IFSTRS\n",
    "  F     IFSPLIT\n",
    "  F     IFMGRID\n",
    "  F     IFMODEL\n",
    "  F     IFKEPS\n",
    "    80.000        80.000      -15.0000      -40.0000"
      "     XFAC,YFAC,XZERO,YZERO\n"
  };
  int i;

  for (i=0; i<31; i++) printf(hdr[i]);
}
