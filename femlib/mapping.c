/*****************************************************************************
 * mapping.c: build/return integer mapping vectors used to gather/scatter
 * element storage formats from one to another.
 *
 * Copyright (c) 1994 <--> $Date: 2015/04/20 11:14:14 $, Hugh Blackburn
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
 *
 * $Id: mapping.c,v 8.1 2015/04/20 11:14:14 hmb Exp $
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include <cfemdef.h>
#include <cveclib.h>
#include <cfemlib.h>

typedef struct mapping {
  int_t         np   ;
  int_t         dim  ;
  int_t*        emap ;
  int_t*        pmap ;
  struct mapping* next ;
} Mapping;

static Mapping* mHead = 0;


void edgemaps (const int_t nk ,
	       const int_t dim,
	       int_t**     map,
	       int_t**     inv)
/* ------------------------------------------------------------------------- *
 * An (e)map is an edge-major permutation of indices, going CCW around the
 * element edges and then traversing the interior in row-major form. 
 * This allows access element storage by traversing edges, tying in with
 * edge-based global numbering.
 *
 * Pmap (the inversion of emap, or inv above) is built afterwards.
 *
 * To convert from row-major storage to boundary-major, then back:
 *   gathr (ntot, value, emap, tmp  );
 *   gathr (ntot, tmp,   pmap, value)  OR  scatr (ntot, tmp, emap, value);
 *
 * If the required maps are "in stock", they are returned from a list,
 * otherwise they are created and added to the list as a first step.
 * NULL input pointers are not assigned a value.
 * ------------------------------------------------------------------------- */
{
  char              routine[] = "edgemaps";
  const int_t     len = (dim == 2) ? nk * nk : nk;
  register int_t  found = 0;
  register Mapping* p;

  if (nk < 2)
    message (routine, "input nk < 2", ERROR);

  if (dim < 1 || dim > 2)
    message (routine, "number of space dimensions must be 1 or 2", ERROR);

  for (p = mHead; p; p = p -> next)
    if (found = nk == p -> np && dim == p -> dim) break;

  if (!found) {
    register int_t i, j, k, n;
    const    int_t nm = nk - 1;
    register int_t *em, *pm;

    p = (Mapping *) calloc (1, sizeof (Mapping));
    if (mHead) p -> next = mHead;
    mHead = p;

    p -> np   = nk;
    p -> dim  = dim;
    p -> emap = ivector (0, len);
    p -> pmap = ivector (0, len);

    em = p -> emap;
    pm = p -> pmap;

    if (dim == 1) {
      
      em[0] = 0;
      em[1] = nk - 1;
      for (i = 2; i < nk; i++) em[i] = i - 1;

    } else if (dim == 2) {

      /* -- Traverse exterior CCW. */

      k = 0;
      n = 0;
      em[0] = 0;
      for (i = 1; i < nk; i++) em[++k] = n += 1;
      for (i = 1; i < nk; i++) em[++k] = n += nk;
      for (i = 1; i < nk; i++) em[++k] = n -= 1;
      for (i = 1; i < nk; i++) em[++k] = n -= nk;
  
      /* -- Traverse interior in row-major order. */

      for (i = 1; i < nm; i++)
	for (j = 1; j < nm; j++)
	  em[k++] = i * nk + j;
    }

    /* -- Build pmap. */
  
    for (i = 0; i < len; i++) pm[em[i]] = i;
  }

  /* -- p now points to valid storage: return requested operators. */

  if (map) *map = p -> emap;
  if (inv) *inv = p -> pmap;
}
