/*****************************************************************************
 * TESTRAND.C: exercise random number generation routines in
 * veclib. The main idea is justto runthem and assess their pdfs
 * independently.
 *****************************************************************************/

static char
RCSid[] = "$Id";

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include <veclib.h>
#include <femlib.h>

#define  VEC_MAX 7500

int main (int argc, char** argv)
/* ------------------------------------------------------------------------- *
 * 
 * ------------------------------------------------------------------------- */
{
  int_t  i, seed;
  real_t x;
  real_t dum[VEC_MAX], vectr[VEC_MAX];

  seed = -(char) time(NULL);

  raninit (seed);

#if 0
  // -- Straight call to underlying C routine.

  for (i = 0; i < VEC_MAX; i++) {
    x = dnormal (0.0, 0.01);
  }

#else
  // -- Call through the parser.

  for (i = 0; i < VEC_MAX; i++) dum[i] = i;

  Femlib::initialize (&argc, &argv);
  Femlib::prepVec  ("dum", "white(0.01)");
  Femlib__parseVec (VEC_MAX, &dum, &vectr);

#endif

  for (i = 0; i < VEC_MAX; i++) {
    printf ("%g\n", vectr[i]);
  }
  return EXIT_SUCCESS;
}

