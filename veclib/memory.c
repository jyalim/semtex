/*****************************************************************************
 *                      MEMORY ALLOCATION UTILITIES
 *
 * $Id: memory.c,v 8.1 2015/04/20 11:14:19 hmb Exp $
 *****************************************************************************/

#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>

#include <cfemdef.h>
#include <cveclib.h>

#define FREE_ARG  void*

  
complex* cvector (int_t nl, int_t nh)
/* ------------------------------------------------------------------------- *
 * Allocates a complex vector with subscript range [nl..nh].
 * ------------------------------------------------------------------------- */
{
  complex* v;

  v = (complex*) malloc((size_t) ((nh-nl+1)*sizeof(complex)));
  if (v) return v-nl;

  message("cvector()", "allocation failure", WARNING);
  return NULL;
}


double* dvector (int_t nl, int_t nh)
/* ------------------------------------------------------------------------- *
 * Allocates a double vector with subscript range [nl..nh].
 * ------------------------------------------------------------------------- */
{
  double* v;

  v = (double*) malloc((size_t) ((nh-nl+1)*sizeof(double)));
  if (v) return v-nl;

  message("dvector()", "allocation failure", WARNING);
  return NULL;
}


float* svector (int_t nl, int_t nh)
/* ------------------------------------------------------------------------- *
 * Allocates a float vector with range [nl..nh].
 * ------------------------------------------------------------------------- */
{
  float* v;

  v = (float*) malloc((size_t) ((nh-nl+1)*sizeof(float)));
  if (v) return v-nl;

  message("svector()", "allocation failure", WARNING);
  return NULL;
}


int_t* ivector (int_t nl, int_t nh)
/* ------------------------------------------------------------------------- *
 * Allocates an int vector with subscript range [nl..nh].
 * ------------------------------------------------------------------------- */
{
  int_t* v;
  
  v = (int_t*) malloc((size_t) ((nh-nl+1)*sizeof(int_t)));
  if (v) return v-nl;

  message("ivector()", "allocation failure", WARNING);
  return NULL;
}


zomplex* zvector(int_t nl, int_t nh)
/* ------------------------------------------------------------------------- *
 * Allocates a zomplex vector with subscript range [nl..nh].
 * ------------------------------------------------------------------------- */
{
  zomplex* v;

  v = (zomplex*) malloc((size_t) ((nh-nl+1)*sizeof(zomplex)));
  if (v) return v-nl;

  message("zvector()", "allocation failure", WARNING);
  return NULL;
}


complex **cmatrix (int_t nrl, int_t nrh,
		   int_t ncl, int_t nch)
/* ------------------------------------------------------------------------- *
 * Allocate a complex matrix with subscript ranges [nrl..nrh][ncl..nch].
 * ------------------------------------------------------------------------- */
{
  int_t i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
  complex **m;

  m = (complex**) malloc((size_t) (nrow*sizeof(complex*)));
  if (!m) {
    message("cmatrix()", "allocation failure 1", WARNING);
    return NULL;
  }
  m -= nrl;

  m[nrl] = (complex*) malloc((size_t) (nrow*ncol*sizeof(complex)));
  if (!m[nrl]) {
    message("cmatrix()", "allocation failure 2", WARNING);
    return NULL;
  }
  m[nrl] -= ncl;

  for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;

  return m;
}


double **dmatrix (int_t nrl, int_t nrh,
		  int_t ncl, int_t nch)
/* ------------------------------------------------------------------------- *
 * Allocate a double  matrix with subscript ranges [nrl..nrh][ncl..nch].
 * ------------------------------------------------------------------------- */
{
  int_t i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
  double  **m;

  m = (double**) malloc((size_t) (nrow*sizeof(double*)));
  if (!m) {
    message("dmatrix()", "allocation failure 1", WARNING);
    return NULL;
  }
  m -= nrl;

  m[nrl] = (double*) malloc((size_t) (nrow*ncol*sizeof(double)));
  if (!m[nrl]) {
    message("dmatrix()", "allocation failure 2", WARNING);
    return NULL;
  }
  m[nrl] -= ncl;

  for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;

  return m;
}


float **smatrix (int_t nrl, int_t nrh,
		 int_t ncl, int_t nch)
/* ------------------------------------------------------------------------- *
 * Allocate a float matrix with subscript ranges [nrl..nrh][ncl..nch].
 * ------------------------------------------------------------------------- */
{
  int_t i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
  float   **m;

  m = (float**) malloc((size_t) (nrow*sizeof(float*)));
  if (!m) {
    message("smatrix()", "allocation failure 1", WARNING);
    return NULL;
  }
  m -= nrl;

  m[nrl] = (float*) malloc((size_t) (nrow*ncol*sizeof(float)));
  if (!m[nrl]) {
    message("smatrix()", "allocation failure 2", WARNING);
    return NULL;
  }
  m[nrl] -= ncl;

  for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;

  return m;
}


int_t **imatrix (int_t nrl, int_t nrh,
		   int_t ncl, int_t nch)
/* ------------------------------------------------------------------------- *
 * Allocate an int matrix with subscript ranges [nrl..nrh][ncl..nch].
 * ------------------------------------------------------------------------- */
{
  int_t i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
  int_t **m;

  m = (int_t**) malloc((size_t) (nrow*sizeof(int_t*)));
  if (!m) {
    message("imatrix()", "allocation failure 1", WARNING);
    return NULL;
  }
  m -= nrl;

  m[nrl] = (int_t*) malloc((size_t) (nrow*ncol*sizeof(int_t)));
  if (!m[nrl]) {
    message("imatrix()", "allocation failure 2", WARNING);
    return NULL;
  }
  m[nrl] -= ncl;

  for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;

  return m;
}


zomplex **zmatrix (int_t nrl, int_t nrh,
		   int_t ncl, int_t nch)
/* ------------------------------------------------------------------------- *
 * Allocate a zomplex matrix with subscript ranges [nrl..nrh][ncl..nch].
 * ------------------------------------------------------------------------- */
{
  int_t i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
  zomplex **m;

  m = (zomplex**) malloc((size_t) (nrow*sizeof(zomplex*)));
  if (!m) {
    message("zmatrix()", "allocation failure 1", WARNING);
    return NULL;
  }
  m -= nrl;

  m[nrl] = (zomplex*) malloc((size_t) (nrow*ncol*sizeof(zomplex)));
  if (!m[nrl]) {
    message("zmatrix()", "allocation failure 2", WARNING);
    return NULL;
  }
  m[nrl] -= ncl;

  for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;

  return m;
}


complex ***c3matrix (int_t nrl, int_t nrh,
		     int_t ncl, int_t nch,
		     int_t ndl, int_t ndh)
/* ------------------------------------------------------------------------- *
 * Allocate a complex 3-matrix with ranges [nrl..nrh][ncl..nch][ndl..ndh].
 * ------------------------------------------------------------------------- */
{
  int_t i, j, nrow = nrh-nrl+1, ncol = nch-ncl+1, ndep = ndh-ndl+1;
  complex ***t;

  t = (complex***) malloc((size_t) (nrow*sizeof(complex**)));
  if (!t) {
    message("c3matrix()", "allocation failure 1", WARNING);
    return NULL;
  }
  t -= nrl;

  t[nrl] = (complex**) malloc((size_t) (nrow*ncol*sizeof(complex*)));
  if (!t[nrl]) {
    message("c3matrix()", "allocation failure 2", WARNING);
    return NULL;
  }
  t[nrl] -= ncl;

  t[nrl][ncl] = (complex*) malloc((size_t) (nrow*ncol*ndep*sizeof(complex)));
  if (!t[nrl][ncl]) {
    message("c3matrix()", "allocation failure 3", WARNING);
    return NULL;
  }
  t[nrl][ncl] -= ndl;

  for (j=ncl+1; j<=nch; j++)
    t[nrl][j] = t[nrl][j-1] + ndep;
  for (i=nrl+1; i<=nrh; i++) {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol*ndep;
    for (j=ncl+1; j<=nch; j++)
      t[i][j] = t[i][j-1] + ndep;
  }

  return t;
}

double ***d3matrix (int_t nrl, int_t nrh,
		    int_t ncl, int_t nch,
		    int_t ndl, int_t ndh)
/* ------------------------------------------------------------------------- *
 * Allocate a double 3-matrix with ranges [nrl..nrh][ncl..nch][ndl..ndh].
 * ------------------------------------------------------------------------- */
{
  int_t i, j, nrow = nrh-nrl+1, ncol = nch-ncl+1, ndep = ndh-ndl+1;
  double  ***t;

  t = (double***) malloc((size_t) (nrow*sizeof(double**)));
  if (!t) {
    message("d3matrix()", "allocation failure 1", WARNING);
    return NULL;
  }
  t -= nrl;

  t[nrl] = (double**) malloc((size_t) (nrow*ncol*sizeof(double*)));
  if (!t[nrl]) {
    message("d3matrix()", "allocation failure 2", WARNING);
    return NULL;
  }
  t[nrl] -= ncl;

  t[nrl][ncl] = (double*) malloc((size_t) (nrow*ncol*ndep*sizeof(double)));
  if (!t[nrl][ncl]) {
    message("d3matrix()", "allocation failure 3", WARNING);
    return NULL;
  }
  t[nrl][ncl] -= ndl;

  for (j=ncl+1; j<=nch; j++)
    t[nrl][j] = t[nrl][j-1] + ndep;
  for (i=nrl+1; i<=nrh; i++) {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol*ndep;
    for (j=ncl+1; j<=nch; j++)
      t[i][j] = t[i][j-1] + ndep;
  }

  return t;
}


float ***s3matrix (int_t nrl, int_t nrh,
		   int_t ncl, int_t nch,
		   int_t ndl, int_t ndh)
/* ------------------------------------------------------------------------- *
 * Allocate a float 3-matrix with ranges [nrl..nrh][ncl..nch][ndl..ndh].
 * ------------------------------------------------------------------------- */
{
  int_t i, j, nrow = nrh-nrl+1, ncol = nch-ncl+1, ndep = ndh-ndl+1;
  float   ***t;

  t = (float***) malloc((size_t) (nrow*sizeof(float**)));
  if (!t) {
    message("s3matrix()", "allocation failure 1", WARNING);
    return NULL;
  }
  t -= nrl;

  t[nrl] = (float**) malloc((size_t) (nrow*ncol*sizeof(float*)));
  if (!t[nrl]) {
    message("s3matrix()", "allocation failure 2", WARNING);
    return NULL;
  }
  t[nrl] -= ncl;

  t[nrl][ncl] = (float*) malloc((size_t) (nrow*ncol*ndep*sizeof(float)));
  if (!t[nrl][ncl]) {
    message("s3matrix()", "allocation failure 3", WARNING);
    return NULL;
  }
  t[nrl][ncl] -= ndl;

  for (j=ncl+1; j<=nch; j++)
    t[nrl][j] = t[nrl][j-1] + ndep;
  for (i=nrl+1; i<=nrh; i++) {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol*ndep;
    for (j=ncl+1; j<=nch; j++)
      t[i][j] = t[i][j-1] + ndep;
  }

  return t;
}


int_t ***i3matrix (int_t nrl, int_t nrh,
		     int_t ncl, int_t nch,
		     int_t ndl, int_t ndh)
/* ------------------------------------------------------------------------- *
 * Allocate an int 3-matrix with ranges [nrl..nrh][ncl..nch][ndl..ndh].
 * ------------------------------------------------------------------------- */
{
  int_t i, j, nrow = nrh-nrl+1, ncol = nch-ncl+1, ndep = ndh-ndl+1;
  int_t ***t;

  t = (int_t***) malloc((size_t) (nrow*sizeof(int_t**)));
  if (!t) {
    message("i3matrix()", "allocation failure 1", WARNING);
    return NULL;
  }
  t -= nrl;

  t[nrl] = (int_t**) malloc((size_t) (nrow*ncol*sizeof(int_t*)));
  if (!t[nrl]) {
    message("i3matrix()", "allocation failure 2", WARNING);
    return NULL;
  }
  t[nrl] -= ncl;

  t[nrl][ncl] = (int_t*) malloc((size_t) (nrow*ncol*ndep*sizeof(int_t)));
  if (!t[nrl][ncl]) {
    message("i3matrix()", "allocation failure 3", WARNING);
    return NULL;
  }
  t[nrl][ncl] -= ndl;

  for (j=ncl+1; j<=nch; j++)
    t[nrl][j] = t[nrl][j-1] + ndep;
  for (i=nrl+1; i<=nrh; i++) {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol*ndep;
    for (j=ncl+1; j<=nch; j++)
      t[i][j] = t[i][j-1] + ndep;
  }

  return t;
}


zomplex ***z3matrix (int_t nrl, int_t nrh,
		     int_t ncl, int_t nch,
		     int_t ndl, int_t ndh)
/* ------------------------------------------------------------------------- *
 * Allocate a zomplex 3-matrix with ranges [nrl..nrh][ncl..nch][ndl..ndh].
 * ------------------------------------------------------------------------- */
{
  int_t i, j, nrow = nrh-nrl+1, ncol = nch-ncl+1, ndep = ndh-ndl+1;
  zomplex ***t;

  t = (zomplex***) malloc((size_t) (nrow*sizeof(zomplex**)));
  if (!t) {
    message("z3matrix()", "allocation failure 1", WARNING);
    return NULL;
  }
  t -= nrl;

  t[nrl] = (zomplex**) malloc((size_t) (nrow*ncol*sizeof(zomplex*)));
  if (!t[nrl]) {
    message("z3matrix()", "allocation failure 2", WARNING);
    return NULL;
  }
  t[nrl] -= ncl;

  t[nrl][ncl] = (zomplex*) malloc((size_t) (nrow*ncol*ndep*sizeof(zomplex)));
  if (!t[nrl][ncl]) {
    message("z3matrix()", "allocation failure 3", WARNING);
    return NULL;
  }
  t[nrl][ncl] -= ndl;

  for (j=ncl+1; j<=nch; j++)
    t[nrl][j] = t[nrl][j-1] + ndep;
  for (i=nrl+1; i<=nrh; i++) {
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol*ndep;
    for (j=ncl+1; j<=nch; j++)
      t[i][j] = t[i][j-1] + ndep;
  }

  return t;
}


void freeCvector (complex *v, int_t nl)
/* ------------------------------------------------------------------------- *
 * Frees a complex vector allocated by cvector().
 * ------------------------------------------------------------------------- */
{
  free((FREE_ARG) (v+nl));
}


void freeDvector (double *v, int_t nl)
/* ------------------------------------------------------------------------- *
 * Frees a double vector allocated by dvector().
 * ------------------------------------------------------------------------- */
{
  free((FREE_ARG) (v+nl));
}


void freeSvector (float *v, int_t nl)
/* ------------------------------------------------------------------------- *
 * Frees a float vector allocated by svector().
 * ------------------------------------------------------------------------- */
{
  free((FREE_ARG) (v+nl));
}


void freeIvector (int_t *v, int_t nl)
/* ------------------------------------------------------------------------- *
 * Frees an int vector allocated by ivector().
 * ------------------------------------------------------------------------- */
{
  free((FREE_ARG) (v+nl));
}


void freeZvector (zomplex *v, int_t nl)
/* ------------------------------------------------------------------------- *
 * Frees a complex vector allocated by zvector().
 * ------------------------------------------------------------------------- */
{
  free((FREE_ARG) (v+nl));
}


void freeCmatrix (complex **m, int_t nrl, int_t ncl)
/* ------------------------------------------------------------------------- *
 * Frees a complex matrix allocated with cmatrix().
 * ------------------------------------------------------------------------- */
{
  free((FREE_ARG) (m[nrl]+ncl));
  free((FREE_ARG) (m+nrl));
}


void freeDmatrix (double **m, int_t nrl, int_t ncl)
/* ------------------------------------------------------------------------- *
 * Frees a double matrix allocated with dmatrix().
 * ------------------------------------------------------------------------- */
{
  free((FREE_ARG) (m[nrl]+ncl));
  free((FREE_ARG) (m+nrl));
}


void freeSmatrix (float **m, int_t nrl, int_t ncl)
/* ------------------------------------------------------------------------- *
 * Frees a float matrix allocated with smatrix().
 * ------------------------------------------------------------------------- */
{
  free((FREE_ARG) (m[nrl]+ncl));
  free((FREE_ARG) (m+nrl));
}


void freeImatrix (int_t **m, int_t nrl, int_t ncl)
/* ------------------------------------------------------------------------- *
 * Frees an int matrix allocated with imatrix().
 * ------------------------------------------------------------------------- */
{
  free((FREE_ARG) (m[nrl]+ncl));
  free((FREE_ARG) (m+nrl));
}


void freeZmatrix (zomplex **m, int_t nrl, int_t ncl)
/* ------------------------------------------------------------------------- *
 * Frees a zomplex matrix allocated with zmatrix().
 * ------------------------------------------------------------------------- */
{
  free((FREE_ARG) (m[nrl]+ncl));
  free((FREE_ARG) (m+nrl));
}


void freeC3matrix (complex ***t, int_t nrl, int_t ncl, int_t ndl)
/* ------------------------------------------------------------------------- *
 * Frees a complex 3-matrix allocated with c3matrix().
 * ------------------------------------------------------------------------- */
{
  free ((FREE_ARG) (t[nrl][ncl]+ndl));
  free ((FREE_ARG) (t[nrl]+ncl));
  free ((FREE_ARG) (t+ndl));
}


void freeD3matrix (double ***t, int_t nrl, int_t ncl, int_t ndl)
/* ------------------------------------------------------------------------- *
 * Frees a double 3-matrix allocated with d3matrix().
 * ------------------------------------------------------------------------- */
{
  free ((FREE_ARG) (t[nrl][ncl]+ndl));
  free ((FREE_ARG) (t[nrl]+ncl));
  free ((FREE_ARG) (t+ndl));
}


void freeS3matrix (float ***t, int_t nrl, int_t ncl, int_t ndl)
/* ------------------------------------------------------------------------- *
 * Frees a float 3-matrix allocated with s3matrix().
 * ------------------------------------------------------------------------- */
{
  free ((FREE_ARG) (t[nrl][ncl]+ndl));
  free ((FREE_ARG) (t[nrl]+ncl));
  free ((FREE_ARG) (t+ndl));
}


void freeI3matrix (int_t ***t, int_t nrl, int_t ncl, int_t ndl)
/* ------------------------------------------------------------------------- *
 * Frees an int 3-matrix allocated with i3matrix().
 * ------------------------------------------------------------------------- */
{
  free ((FREE_ARG) (t[nrl][ncl]+ndl));
  free ((FREE_ARG) (t[nrl]+ncl));
  free ((FREE_ARG) (t+ndl));
}


void freeZ3matrix (zomplex ***t, int_t nrl, int_t ncl, int_t ndl)
/* ------------------------------------------------------------------------- *
 * Frees a zomplex 3-matrix allocated with z3matrix().
 * ------------------------------------------------------------------------- */
{
  free ((FREE_ARG) (t[nrl][ncl]+ndl));
  free ((FREE_ARG) (t[nrl]+ncl));
  free ((FREE_ARG) (t+ndl));
}
