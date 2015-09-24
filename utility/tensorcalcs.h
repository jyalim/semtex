#ifndef TENSORCALCS_H
#define TENSORCALCS_H

/*****************************************************************************
 *          PROTOTYPES FOR 3^i TENSOR & VECTOR CALCULATIONS
 *****************************************************************************/

#ifdef __cplusplus
extern "C" {
  namespace tensor3 {
#endif

enum {UFC, UNS, SFS, SNS};

void   transpose           (const double  T[9], double O[9]);
void   transpose_in_place  (double  T[9]);
void   symmetric_part      (const double  T[9], double S[9]);
void   anti_symmetric_part (const double  T[9], double A[9]);
double contracted_product  (const double  T[9]);
void   S2plusA2            (const double V[9], double O[9]);
double enstrophy           (const double VG[9]);
double dissipation         (const double VG[9]);
double strainrate          (const double VG[9]);
double trace               (const double  T[9]);
void   remove_trace        (double  T[9]);
double det                 (const double  T[9]);
void   invariants          (const double  T[9], double*, double*, double*);
void   real_eigenvalues    (double a2, double a1, double a0, 
                            double  *z1, double *z2, double *z3);
void   eigenvector         (const double  T[9], double z, double e[3]);
void   vec                 (const double  A[9], double V[3]);
void   scale_vect          (double  v[3], double s);
void   vorticity           (const double  T[9], double w[3]);
double dot                 (const double  a[3], const double b[3]);
void   cross               (const double  a[3], const double b[3], double c[3]);
void   normalize           (double  a[3]);
double classifyi           (const double VG[9], int *topology);
double discrimi            (const double VG[9]);
double lambda2             (const double VG[9]);
void   lambvector          (const double VG[9], const double V[3], double L[3]);
double helicity            (const double VG[9], const double V[3]);

#ifdef __cplusplus
  }
}
#endif

#endif
