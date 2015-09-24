#ifndef INTEGRATION_H
#define INTEGRATION_H

class Integration
// ===========================================================================
// Return coefficients for time-integration schemes.
// ===========================================================================
{
public:
  static const int_t OrderMax;

  static void AdamsBashforth (const int_t, real_t*);
  static void AdamsMoulton   (const int_t, real_t*);
  static void StifflyStable  (const int_t, real_t*);
  static void Extrapolation  (const int_t, real_t*);
};

#endif
