#ifndef FLOWRATE_H
#define FLOWRATE_H


class Flowrate
// ===========================================================================
// Routines to measure and set volumetric flow rates using Stokes
// Green's functions.
// ===========================================================================
{
public:
  Flowrate (Domain*, FEML*);
  ~Flowrate () { }

  real_t getQ () const;
  real_t setQ (AuxField*, AuxField*) const;

protected:
  Domain*           _src;
  vector<Edge*>     _curve;
  vector<AuxField*> _green;
  real_t            _gamma;
  real_t            _refQ ;

};

#endif
