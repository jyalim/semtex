#ifndef STATISTICS_H
#define STATISTICS_H

class Statistics
// ===========================================================================
// Routines for statistical analysis of AuxFields.
// ===========================================================================
{
friend ifstream& operator >> (ifstream&, Statistics&);
friend ofstream& operator << (ofstream&, Statistics&);
public:
  Statistics (Domain*);
  
  void initialise  (const char*);
  void update      (AuxField**, AuxField**); // -- Two work arrays supplied.
  void dump        (const char*);

  void phaseUpdate (const int_t, AuxField**, AuxField**);

protected:
  const char*          _name;
  Domain*              _base;	// -- Local pointer to external Domain.
  map<char, AuxField*> _raw ;	// -- Pointers to the base/raw storage areas.
  map<char, AuxField*> _avg ;	// -- Map names to averaging buffers.
  int_t                _navg;	// -- Number of averages so far.
  int_t                _iavg;	// -- Same as value of token "AVERAGE".
  int_t                _nraw;	// -- No. of raw field variables in _base.
  int_t                _nvel;	// -- No. of velocity components.
  int_t                _nrey;	// -- No. of Reynolds stress correlations.
  int_t                _neng;	// -- No. additional correlations for energy.
};

#endif
