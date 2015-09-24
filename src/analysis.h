#ifndef ANALYSIS_H
#define ANALYSIS_H

class Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow
// solver.  This is designed to be overridden at implementation level
// if needed.
// ===========================================================================
{
public:
  Analyser  (Domain*, FEML*);
  ~Analyser () { }

  void analyse (AuxField**, AuxField**);

protected:
  Domain*               _src      ; // Source information.
  ofstream              _par_strm ; // File for particle tracking.
  ofstream              _his_strm ; // File for history points.
  ofstream              _mdl_strm ; // File for modal energies.
  vector<HistoryPoint*> _history  ; // Locations, etc. of history points.
  list<FluidParticle*>  _particle ; // List of fluid particles.
  vector<Point*>        _initial  ; // Starting locations of particles.
  Statistics*           _stats    ; // Field average statistics.
  Statistics*           _ph_stats ; // Phase-average field statistics.

  void modalEnergy ();
  void divergence  (AuxField**) const;
  void estimateCFL ()           const;
};

#endif
