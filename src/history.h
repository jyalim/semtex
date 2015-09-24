#ifndef HISTORY_H
#define HISTORY_H


class HistoryPoint
// ===========================================================================
// Class used to provide x,y,z history information.
// ===========================================================================
{
public:
  HistoryPoint (const int_t id, const Element* e, const real_t r, 
		const real_t s, const real_t x, const real_t y, const real_t z):
    _id (id), _E (e), _r (r), _s (s), _x (x), _y (y), _z (z) { }

  int_t                 ID      () const { return _id; } 
  void                  extract (vector<AuxField*>&, real_t*) const;
  static const Element* locate  (const real_t, const real_t,
				 vector<Element*>&, real_t&, real_t&);

private:
  const int_t    _id;		// Numeric identifier.
  const Element* _E ;		// Pointer to element.
  const real_t   _r ;		// Canonical-space r-location.
  const real_t   _s ;		// Canonical-space s-location.
  const real_t   _x ;		// x location.
  const real_t   _y ;		// y location.
  const real_t   _z ;		// Location in homogeneous direction.
};

#endif
