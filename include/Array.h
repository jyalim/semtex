#ifndef ARRAY_H
#define ARRAY_H
//////////////////////////////////////////////////////////////////////////////
// Array.h: templated 1D (vector) and 2D (matrix) array classes.
// 
// No subscript checking, and direct access to the underlying pointers
// can be obtained if desired using operator ().
//
// Matrix storage is row-major contiguous.
//
// $Id: Array.h,v 8.1 2015/04/20 11:14:14 hmb Exp $
//////////////////////////////////////////////////////////////////////////////


template<class T>
class vector
// ===========================================================================
// 1D array class.
// ===========================================================================
{
public:
  // -- Creation/destruction.

  vector (const long n)         { num_elts = n;
                                  data = (T*) (n) ? new T[(size_t) n] : 0; }
  vector ()                     { num_elts = 0; data = 0; } 
  vector (const vector<T>& src) { num_elts = src.num_elts;
				  data = (T*) new T[(size_t) num_elts];
				  copy (src); } 
  ~vector ()                    { delete [] data; }

  // -- Assignment.

  vector<T>& operator = (const vector<T>& src) {
    if (data != src.data ) { setSize (src.num_elts); copy (src); }
    return *this;
  }

  vector<T>& operator = (const T& src) {
    register long i;
    for (i = 0; i < num_elts; i++) data[i] = src;
    return *this;
  }

  // -- Subscripting.

        T* operator () ()       { return data; }
  const T* operator () () const { return data; } 

        T& operator [] (const long i)       { return data[i]; }
  const T& operator [] (const long i) const { return data[i]; }
        T& operator () (const long i)       { return data[i]; }
  const T& operator () (const long i) const { return data[i]; } 

  // -- Size operators/information.

  long getSize () const       { return num_elts; }
  void setSize (const long n) {
    if (n != num_elts) 
      {
	delete [] data;
	num_elts = n;
	data     = (T*) (n) ? new T[(size_t) n] : 0; 
      }
  }

  // -- Provided for compatibility/interchangability with STL.

  long size() const          { return getSize(); }
  void resize (const long n) { setSize(n); }

private:
  long num_elts;
  T*   data;
  void copy (const vector<T>& src) {
    register long i;
    for (i = 0; i < num_elts; i++) data[i] = src.data[i];
  }
};


template<class T> class matrix
// ==========================================================================
// 2D contiguous-storage row-major array.
// ==========================================================================
{
public:
  // -- Creation/destruction.

  matrix (const long n_rows, const long n_cols) {
    nr   = n_rows;
    nc   = n_cols;
    row  = (T**) new T* [(size_t) nr];
    data = (T*)  new T  [(size_t) (nr * nc)];
    for (register long i = 0; i < nr; i++) row[i] = data + i * nc;
  }
  matrix () { 
    nr   = nc = 0;
    row  = 0;
    data = 0;
  } 
  matrix (const matrix<T>& src) { 
    nr   = src.nr;
    nc   = src.nc;
    row  = (T**) new T* [(size_t) nr];
    data = (T*)  new T  [(size_t) (nr * nc)];
    for (register long i = 0; i < nr; i++) row[i] = data + i * nc;
    copy (src); 
  } 
  ~matrix () {
    delete [] data;
    delete [] row;
  }

  // -- Assignment.

  matrix<T>& operator = (const matrix<T>& src) {
    if (data != src.data ) {
      setSize (src.nr, src.nc);
      copy    (src); }
    return *this;
  }

  matrix<T>& operator = (const T& src) {
    const long nelts = nr * nc;
    for (register long i = 0; i < nelts; i++) data[i] = src;
    return *this;
  }

  // -- Subscripting.

        T** operator () ()                          { return             row; }
  const T** operator () () const                    { return (const T**) row; }
  
        T* operator () (const long i)                    { return row[i]; }
  const T* operator () (const long i) const              { return row[i]; }
  
        T& operator () (const long i, const long j)      {return data[j+i*nc];}
  const T& operator () (const long i, const long j) const{return data[j+i*nc];}

  // -- Size operators/information.

  long nRows   () const      { return nr;      }
  long nCols   () const      { return nc;      }
  long getSize () const      { return nr * nc; }
  void setSize (const long n_rows, const long n_cols) { 
    if (n_rows != nr && n_cols != nc) 
      {
	delete [] data;
	delete [] row;
	nr   = n_rows;
	nc   = n_cols;
	row  = (T**) (nr*nc) ? new T* [(size_t) nr]      : 0;
	data = (T*)  (nr*nc) ? new T  [(size_t) (nr*nc)] : 0;
	for (register long i = 0; i < nr; i++) row[i] = data + i * nc;
      }
  }

private:
  long nr;
  long nc;
  T**  row;
  T*   data;
  void copy (const matrix<T>& src) {
    const long nelts = nr * nc;
    for (register long i = 0; i < nelts; i++) data[i] = src.data[i];
  }
};

#endif
