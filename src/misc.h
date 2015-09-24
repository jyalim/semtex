#ifndef MISC_H
#define MISC_H


// -- Routines from misc.C:

ostream& printVector (ostream&, const char*, const int_t, ... );
char*    upperCase   (char *);
void     writeField  (ostream&, const char*, const int_t, const real_t,
		      vector<AuxField*>&);
void     readField   (istream&, vector<AuxField*>&);
#endif
