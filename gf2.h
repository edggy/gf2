#ifndef GF2_H_   /* Include guard */
#define GF2_H_

#include <stdbool.h>

typedef unsigned long long GF2;

unsigned int gf2bitlength(GF2 a);
GF2* gf2divmodvect(GF2 avec[], GF2 dvec[], size_t len);
GF2* gf2exteuc(GF2 a, GF2 b, GF2* brow);

GF2 gf2add(GF2 a, GF2 b, unsigned int lgsize);
GF2 gf2sub(GF2 a, GF2 b, unsigned int lgsize);
GF2 gf2mulmod(GF2 a, GF2 b, GF2 m);
GF2 gf2powmod(GF2 a, GF2 k, GF2 m);
GF2 gf2divmod(GF2 a, GF2 d, GF2 m, bool* err);

GF2 gf2mod(GF2 a, GF2 m);
GF2 gf2modinv(GF2 x, GF2 m, bool* err);

#endif
