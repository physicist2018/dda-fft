
#ifndef _CPLX_ESZ_DEFINED_

#include <math.h>
#include <stdlib.h>
#define _CPLX_ESZ_DEFINED_

//typedef double real;
//typedef long double real;
typedef float real;

typedef
	struct
	{   long double r;
		long double i;
	}
	ld_cplx;

typedef
	struct
	{   real r;
		real i;
	}
	cplx;

#define PI 3.1415926535897932384626433832795


// Float type of components of complex numbers

cplx operator+(cplx c0, cplx c1);
cplx operator+(real r0, cplx c0);
cplx operator+(cplx c0, real r0);
cplx operator-(cplx c0, cplx c1);
cplx operator-(real r0, cplx c0);
cplx operator-(cplx c0, real r0);
cplx operator-(cplx c0);
cplx operator*(cplx c0, cplx c1);
cplx operator*(real r0, cplx c0);
cplx operator*(cplx c0, real r0);
cplx operator/(cplx c0, cplx c1);
cplx operator/(real r0, cplx c0);
cplx operator/(cplx c0, real r0);
cplx conj(cplx c0);
cplx root(real n, cplx c0);
void root(real n, cplx c0, cplx * cRootArray);
real re(cplx c0);
real im(cplx c0);
real abs(cplx c0);
cplx exp(cplx c0);
cplx log(cplx c0);
cplx exp_i(real r0);
cplx cos(cplx c0);
cplx sin(cplx c0);
cplx sqrt(cplx c0);
cplx zero();
int rounding(real r0);


// Fixed type (long double) of components of complex numbers

ld_cplx operator+(ld_cplx c0, ld_cplx c1);
ld_cplx operator-(ld_cplx c0, ld_cplx c1);
ld_cplx operator-(ld_cplx c0);
ld_cplx operator*(ld_cplx c0, ld_cplx c1);
ld_cplx operator/(ld_cplx c0, ld_cplx c1);

#endif
