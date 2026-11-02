
#ifndef _MJONES_ESZ_DEFINED_

#include "cvector.hpp"
#define _MJONES_ESZ_DEFINED_

typedef
	struct
	{   cplx e[2][2];
	}   mjones;

typedef
	struct
	{   real e[4][4];
	}   mmuller;

/********************************* Operators ********************************/

// for Muller matrix

mmuller operator+(mmuller mm0, mmuller mm1);
mmuller operator-(mmuller mm0, mmuller mm1);
mmuller operator-(mmuller mm0);
mmuller operator*(real r0, mmuller mm0);
mmuller operator*(mmuller mm0, real r0);
mmuller operator/(mmuller mm0, real r0);
mmuller mmNull();

// for Jones matrix

mjones operator+(mjones mj0, mjones mj1);
mjones operator-(mjones mj0, mjones mj1);
mjones operator-(mjones mj0);
mjones operator*(real r0, mjones mj0);
mjones operator*(mjones mj0, real r0);
mjones operator*(cplx c0, mjones mj0);
mjones operator*(mjones mj0, cplx c0);
mjones mjNull();
mjones mjUnity();
mjones BuildMatrixJones(vector kinc, vector ksc, vector vNormal,
						cvector Einc0, cvector Einc1, cvector Esc0, cvector Esc1);
mmuller JonesToMuller(mjones mJones);

#endif
