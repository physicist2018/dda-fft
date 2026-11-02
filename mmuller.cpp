
#include "cplx.h"
#include "mmuller.hpp"
#include <stdio.h>

/********************************* Operators ********************************/

// for Muller matrix

mmuller operator+(mmuller mm0, mmuller mm1)
{   mmuller mmTmp;

	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			mmTmp.e[i][j]=mm0.e[i][j]+mm1.e[i][j];

	return mmTmp;
}

mmuller operator-(mmuller mm0, mmuller mm1)
{   mmuller mmTmp;

	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			mmTmp.e[i][j]=mm0.e[i][j]-mm1.e[i][j];

	return mmTmp;
}

mmuller operator-(mmuller mm0)
{   mmuller mmTmp;

	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			mmTmp.e[i][j]=-mm0.e[i][j];

	return mmTmp;
}

mmuller operator*(real r0, mmuller mm0)
{   mmuller mmTmp;

	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			mmTmp.e[i][j]=r0*mm0.e[i][j];

	return mmTmp;
}

mmuller operator*(mmuller mm0, real r0)
{   mmuller mmTmp;

	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			mmTmp.e[i][j]=r0*mm0.e[i][j];

	return mmTmp;
}

mmuller operator/(mmuller mm0, real r0)
{   mmuller mmTmp;

	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			mmTmp.e[i][j]=mm0.e[i][j]/r0;

	return mmTmp;
}

mmuller mmNull()
{   mmuller mmTmp;

	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			mmTmp.e[i][j]=0.;

	return mmTmp;
}

// for Jones matrix

mjones operator+(mjones mj0, mjones mj1)
{   mjones mjTmp;

	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++)
			mjTmp.e[i][j]=mj0.e[i][j]+mj1.e[i][j];

	return mjTmp;
}

mjones operator-(mjones mj0, mjones mj1)
{   mjones mjTmp;

	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++)
			mjTmp.e[i][j]=mj0.e[i][j]-mj1.e[i][j];

	return mjTmp;
}

mjones operator-(mjones mj0)
{   mjones mjTmp;

	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++)
			mjTmp.e[i][j]=-mj0.e[i][j];

	return mjTmp;
}

mjones operator*(real r0, mjones mj0)
{   mjones mjTmp;

	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++)
			mjTmp.e[i][j]=r0*mj0.e[i][j];

	return mjTmp;
}

mjones operator*(mjones mj0, real r0)
{   mjones mjTmp;

	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++)
			mjTmp.e[i][j]=r0*mj0.e[i][j];

	return mjTmp;
}

mjones operator*(cplx c0, mjones mj0)
{   mjones mjTmp;

	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++)
			mjTmp.e[i][j]=c0*mj0.e[i][j];

	return mjTmp;
}

mjones operator*(mjones mj0, cplx c0)
{   mjones mjTmp;

	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++)
			mjTmp.e[i][j]=c0*mj0.e[i][j];

	return mjTmp;
}

mjones mjNull()
{   mjones mjTmp;

	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++)
		{   mjTmp.e[i][j].r=0.;
			mjTmp.e[i][j].i=0.;
		}

	return mjTmp;
}

mjones mjUnity()
{   mjones mjTmp;

	mjTmp.e[0][0].r=1.; mjTmp.e[0][0].i=0.;
	mjTmp.e[0][1].r=0.; mjTmp.e[0][1].i=0.;
	mjTmp.e[1][0].r=0.; mjTmp.e[1][0].i=0.;
	mjTmp.e[1][1].r=1.; mjTmp.e[1][1].i=0.;

	return mjTmp;
}

mjones BuildMatrixJones(vector kinc, vector ksc, vector vNormal,
						cvector Einc0, cvector Einc1, cvector Esc0, cvector Esc1)
{   mjones mjTmp;
	cplx i={0.,1.},cDen;
	cplx c0,c1,c2,c3;
	cplx cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc7;
	vector v0,vPlaneIn,vPlaneSc;

	v0=kinc^vNormal;    vPlaneIn=(1./abs(v0))*(v0);
	v0=ksc^vNormal;     vPlaneSc=(-1./abs(v0))*(v0);

	DCmul(vNormal.x,Einc0.x,c0);
	DCmul(vNormal.y,Einc0.y,c1);
	DCmul(vNormal.z,Einc0.z,c2);
	Cadd(c0,c1,c3); Cadd(c2,c3,cc0);

	DCmul(vNormal.x,Einc1.x,c0);
	DCmul(vNormal.y,Einc1.y,c1);
	DCmul(vNormal.z,Einc1.z,c2);
	Cadd(c0,c1,c3); Cadd(c2,c3,cc1);

	DCmul(vPlaneIn.x,Einc0.x,c0);
	DCmul(vPlaneIn.y,Einc0.y,c1);
	DCmul(vPlaneIn.z,Einc0.z,c2);
	Cadd(c0,c1,c3); Cadd(c2,c3,cc2);

	DCmul(vPlaneIn.x,Einc1.x,c0);
	DCmul(vPlaneIn.y,Einc1.y,c1);
	DCmul(vPlaneIn.z,Einc1.z,c2);
	Cadd(c0,c1,c3); Cadd(c2,c3,cc3);

	Cmul(cc0,cc3,c0);
	Cmul(cc1,cc2,c1);
	Csub(c0,c1,c3);
	cDen.r=-c3.i;   cDen.i=c3.r;

	DCmul(vNormal.x,Esc0.x,c0);
	DCmul(vNormal.y,Esc0.y,c1);
	DCmul(vNormal.z,Esc0.z,c2);
	Cadd(c0,c1,c3); Cadd(c2,c3,cc4);

	DCmul(vNormal.x,Esc1.x,c0);
	DCmul(vNormal.y,Esc1.y,c1);
	DCmul(vNormal.z,Esc1.z,c2);
	Cadd(c0,c1,c3); Cadd(c2,c3,cc5);

	DCmul(vPlaneSc.x,Esc0.x,c0);
	DCmul(vPlaneSc.y,Esc0.y,c1);
	DCmul(vPlaneSc.z,Esc0.z,c2);
	Cadd(c0,c1,c3); Cadd(c2,c3,cc6);

	DCmul(vPlaneSc.x,Esc1.x,c0);
	DCmul(vPlaneSc.y,Esc1.y,c1);
	DCmul(vPlaneSc.z,Esc1.z,c2);
	Cadd(c0,c1,c3); Cadd(c2,c3,cc7);

	Cmul(cc7,cc0,c0);
	Cmul(cc6,cc1,c1);
	Csub(c0,c1,c2);
	Cdiv(c2,cDen,mjTmp.e[0][0]);

	Cmul(cc6,cc3,c0);
	Cmul(cc7,cc2,c1);
	Csub(c0,c1,c2);
	Cdiv(c2,cDen,mjTmp.e[0][1]);

	Cmul(cc5,cc0,c0);
	Cmul(cc4,cc1,c1);
	Csub(c0,c1,c2);
	Cdiv(c2,cDen,mjTmp.e[1][0]);

	Cmul(cc4,cc3,c0);
	Cmul(cc5,cc2,c1);
	Csub(c0,c1,c2);
	Cdiv(c2,cDen,mjTmp.e[1][1]);

	return mjTmp;
}

mmuller JonesToMuller(mjones mJones)
{   mmuller mmTmp;

	mmTmp.e[0][0]=.5*re(mJones.e[0][0]*conj(mJones.e[0][0])
						+mJones.e[0][1]*conj(mJones.e[0][1])
						+mJones.e[1][0]*conj(mJones.e[1][0])
						+mJones.e[1][1]*conj(mJones.e[1][1]));

	mmTmp.e[0][1]=.5*re(mJones.e[0][0]*conj(mJones.e[0][0])
						-mJones.e[0][1]*conj(mJones.e[0][1])
						+mJones.e[1][0]*conj(mJones.e[1][0])
						-mJones.e[1][1]*conj(mJones.e[1][1]));

	mmTmp.e[0][2]=re(mJones.e[0][0]*conj(mJones.e[0][1])
					 +mJones.e[1][1]*conj(mJones.e[1][0]));

	mmTmp.e[0][3]=im(mJones.e[0][0]*conj(mJones.e[0][1])
					 -mJones.e[1][1]*conj(mJones.e[1][0]));

	mmTmp.e[1][0]=.5*re(mJones.e[0][0]*conj(mJones.e[0][0])
						+mJones.e[0][1]*conj(mJones.e[0][1])
						-mJones.e[1][0]*conj(mJones.e[1][0])
						-mJones.e[1][1]*conj(mJones.e[1][1]));

	mmTmp.e[1][1]=.5*re(mJones.e[0][0]*conj(mJones.e[0][0])
						-mJones.e[0][1]*conj(mJones.e[0][1])
						-mJones.e[1][0]*conj(mJones.e[1][0])
						+mJones.e[1][1]*conj(mJones.e[1][1]));

	mmTmp.e[1][2]=re(mJones.e[0][0]*conj(mJones.e[0][1])
					 -mJones.e[1][1]*conj(mJones.e[1][0]));

	mmTmp.e[1][3]=im(mJones.e[0][0]*conj(mJones.e[0][1])
					 +mJones.e[1][1]*conj(mJones.e[1][0]));

	mmTmp.e[2][0]=re(mJones.e[0][0]*conj(mJones.e[1][0])
					 +mJones.e[1][1]*conj(mJones.e[0][1]));

	mmTmp.e[2][1]=re(mJones.e[0][0]*conj(mJones.e[1][0])
					 -mJones.e[1][1]*conj(mJones.e[0][1]));

	mmTmp.e[2][2]=re(mJones.e[1][1]*conj(mJones.e[0][0])
					 +mJones.e[0][1]*conj(mJones.e[1][0]));

	mmTmp.e[2][3]=im(mJones.e[0][0]*conj(mJones.e[1][1])
					 +mJones.e[1][0]*conj(mJones.e[0][1]));

	mmTmp.e[3][0]=im(mJones.e[1][0]*conj(mJones.e[0][0])
					 +mJones.e[1][1]*conj(mJones.e[0][1]));

	mmTmp.e[3][1]=im(mJones.e[1][0]*conj(mJones.e[0][0])
					 -mJones.e[1][1]*conj(mJones.e[0][1]));

	mmTmp.e[3][2]=im(mJones.e[1][1]*conj(mJones.e[0][0])
					 -mJones.e[0][1]*conj(mJones.e[1][0]));

	mmTmp.e[3][3]=re(mJones.e[1][1]*conj(mJones.e[0][0])
					 -mJones.e[0][1]*conj(mJones.e[1][0]));

	return mmTmp;
}
