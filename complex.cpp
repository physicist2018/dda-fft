
#include "complex.hpp"
#include "math.h"


// Float type of components of complex numbers

/******************************* Operators **********************************/

/******************************* operator+ **********************************/

cplx operator+(cplx c0, cplx c1)
{   cplx cTmp;
    cTmp.r=c0.r+c1.r;    cTmp.i=c0.i+c1.i;
    return cTmp;
}

cplx operator+(real r0, cplx c0)
{   cplx cTmp;
    cTmp.r=r0+c0.r;    cTmp.i=c0.i;
    return cTmp;
}

cplx operator+(cplx c0, real r0)
{   cplx cTmp;
    cTmp.r=r0+c0.r;    cTmp.i=c0.i;
    return cTmp;
}

/******************************* operator- **********************************/

cplx operator-(cplx c0, cplx c1)
{   cplx cTmp;
    cTmp.r=c0.r-c1.r;    cTmp.i=c0.i-c1.i;
    return cTmp;
}

cplx operator-(real r0, cplx c0)
{   cplx cTmp;
    cTmp.r=r0-c0.r;    cTmp.i=-c0.i;
    return cTmp;
}

cplx operator-(cplx c0, real r0)
{   cplx cTmp;
    cTmp.r=c0.r-r0;    cTmp.i=c0.i;
    return cTmp;
}

cplx operator-(cplx c0)
{   cplx cTmp;
    cTmp.r=-c0.r;    cTmp.i=-c0.i;
    return cTmp;
}

/******************************* operator* **********************************/

cplx operator*(cplx c0, cplx c1)
{   cplx cTmp;
    cTmp.r=c0.r*c1.r-c0.i*c1.i;    cTmp.i=c0.r*c1.i+c0.i*c1.r;
    return cTmp;
}

cplx operator*(real r0, cplx c0)
{   cplx cTmp;
    cTmp.r=r0*c0.r;    cTmp.i=r0*c0.i;
    return cTmp;
}

cplx operator*(cplx c0, real r0)
{   cplx cTmp;
    cTmp.r=r0*c0.r;    cTmp.i=r0*c0.i;
    return cTmp;
}

/********************************* operator/ ********************************/

cplx operator/(cplx c0, cplx c1)
{   real c1_norma;
    cplx cTmp;

    c1_norma=c1.r*c1.r+c1.i*c1.i;
    cTmp.r=(c0.r*c1.r+c0.i*c1.i)/c1_norma;
    cTmp.i=(c0.i*c1.r-c0.r*c1.i)/c1_norma;
    return cTmp;
}

cplx operator/(real r0, cplx c0)
{   real c0_norma;
    cplx cTmp;

    c0_norma=c0.r*c0.r+c0.i*c0.i;
    cTmp.r=( r0*c0.r)/c0_norma;
    cTmp.i=(-r0*c0.i)/c0_norma;
    return cTmp;
}

cplx operator/(cplx c0, real r0)
{   cplx cTmp;

    cTmp.r=c0.r/r0;
    cTmp.i=c0.i/r0;
    return cTmp;
}

/************************* Commonplace functions ****************************/

cplx conj(cplx c0)
{   cplx cTmp;

    cTmp.r= c0.r;
    cTmp.i=-c0.i;
    return cTmp;
}

cplx root(real n, cplx c0)
{   real length=pow((real)abs(c0),(real)(1./n)),phase=atan2(c0.i,c0.r)/n;
	cplx cTmp;

	cTmp.r=length*cos(phase);
	cTmp.i=length*sin(phase);
	return cTmp;
}

void root(real n, cplx c0, cplx * cRootArray)
{   int i;
	real length=pow((real)abs(c0),(real)(1./n)),phase=atan2(c0.i,c0.r)/n,r=2.*PI/n;

	for(i=0;i<(int)n;i++)
	{   cRootArray[i].r=length*cos(phase+r*i);
		cRootArray[i].i=length*sin(phase+r*i);
	}

	return;
}

real re(cplx c0)
{   real rTmp;
    rTmp=c0.r;
    return rTmp;
}

real im(cplx c0)
{   real rTmp;
    rTmp=c0.i;
    return rTmp;
}

real abs(cplx c0)
{   real rTmp;
    rTmp=sqrt(c0.r*c0.r+c0.i*c0.i);
    return rTmp;
}

cplx exp(cplx c0)
{   real r0;
    cplx cTmp;

    r0=exp(c0.r);
    cTmp.r=r0*cos(c0.i);
    cTmp.i=r0*sin(c0.i);
    return cTmp;
}

cplx log(cplx c0)
{   real r0;
    cplx cTmp;

	r0=abs(c0);

    cTmp.r=log(r0);
    cTmp.i=atan2(c0.i,c0.r);
    return cTmp;
}

cplx exp_i(real r0)
{   cplx cTmp;
    cTmp.r=cos(r0);
    cTmp.i=sin(r0);
    return cTmp;
}

cplx cos(cplx c0)
{   cplx cTmp;
    cTmp.r= cos(c0.r)*cosh(c0.i);
    cTmp.i=-sin(c0.r)*sinh(c0.i);
    return cTmp;
}

cplx sin(cplx c0)
{   cplx cTmp;
    cTmp.r=sin(c0.r)*cosh(c0.i);
    cTmp.i=cos(c0.r)*sinh(c0.i);
    return cTmp;
}

cplx sqrt(cplx c0)
{   real r0,r1;
    cplx cTmp;

	r0=sqrt(abs(c0));
    r1=atan2(c0.i,c0.r)/2.;
    cTmp.r=cos(r1)*r0;
    cTmp.i=sin(r1)*r0;
    return cTmp;
}

cplx zero()
{   cplx cTmp={0.,0.};
	return cTmp;
}

int rounding(real r0)
{   int iTmp;
	double r1;

	if(modf(r0,&r1)>0.5)
		iTmp=(int)r1+1;
	else
		iTmp=(int)r1;

	return iTmp;
}


// Fixed type (long double) of components of complex numbers

ld_cplx operator+(ld_cplx c0, ld_cplx c1)
{   ld_cplx cTmp;
	cTmp.r=c0.r+c1.r;    cTmp.i=c0.i+c1.i;
	return cTmp;
}

ld_cplx operator-(ld_cplx c0, ld_cplx c1)
{   ld_cplx cTmp;
	cTmp.r=c0.r-c1.r;    cTmp.i=c0.i-c1.i;
	return cTmp;
}

ld_cplx operator-(ld_cplx c0)
{   ld_cplx cTmp;
	cTmp.r=-c0.r;    cTmp.i=-c0.i;
	return cTmp;
}

ld_cplx operator*(ld_cplx c0, ld_cplx c1)
{   ld_cplx cTmp;
	cTmp.r=c0.r*c1.r-c0.i*c1.i;    cTmp.i=c0.r*c1.i+c0.i*c1.r;
	return cTmp;
}

ld_cplx operator/(ld_cplx c0, ld_cplx c1)
{   long double c1_norma;
	ld_cplx cTmp;

	c1_norma=c1.r*c1.r+c1.i*c1.i;
	cTmp.r=(c0.r*c1.r+c0.i*c1.i)/c1_norma;
	cTmp.i=(c0.i*c1.r-c0.r*c1.i)/c1_norma;
	return cTmp;
}
