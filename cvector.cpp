
#include "cvector.hpp"

/********************************* Operators ********************************/

/********************************* operator+ ********************************/

vector operator+(vector v0, vector v1)
{   vector vTmp;
    vTmp.x=v0.x+v1.x;
    vTmp.y=v0.y+v1.y;
    vTmp.z=v0.z+v1.z;
    return vTmp;
}

cvector operator+(cvector cv0, vector v0)
{   cvector cvTmp;
    cvTmp.x=cv0.x+v0.x;
    cvTmp.y=cv0.y+v0.y;
    cvTmp.z=cv0.z+v0.z;
    return cvTmp;
}

cvector operator+(vector v0, cvector cv0)
{   cvector cvTmp;
    cvTmp.x=cv0.x+v0.x;
    cvTmp.y=cv0.y+v0.y;
    cvTmp.z=cv0.z+v0.z;
    return cvTmp;
}

cvector operator+(cvector cv0, cvector cv1)
{   cvector cvTmp;
    cvTmp.x=cv0.x+cv1.x;
    cvTmp.y=cv0.y+cv1.y;
    cvTmp.z=cv0.z+cv1.z;
    return cvTmp;
}

/********************************* operator- ********************************/

vector operator-(vector v0, vector v1)
{   vector vTmp;
    vTmp.x=v0.x-v1.x;
    vTmp.y=v0.y-v1.y;
    vTmp.z=v0.z-v1.z;
    return vTmp;
}

cvector operator-(cvector cv0, vector v0)
{   cvector cvTmp;
    cvTmp.x=cv0.x-v0.x;
    cvTmp.y=cv0.y-v0.y;
    cvTmp.z=cv0.z-v0.z;
    return cvTmp;
}

cvector operator-(vector v0, cvector cv0)
{   cvector cvTmp;
    cvTmp.x=v0.x-cv0.x;
    cvTmp.y=v0.y-cv0.y;
    cvTmp.z=v0.z-cv0.z;
    return cvTmp;
}

cvector operator-(cvector cv0, cvector cv1)
{   cvector cvTmp;
    cvTmp.x=cv0.x-cv1.x;
    cvTmp.y=cv0.y-cv1.y;
    cvTmp.z=cv0.z-cv1.z;
    return cvTmp;
}

vector operator-(vector v0)
{   vector vTmp;
    vTmp.x=-v0.x;
    vTmp.y=-v0.y;
    vTmp.z=-v0.z;
    return vTmp;
}

cvector operator-(cvector cv0)
{   cvector cvTmp;
    cvTmp.x=-cv0.x;
    cvTmp.y=-cv0.y;
    cvTmp.z=-cv0.z;
    return cvTmp;
}

/********************************* operator* ********************************/

real operator*(vector v0, vector v1)
{   real rTmp;
    rTmp=v0.x*v1.x+v0.y*v1.y+v0.z*v1.z;
    return rTmp;
}

vector operator*(real r0, vector v0)
{   vector vTmp;
    vTmp.x=r0*v0.x;
    vTmp.y=r0*v0.y;
    vTmp.z=r0*v0.z;
    return vTmp;
}

vector operator*(vector v0, real r0)
{   vector vTmp;
    vTmp.x=r0*v0.x;
    vTmp.y=r0*v0.y;
    vTmp.z=r0*v0.z;
    return vTmp;
}

cvector operator*(cplx c0, vector v0)
{   cvector cvTmp;
    cvTmp.x=c0*v0.x;
    cvTmp.y=c0*v0.y;
    cvTmp.z=c0*v0.z;
    return cvTmp;
}

cvector operator*(vector v0, cplx c0)
{   cvector cvTmp;
    cvTmp.x=c0*v0.x;
    cvTmp.y=c0*v0.y;
    cvTmp.z=c0*v0.z;
    return cvTmp;
}

cvector operator*(real r0, cvector cv0)
{   cvector cvTmp;
    cvTmp.x=r0*cv0.x;
    cvTmp.y=r0*cv0.y;
    cvTmp.z=r0*cv0.z;
    return cvTmp;
}

cvector operator*(cvector cv0, real r0)
{   cvector cvTmp;
    cvTmp.x=r0*cv0.x;
    cvTmp.y=r0*cv0.y;
    cvTmp.z=r0*cv0.z;
    return cvTmp;
}

cvector operator*(cplx c0, cvector cv0)
{   cvector cvTmp;
    cvTmp.x=c0*cv0.x;
    cvTmp.y=c0*cv0.y;
    cvTmp.z=c0*cv0.z;
    return cvTmp;
}

cvector operator*(cvector cv0, cplx c0)
{   cvector cvTmp;
    cvTmp.x=c0*cv0.x;
    cvTmp.y=c0*cv0.y;
    cvTmp.z=c0*cv0.z;
    return cvTmp;
}

cplx operator*(cvector cv0, cvector cv1)
{   cplx cTmp;
    cTmp=cv0.x*cv1.x+cv0.y*cv1.y+cv0.z*cv1.z;
    return cTmp;
}

cplx operator*(cvector cv0, vector v0)
{   cplx cTmp;
    cTmp=cv0.x*v0.x+cv0.y*v0.y+cv0.z*v0.z;
    return cTmp;
}

cplx operator*(vector v0, cvector cv0)
{   cplx cTmp;
    cTmp=cv0.x*v0.x+cv0.y*v0.y+cv0.z*v0.z;
    return cTmp;
}

vector operator^(vector v0, vector v1)
{   vector vTmp;
    vTmp.x=v0.y*v1.z-v0.z*v1.y;
    vTmp.y=v0.z*v1.x-v0.x*v1.z;
    vTmp.z=v0.x*v1.y-v0.y*v1.x;
    return vTmp;
}

cvector operator^(vector v0, cvector cv0)
{   cvector cvTmp;
    cvTmp.x=v0.y*cv0.z-v0.z*cv0.y;
    cvTmp.y=v0.z*cv0.x-v0.x*cv0.z;
    cvTmp.z=v0.x*cv0.y-v0.y*cv0.x;
    return cvTmp;
}

cvector operator^(cvector cv0, vector v0)
{   cvector cvTmp;
    cvTmp.x=v0.z*cv0.y-v0.y*cv0.z;
    cvTmp.y=v0.x*cv0.z-v0.z*cv0.x;
    cvTmp.z=v0.y*cv0.x-v0.x*cv0.y;
    return cvTmp;
}

cvector operator^(cvector cv0, cvector cv1)
{   cvector cvTmp;
    cvTmp.x=cv1.z*cv0.y-cv1.y*cv0.z;
    cvTmp.y=cv1.x*cv0.z-cv1.z*cv0.x;
    cvTmp.z=cv1.y*cv0.x-cv1.x*cv0.y;
    return cvTmp;
}

/************************* Commonplace functions ****************************/

cvector conj(cvector cv0)
{   cvector cvTmp;
    cvTmp.x=conj(cv0.x);
    cvTmp.y=conj(cv0.y);
    cvTmp.z=conj(cv0.z);
    return cvTmp;
}

vector re(cvector cv0)
{   vector vTmp;
    vTmp.x=re(cv0.x);
    vTmp.y=re(cv0.y);
    vTmp.z=re(cv0.z);
    return vTmp;
}

vector im(cvector cv0)
{   vector vTmp;
    vTmp.x=im(cv0.x);
    vTmp.y=im(cv0.y);
    vTmp.z=im(cv0.z);
    return vTmp;
}

real abs(vector v0)
{   real rTmp;
    rTmp=sqrt(v0.x*v0.x+v0.y*v0.y+v0.z*v0.z);
    return rTmp;
}

cplx abs(cvector cv0)
{   cplx cTmp;
    cTmp=sqrt(cv0.x*cv0.x+cv0.y*cv0.y+cv0.z*cv0.z);
    return cTmp;
}

/********************************* Matrix 3x3 *******************************/

matrix dyad(vector v0)
{   matrix mTmp;

	mTmp.e00=v0.x*v0.x;
	mTmp.e11=v0.y*v0.y;
	mTmp.e22=v0.z*v0.z;

	mTmp.e01=v0.x*v0.y;
	mTmp.e02=v0.x*v0.z;
	mTmp.e12=v0.y*v0.z;

	return mTmp;
}

cmatrix dyad(cvector cv0)
{   cmatrix cmTmp;

	cmTmp.e00=cv0.x*cv0.x;
	cmTmp.e11=cv0.y*cv0.y;
	cmTmp.e22=cv0.z*cv0.z;

	cmTmp.e01=cv0.x*cv0.y;
	cmTmp.e02=cv0.x*cv0.z;
	cmTmp.e12=cv0.y*cv0.z;

	return cmTmp;
}

void unit(matrix * m0)
{   matrix mTmp;

	mTmp.e00=1.;
	mTmp.e11=1.;
	mTmp.e22=1.;

	mTmp.e01=0.;
	mTmp.e02=0.;
	mTmp.e12=0.;

	* m0=mTmp;
}

void unit(cmatrix * cm0)
{   cmatrix cmTmp;

	cmTmp.e00.r=1.; cmTmp.e00.i=0.;
	cmTmp.e11.r=1.; cmTmp.e11.i=0.;
	cmTmp.e22.r=1.; cmTmp.e22.i=0.;

	cmTmp.e01.r=0.; cmTmp.e01.i=0.;
	cmTmp.e02.r=0.; cmTmp.e02.i=0.;
	cmTmp.e12.r=0.; cmTmp.e12.i=0.;

	* cm0=cmTmp;
}

/********************************* operator+ ********************************/

matrix operator+(matrix m0, matrix m1)
{   matrix mTmp;

	mTmp.e00=m0.e00+m1.e00;
	mTmp.e11=m0.e11+m1.e11;
	mTmp.e22=m0.e22+m1.e22;

	mTmp.e01=m0.e01+m1.e01;
	mTmp.e02=m0.e02+m1.e02;
	mTmp.e12=m0.e12+m1.e12;

	return mTmp;
}

cmatrix operator+(cmatrix cm0, matrix m1)
{   cmatrix cmTmp;

	cmTmp.e00=cm0.e00+m1.e00;
	cmTmp.e11=cm0.e11+m1.e11;
	cmTmp.e22=cm0.e22+m1.e22;

	cmTmp.e01=cm0.e01+m1.e01;
	cmTmp.e02=cm0.e02+m1.e02;
	cmTmp.e12=cm0.e12+m1.e12;

	return cmTmp;
}

cmatrix operator+(matrix m1, cmatrix cm0)
{   cmatrix cmTmp;

	cmTmp.e00=cm0.e00+m1.e00;
	cmTmp.e11=cm0.e11+m1.e11;
	cmTmp.e22=cm0.e22+m1.e22;

	cmTmp.e01=cm0.e01+m1.e01;
	cmTmp.e02=cm0.e02+m1.e02;
	cmTmp.e12=cm0.e12+m1.e12;

	return cmTmp;
}

cmatrix operator+(cmatrix cm0, cmatrix cm1)
{   cmatrix cmTmp;

	cmTmp.e00=cm0.e00+cm1.e00;
	cmTmp.e11=cm0.e11+cm1.e11;
	cmTmp.e22=cm0.e22+cm1.e22;

	cmTmp.e01=cm0.e01+cm1.e01;
	cmTmp.e02=cm0.e02+cm1.e02;
	cmTmp.e12=cm0.e12+cm1.e12;

	return cmTmp;
}

/********************************* operator- ********************************/

matrix operator-(matrix m0)
{   matrix mTmp;

	mTmp.e00=-m0.e00;
	mTmp.e11=-m0.e11;
	mTmp.e22=-m0.e22;

	mTmp.e01=-m0.e01;
	mTmp.e02=-m0.e02;
	mTmp.e12=-m0.e12;

	return mTmp;
}

matrix operator-(matrix m0, matrix m1)
{   matrix mTmp;

	mTmp.e00=m0.e00-m1.e00;
	mTmp.e11=m0.e11-m1.e11;
	mTmp.e22=m0.e22-m1.e22;

	mTmp.e01=m0.e01-m1.e01;
	mTmp.e02=m0.e02-m1.e02;
	mTmp.e12=m0.e12-m1.e12;

	return mTmp;
}

cmatrix operator-(cmatrix cm0)
{   cmatrix cmTmp;

	cmTmp.e00=-cm0.e00;
	cmTmp.e11=-cm0.e11;
	cmTmp.e22=-cm0.e22;

	cmTmp.e01=-cm0.e01;
	cmTmp.e02=-cm0.e02;
	cmTmp.e12=-cm0.e12;

	return cmTmp;
}

cmatrix operator-(matrix m0, cmatrix cm1)
{   cmatrix cmTmp;

	cmTmp.e00=m0.e00-cm1.e00;
	cmTmp.e11=m0.e11-cm1.e11;
	cmTmp.e22=m0.e22-cm1.e22;

	cmTmp.e01=m0.e01-cm1.e01;
	cmTmp.e02=m0.e02-cm1.e02;
	cmTmp.e12=m0.e12-cm1.e12;

	return cmTmp;
}

cmatrix operator-(cmatrix cm0, matrix m1)
{   cmatrix cmTmp;

	cmTmp.e00=cm0.e00-m1.e00;
	cmTmp.e11=cm0.e11-m1.e11;
	cmTmp.e22=cm0.e22-m1.e22;

	cmTmp.e01=cm0.e01-m1.e01;
	cmTmp.e02=cm0.e02-m1.e02;
	cmTmp.e12=cm0.e12-m1.e12;

	return cmTmp;
}

cmatrix operator-(cmatrix cm0, cmatrix cm1)
{   cmatrix cmTmp;

	cmTmp.e00=cm0.e00-cm1.e00;
	cmTmp.e11=cm0.e11-cm1.e11;
	cmTmp.e22=cm0.e22-cm1.e22;

	cmTmp.e01=cm0.e01-cm1.e01;
	cmTmp.e02=cm0.e02-cm1.e02;
	cmTmp.e12=cm0.e12-cm1.e12;

	return cmTmp;
}

/********************************* operator* ********************************/

matrix operator*(real r0, matrix m0)
{   matrix mTmp;

	mTmp.e00=r0*m0.e00;
	mTmp.e11=r0*m0.e11;
	mTmp.e22=r0*m0.e22;

	mTmp.e01=r0*m0.e01;
	mTmp.e02=r0*m0.e02;
	mTmp.e12=r0*m0.e12;

	return mTmp;
}

cmatrix operator*(cplx c0, matrix m0)
{   cmatrix cmTmp;

	cmTmp.e00=c0*m0.e00;
	cmTmp.e11=c0*m0.e11;
	cmTmp.e22=c0*m0.e22;

	cmTmp.e01=c0*m0.e01;
	cmTmp.e02=c0*m0.e02;
	cmTmp.e12=c0*m0.e12;

	return cmTmp;
}

cmatrix operator*(real r0, cmatrix cm0)
{   cmatrix cmTmp;

	cmTmp.e00=r0*cm0.e00;
	cmTmp.e11=r0*cm0.e11;
	cmTmp.e22=r0*cm0.e22;

	cmTmp.e01=r0*cm0.e01;
	cmTmp.e02=r0*cm0.e02;
	cmTmp.e12=r0*cm0.e12;

	return cmTmp;
}

cmatrix operator*(cplx c0, cmatrix cm0)
{   cmatrix cmTmp;

	cmTmp.e00=c0*cm0.e00;
	cmTmp.e11=c0*cm0.e11;
	cmTmp.e22=c0*cm0.e22;

	cmTmp.e01=c0*cm0.e01;
	cmTmp.e02=c0*cm0.e02;
	cmTmp.e12=c0*cm0.e12;

	return cmTmp;
}

vector operator*(matrix m0, vector v0)
{   vector vTmp;

	vTmp.x=m0.e00*v0.x+m0.e01*v0.y+m0.e02*v0.z;
	vTmp.y=m0.e01*v0.x+m0.e11*v0.y+m0.e12*v0.z;
	vTmp.z=m0.e02*v0.x+m0.e12*v0.y+m0.e22*v0.z;

	return vTmp;
}

cvector operator*(cmatrix cm0, vector v0)
{   cvector cvTmp;

	cvTmp.x=cm0.e00*v0.x+cm0.e01*v0.y+cm0.e02*v0.z;
	cvTmp.y=cm0.e01*v0.x+cm0.e11*v0.y+cm0.e12*v0.z;
	cvTmp.z=cm0.e02*v0.x+cm0.e12*v0.y+cm0.e22*v0.z;

	return cvTmp;
}

cvector operator*(matrix m0, cvector cv0)
{   cvector cvTmp;

	cvTmp.x=m0.e00*cv0.x+m0.e01*cv0.y+m0.e02*cv0.z;
	cvTmp.y=m0.e01*cv0.x+m0.e11*cv0.y+m0.e12*cv0.z;
	cvTmp.z=m0.e02*cv0.x+m0.e12*cv0.y+m0.e22*cv0.z;

	return cvTmp;
}

cvector operator*(cmatrix cm0, cvector cv0)
{   cvector cvTmp;

	cvTmp.x=cm0.e00*cv0.x+cm0.e01*cv0.y+cm0.e02*cv0.z;
	cvTmp.y=cm0.e01*cv0.x+cm0.e11*cv0.y+cm0.e12*cv0.z;
	cvTmp.z=cm0.e02*cv0.x+cm0.e12*cv0.y+cm0.e22*cv0.z;

	return cvTmp;
}

cmatrix conj(cmatrix cm0)
{   cmatrix cmTmp;

	cmTmp.e00.r=cm0.e00.r;  cmTmp.e00.i=-cm0.e00.i;
	cmTmp.e11.r=cm0.e11.r;  cmTmp.e11.i=-cm0.e11.i;
	cmTmp.e22.r=cm0.e22.r;  cmTmp.e22.i=-cm0.e22.i;

	cmTmp.e01.r=cm0.e01.r;  cmTmp.e01.i=-cm0.e01.i;
	cmTmp.e02.r=cm0.e02.r;  cmTmp.e02.i=-cm0.e02.i;
	cmTmp.e12.r=cm0.e12.r;  cmTmp.e12.i=-cm0.e12.i;

	return cmTmp;
}
