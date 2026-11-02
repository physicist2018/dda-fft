
#ifndef _VECTOR_ESZ_DEFINED_

#include "complex.hpp"
#define _VECTOR_ESZ_DEFINED_

typedef
	struct
	{   real x;   real y;   real z;
	}   vector;

typedef
	struct
	{   cplx x;     cplx y;     cplx z;
	}   cvector;

typedef
	struct
	{   real e00;   real e01;   real e02;
		real e11;   real e12;   real e22;
	}   matrix;

typedef
	struct
	{   cplx e00;   cplx e01;   cplx e02;
		cplx e11;   cplx e12;   cplx e22;
	}   cmatrix;

vector  operator+(vector v0, vector v1);
cvector operator+(cvector cv0, vector v0);
cvector operator+(vector v0, cvector cv0);
cvector operator+(cvector cv0, cvector cv1);
vector  operator-(vector v0, vector v1);
cvector operator-(cvector cv0, vector v0);
cvector operator-(vector v0, cvector cv0);
cvector operator-(cvector cv0, cvector cv1);
vector  operator-(vector v0);
cvector operator-(cvector cv0);
real    operator*(vector v0, vector v1);
vector  operator*(real r0, vector v0);
vector  operator*(vector v0, real r0);
cvector operator*(cplx c0, vector v0);
cvector operator*(vector v0, cplx c0);
cvector operator*(real r0, cvector cv0);
cvector operator*(cvector cv0, real r0);
cvector operator*(cplx c0, cvector cv0);
cvector operator*(cvector cv0, cplx c0);
cplx    operator*(cvector cv0, cvector cv1);
cplx    operator*(cvector cv0, vector v0);
cplx    operator*(vector v0, cvector cv0);
vector  operator^(vector v0, vector v1);
cvector operator^(vector v0, cvector cv0);
cvector operator^(cvector cv0, vector v0);
cvector operator^(cvector cv0, cvector cv1);
cvector conj(cvector cv0);
vector  re(cvector cv0);
vector  im(cvector cv0);
real    abs(vector v0);
cplx    abs(cvector cv0);

/********************************* Matrix 3x3 *******************************/

matrix  dyad(vector v0);
cmatrix dyad(cvector cv0);
void unit(matrix * m0);
void unit(cmatrix * cm0);
cmatrix conj(cmatrix cm0);

matrix operator+(matrix m0, matrix m1);
cmatrix operator+(cmatrix cm0, matrix m1);
cmatrix operator+(matrix m1, cmatrix cm0);
cmatrix operator+(cmatrix cm0, cmatrix cm1);

matrix operator-(matrix m0);
matrix operator-(matrix m0, matrix m1);
cmatrix operator-(cmatrix cm0);
cmatrix operator-(matrix m0, cmatrix cm1);
cmatrix operator-(cmatrix cm0, matrix m1);
cmatrix operator-(cmatrix cm0, cmatrix cm1);

matrix operator*(real r0, matrix m0);
cmatrix operator*(cplx c0, matrix m0);
cmatrix operator*(real r0, cmatrix cm0);
cmatrix operator*(cplx c0, cmatrix cm0);

vector operator*(matrix m0, vector v0);
cvector operator*(cmatrix cm0, vector v0);
cvector operator*(matrix m0, cvector cv0);
cvector operator*(cmatrix cm0, cvector cv0);

#endif
