
#ifndef _CPLX_MOD_DEFINED_

#include <math.h>
#define _CPLX_MOD_DEFINED_

/************************ Complex operation *********************************/

/*  void Ctopa(cplx x,* cplx c)  å=(* c)  */
#define  Ctopa(x,c)  (c)->r=(x).r,(c)->i=(x).i

/*  void Cfropa(* cplx x, cplx c)   (* c)=x */
#define  Cfropa(x,c)  (c).r=(x)->r,(c).i=(x)->i

/*  void Ccopy(cplx x, cplx c)   c=x  */
#define  Ccopy(x,c)  (c).r=(x).r,(c).i=(x).i

/* real  CIm(cplx x) { return (  Im(x) ); }  Im */
#define  CIm(x) ((x).i)

/* real  CRe(cplx x) { return (  Re(x) ); }  Re */
#define  CRe(x) ((x).r)

/* void  Ccplx(real x, real y, cplx c) { c = (x+y*i); } */
#define  Ccplx(x,y,c)  (c).r=(x),(c).i=(y)

/* real  Cabs(cplx x) { return ( fabs(x) ); }  |å|  */
#define  Cabs(x) (hypot((x).r,(x).i))

/* void  Cneg(cplx x) { x = -x; }    */
#define  Cneg(x)  (x).r=-(x).r,(x).i=-(x).i

/* void  Cadd(cplx x, cplx y, cplx c) { c = x + y; }  */
#define  Cadd(x,y,c)  (c).r=(x).r+(y).r,(c).i=(x).i+(y).i

/* void  Csub(cplx x, cplx y, cplx c) { c = x - y; }  */
#define  Csub(x,y,c)  (c).r=(x).r-(y).r,(c).i=(x).i-(y).i

/* void  Cmul(cplx x, cplx y, cplx c) { c = x * y; }  */
#define  Cmul(x,y,c)  (c).r=(x).r*(y).r - (x).i*(y).i,\
							 (c).i=(x).r*(y).i + (x).i*(y).r

/* void  Cmuls(cplx x, cplx y, cplx c) { c = x * y(sopr); } */
#define  Cmuls(x,y,c) (c).r=(x).r*(y).r + (x).i*(y).i,\
							 (c).i=(x).i*(y).r - (x).r*(y).i

/* void  Cdiv(cplx x, cplx y, cplx c) { c = x / y; }  */
#define  Cdiv(x,y,c) (c).i=(y).r*(y).r + (y).i*(y).i,        \
							(c).r=((x).r*(y).r + (x).i*(y).i)/(c).i,\
							(c).i=(-(x).r*(y).i + (x).i*(y).r)/(c).i

/* void  Cexp(cplx x, cplx c) { c = exp(x); } */
#define  Cexp(x,c)   (c).i=exp((x).r),      \
							(c).r=cos((x).i)*(c).i,\
							(c).i=sin((x).i)*(c).i

/* void  Csqrt(cplx x, cplx c) { c = sqrt(x); }   */
#define  Csqrt(x,c)  (fabs((x).i)<1e-7 && fabs((x).r)<1e-7 ) ?\
							((c).r=0.,(r).i=0.) :                    \
							((c).i=atan2((x).i,(x).r)/2.,            \
							 (c).r=(cos((c).i))*sqrt(Cabs(x)),       \
							 (c).i=(sin((c).i))*sqrt(Cabs(x)))

/*********************** real_Complex operation ***************************/

/* void  Cexp_im(real x, cplx c) { c = exp(x); }        */
#define  Cexp_im(x,c) (c).r=cos((x)) , (c).i=sin((x))

/* void  DCadd(real x, cplx y, cplx c) { c = x + y; }   */
#define  DCadd(x,y,c)  (c).r=(x)+(y).r , (c).i=(y).i

/* void  DCsub(real x, cplx y, cplx c) { c = x - y; }   */
#define  DCsub(x,y,c)  (c).r=(x)-(y).r , (c).i=-(y).i

/* void  CDsub(cplx x, real y, cplx c) { c = x - y; }   */
#define  CDsub(x,y,c)  (c).r=(x).r-(y) , (c).i=(x).i

/* void  DCmul(real x, cplx y, cplx c) { c = x * y; }   */
#define  DCmul(x,y,c)  (c).r=(x)*(y).r , (c).i=(x)*(y).i

/* void  DCdiv(real x, cplx y, cplx c) { c = x / y; }   */
#define  DCdiv(x,y,c)  (c).i=(y).r*(y).r + (y).i*(y).i,\
							  (c).r=(x)*(y).r/(c).i , (c).i=-(x)*(y).i/(c).i

/* void  CDdiv(cplx x, real y, cplx c) { c = x / y; }   */
#define  CDdiv(x,y,c)  (c).r=(x).r/y, (c).i=(x).i/y

/****************************************************************************/

#endif
