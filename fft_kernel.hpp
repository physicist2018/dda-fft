/* Mixed-radix (4/8) in-place 1-D FFT kernel for the DDA convolution.
 *
 * The old kernel in fft_func.cpp is a radix-2 decimation-in-time transform
 * that (a) recomputes the bit reversal on every call and (b) copies the whole
 * buffer between two ping-pong arrays on every stage.  For a line of N points
 * that is log2(N) full copies per transform, which dominates at N = 64..128.
 *
 * This kernel keeps exactly the same mathematical contract -- OUT[k] =
 * sum_j IN[j] * exp(-2*pi*i*j*k/N), natural order -- so the "half" convention
 * of the zero-padded convolution (see FFT_3D: the padded line of length 2N has
 * only its first N samples non-zero, and the driver combines two N-point
 * transforms into even/odd bins) is preserved unchanged.
 *
 * What is different:
 *   - the permutation (mixed-radix digit reversal) is precomputed once per
 *     size, together with a full table of W_N^k;
 *   - the transform runs in place, so there are no per-stage buffer copies;
 *   - radix-4 / radix-8 butterflies with the usual trivial-twiddle paths.
 *
 * The input buffer is read only; OUT may not alias IN.
 */

#ifndef DDA_FFT_KERNEL_HPP
#define DDA_FFT_KERNEL_HPP

#include <cmath>
#include <cstdlib>
#include <mutex>

#include "complex.hpp"
#include "cvector.hpp"

namespace ddafft {

struct plan
{   int N;
    int nstages;
    int radix[24];
    int * perm;     /* perm[k]: source index of output slot k */
    float * twr;    /* W_N^k = twr[k] + i*twi[k], k = 0 .. N-1 */
    float * twi;
};

/****************************************************************************/
/************************** small complex helpers ***************************/
/****************************************************************************/

static inline void bf2c(cplx & a, cplx & b)
{   float tr=a.r+b.r, ti=a.i+b.i;
    b.r=a.r-b.r;  b.i=a.i-b.i;
    a.r=tr;       a.i=ti;
}

/* 4-point DFT: (x0,x1,x2,x3) -> (X0,X1,X2,X3), W_4 = -i */

static inline void bf4c(cplx & x0, cplx & x1, cplx & x2, cplx & x3)
{   float t0r=x0.r+x2.r, t0i=x0.i+x2.i;
    float t1r=x0.r-x2.r, t1i=x0.i-x2.i;
    float t2r=x1.r+x3.r, t2i=x1.i+x3.i;
    float t3r=x1.r-x3.r, t3i=x1.i-x3.i;

    x0.r=t0r+t2r; x0.i=t0i+t2i;
    x2.r=t0r-t2r; x2.i=t0i-t2i;
    x1.r=t1r+t3i; x1.i=t1i-t3r;      /* t1 - i t3 */
    x3.r=t1r-t3i; x3.i=t1i+t3r;      /* t1 + i t3 */
}

/* 8-point DFT built from two 4-point transforms (Burrus, radix-8 core). */

static inline void bf8c(cplx & z0, cplx & z1, cplx & z2, cplx & z3,
                        cplx & z4, cplx & z5, cplx & z6, cplx & z7)
{   const float s=0.70710678118654752440f;

    /* even bins: DFT_4 of A_d = z_d + z_{d+4} */
    cplx a0,a1,a2,a3,b0,b1,b2,b3;
    a0.r=z0.r+z4.r; a0.i=z0.i+z4.i;
    a1.r=z1.r+z5.r; a1.i=z1.i+z5.i;
    a2.r=z2.r+z6.r; a2.i=z2.i+z6.i;
    a3.r=z3.r+z7.r; a3.i=z3.i+z7.i;

    /* odd bins: DFT_4 of B_d = (z_d - z_{d+4}) * W_8^d */
    b0.r=z0.r-z4.r; b0.i=z0.i-z4.i;

    {   float dr=z1.r-z5.r, di=z1.i-z5.i;          /* W_8^1 */
        b1.r=(dr+di)*s; b1.i=(di-dr)*s;
    }
    {   float dr=z2.r-z6.r, di=z2.i-z6.i;          /* W_8^2 = -i */
        b2.r=di; b2.i=-dr;
    }
    {   float dr=z3.r-z7.r, di=z3.i-z7.i;          /* W_8^3 = -(1+i)/sqrt2 */
        b3.r=(di-dr)*s; b3.i=-(dr+di)*s;
    }

    bf4c(a0,a1,a2,a3);
    bf4c(b0,b1,b2,b3);

    /* a_p -> bin 2p, b_p -> bin 2p+1 */
    z0=a0; z1=b0; z2=a1; z3=b1;
    z4=a2; z5=b2; z6=a3; z7=b3;
}

static inline void cmulv(cvector & v, float wr, float wi)
{   float tr;
    tr=v.x.r*wr-v.x.i*wi; v.x.i=v.x.r*wi+v.x.i*wr; v.x.r=tr;
    tr=v.y.r*wr-v.y.i*wi; v.y.i=v.y.r*wi+v.y.i*wr; v.y.r=tr;
    tr=v.z.r*wr-v.z.i*wi; v.z.i=v.z.r*wi+v.z.i*wr; v.z.r=tr;
}

/****************************************************************************/
/****************************** stage kernels *******************************/
/****************************************************************************/

static inline void stage2(cvector * v, int N, int span, int h, int step,
                          const float * twr, const float * twi)
{   (void)step; (void)twr; (void)twi;

    for(int base=0;base<N;base+=span)
        for(int i=0;i<h;i++)
        {   cvector & a=v[base+i], & b=v[base+i+h];
            bf2c(a.x,b.x); bf2c(a.y,b.y); bf2c(a.z,b.z);
        }
}

static inline void stage4(cvector * v, int N, int span, int h, int step,
                          const float * twr, const float * twi)
{   for(int base=0;base<N;base+=span)
    {   /* i == 0: all twiddles are 1 */
        if(h>0)
        {   cvector & a0=v[base], & a1=v[base+h];
            cvector & a2=v[base+2*h], & a3=v[base+3*h];
            bf4c(a0.x,a1.x,a2.x,a3.x);
            bf4c(a0.y,a1.y,a2.y,a3.y);
            bf4c(a0.z,a1.z,a2.z,a3.z);
        }

        for(int i=1;i<h;i++)
        {   cvector & a0=v[base+i];
            cvector & a1=v[base+i+h];
            cvector & a2=v[base+i+2*h];
            cvector & a3=v[base+i+3*h];

            cmulv(a1,twr[i*step],twi[i*step]);
            cmulv(a2,twr[2*i*step],twi[2*i*step]);
            cmulv(a3,twr[3*i*step],twi[3*i*step]);

            bf4c(a0.x,a1.x,a2.x,a3.x);
            bf4c(a0.y,a1.y,a2.y,a3.y);
            bf4c(a0.z,a1.z,a2.z,a3.z);
        }
    }
}

static inline void stage8(cvector * v, int N, int span, int h, int step,
                          const float * twr, const float * twi)
{   for(int base=0;base<N;base+=span)
    {   /* i == 0: all twiddles are 1 */
        {   cvector & a0=v[base];
            cvector & a1=v[base+h],   & a2=v[base+2*h], & a3=v[base+3*h];
            cvector & a4=v[base+4*h], & a5=v[base+5*h], & a6=v[base+6*h];
            cvector & a7=v[base+7*h];
            bf8c(a0.x,a1.x,a2.x,a3.x,a4.x,a5.x,a6.x,a7.x);
            bf8c(a0.y,a1.y,a2.y,a3.y,a4.y,a5.y,a6.y,a7.y);
            bf8c(a0.z,a1.z,a2.z,a3.z,a4.z,a5.z,a6.z,a7.z);
        }

        for(int i=1;i<h;i++)
        {   cvector & a0=v[base+i];
            cvector & a1=v[base+i+h],   & a2=v[base+i+2*h], & a3=v[base+i+3*h];
            cvector & a4=v[base+i+4*h], & a5=v[base+i+5*h], & a6=v[base+i+6*h];
            cvector & a7=v[base+i+7*h];

            cmulv(a1,twr[i*step],twi[i*step]);
            cmulv(a2,twr[2*i*step],twi[2*i*step]);
            cmulv(a3,twr[3*i*step],twi[3*i*step]);
            cmulv(a4,twr[4*i*step],twi[4*i*step]);
            cmulv(a5,twr[5*i*step],twi[5*i*step]);
            cmulv(a6,twr[6*i*step],twi[6*i*step]);
            cmulv(a7,twr[7*i*step],twi[7*i*step]);

            bf8c(a0.x,a1.x,a2.x,a3.x,a4.x,a5.x,a6.x,a7.x);
            bf8c(a0.y,a1.y,a2.y,a3.y,a4.y,a5.y,a6.y,a7.y);
            bf8c(a0.z,a1.z,a2.z,a3.z,a4.z,a5.z,a6.z,a7.z);
        }
    }
}

/****************************************************************************/
/**************************** plan construction *****************************/
/****************************************************************************/

/* log2(N) = 3a + b; the large radix goes into the early (contiguous) stages. */

static inline int choose_radices(int k, int * r)
{   int m=0;

    while(k>6) { r[m++]=8; k-=3; }

    switch(k)
    {   case 0: break;
        case 1: r[m++]=2; break;
        case 2: r[m++]=4; break;
        case 3: r[m++]=8; break;
        case 4: r[m++]=4; r[m++]=4; break;
        case 5: r[m++]=8; r[m++]=4; break;
        case 6: r[m++]=8; r[m++]=8; break;
    }
    return m;
}

static inline bool build_plan(plan & p, int N)
{   int k=0;
    while((1<<k)<N) k++;
    if((1<<k)!=N) return false;

    p.N=N;
    p.nstages=choose_radices(k,p.radix);

    p.perm=(int *)malloc((size_t)N*sizeof(int));
    p.twr=(float *)malloc((size_t)N*sizeof(float));
    p.twi=(float *)malloc((size_t)N*sizeof(float));
    if(p.perm==NULL||p.twr==NULL||p.twi==NULL) return false;

    /* mixed-radix digit reversal: digits are taken in the stage order, and
       the first digit becomes the most significant part of the result. */

    for(int j=0;j<N;j++)
    {   int tmp=j, L=1, val=0;
        for(int t=0;t<p.nstages;t++)
        {   int rt=p.radix[t];
            int d=tmp%rt;
            tmp/=rt;
            L*=rt;
            val+=d*(N/L);
        }
        p.perm[j]=val;
    }

    const double two_pi=6.283185307179586476925286766559;
    for(int q=0;q<N;q++)
    {   double a=-two_pi*(double)q/(double)N;
        p.twr[q]=(float)cos(a);
        p.twi[q]=(float)sin(a);
    }

    return true;
}

/* One plan per size, built once and then read-only, like the cW tables. */

static inline const plan * get_plan(int N)
{   static plan cache[16];
    static int ncache=0;
    static std::mutex m;

    if(N<1 || (N&(N-1))!=0) return NULL;

    std::lock_guard<std::mutex> g(m);

    for(int i=0;i<ncache;i++)
        if(cache[i].N==N) return &cache[i];

    if(ncache>=16) return NULL;
    if(!build_plan(cache[ncache],N)) return NULL;

    return &cache[ncache++];
}

/****************************************************************************/
/******************************** the kernel ********************************/
/****************************************************************************/

/* out[k] = DFT_N(in)[k], natural order.  in is not modified. */

static inline void run_stages(const plan * p, cvector * v)
{   const int N=p->N;

    int span=1;
    for(int t=0;t<p->nstages;t++)
    {   int r=p->radix[t];
        span*=r;

        int h=span/r;
        int step=N/span;

        if(r==8)      stage8(v,N,span,h,step,p->twr,p->twi);
        else if(r==4) stage4(v,N,span,h,step,p->twr,p->twi);
        else          stage2(v,N,span,h,step,p->twr,p->twi);
    }
}

static inline void fft_line(const plan * p, cvector * out, const cvector * in)
{   const int N=p->N;
    const int * perm=p->perm;

    for(int k=0;k<N;k++)
        out[k]=in[perm[k]];

    run_stages(p,out);
}

/* Inverse transform in the same "no 1/N" convention as iFFT() in
   fft_func.cpp: out = conj( DFT( conj(in) ) ). */

static inline void ifft_line(const plan * p, cvector * out, const cvector * in)
{   const int N=p->N;
    const int * perm=p->perm;

    for(int k=0;k<N;k++)
    {   cvector t=in[perm[k]];
        out[k].x.r=t.x.r; out[k].x.i=-t.x.i;
        out[k].y.r=t.y.r; out[k].y.i=-t.y.i;
        out[k].z.r=t.z.r; out[k].z.i=-t.z.i;
    }

    run_stages(p,out);

    for(int k=0;k<N;k++)
    {   out[k].x.i=-out[k].x.i;
        out[k].y.i=-out[k].y.i;
        out[k].z.i=-out[k].z.i;
    }
}

} /* namespace ddafft */

#endif
