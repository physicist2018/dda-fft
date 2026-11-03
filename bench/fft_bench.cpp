/* Standalone harness for the 1-D FFT kernel.
 *
 * It is deliberately built outside the main binary (see bench/build.sh)
 * because `compile` globs every *.cpp in the project root.
 *
 * Modes:
 *   delta  <N>        - old kernel vs a double-precision DFT_N, all deltas
 *   half   <N>        - the zero-padding "half" composite vs a 2N-point DFT
 *   speed  <N> <reps> - ns per line for the old kernel
 *   roundtrip <N>     - FFT then iFFT, relative error
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include <random>
#include <chrono>

#include "complex.hpp"   /* real, cplx, PI */
#include "cplx.h"        /* Cmul and friends */
#include "cvector.hpp"   /* cvector        */
#include "fft_kernel.hpp" /* ddafft::plan, ddafft::fft_line */

/* The old kernel, as declared in function.hpp.  Including function.hpp
   would drag in the whole program, so the two prototypes are repeated. */

void FourierPhase(int i_max, cplx * cW);
void Get_number_bits(int number, int * number_of_bits);
void FFT(int lim_NoB, int i_max, cplx * cW, cvector * cvDATA_ini,
         cvector * cvDATA);
void iFFT(int lim_NoB, int i_max, cplx * cW, cvector * cvDATA_ini,
          cvector * cvDATA);

/****************************************************************************/
/******************************* Helpers ************************************/
/****************************************************************************/

typedef std::complex<double> dcmplx;

/* Reference DFT_N in double precision. */

static void ref_dft(int N, const std::vector<dcmplx> & in,
                    std::vector<dcmplx> & out)
{   out.assign((size_t)N, dcmplx(0.,0.));

    for(int k=0;k<N;k++)
    {   dcmplx s(0.,0.);
        for(int j=0;j<N;j++)
            s += in[(size_t)j] * std::polar(1., -2.*M_PI*(double)j*(double)k/(double)N);
        out[(size_t)k]=s;
    }
}

static double rel_err(dcmplx a, dcmplx b)
{   double d=std::abs(a-b);
    double s=std::abs(a)+std::abs(b);
    if(s==0.) return 0.;
    return d/s;
}

/* Twiddle table used by the kernel for a transform of size N:
   cW[k] = W_N^k = exp(-2*pi*i*k/N), k = 0 .. N/2-1.  This is exactly
   FourierPhase(N/2), as main.cpp builds it. */

static cplx * make_kernel_table(int N)
{   cplx * cW=(cplx *)malloc(((size_t)N/2)*sizeof(cplx));
    if(cW==NULL) { printf("out of memory\n"); exit(1); }
    FourierPhase(N/2,cW);
    return cW;
}

static int bits_for(int N)
{   int b=0;
    Get_number_bits(N-1,&b);
    return b;
}

/* The zero-padding composite of FFT_3D's inner line pass:
   a line of padded length 2N whose second half is zero.
   Forward bin 2l  = DFT_N(x)[l]
   Forward bin 2l+1 = DFT_N(x * W_{2N}^j)[l]
   Reproduced here exactly as fft_func.cpp does it. */

static void half_composite(int N, const std::vector<cvector> & x,
                           std::vector<cvector> & out)
{   int lim=bits_for(N);
    cplx * cW2=make_kernel_table(N);       /* W_N^k,   k < N/2 */
    cplx * cWf=(cplx *)malloc((size_t)N*sizeof(cplx));
    if(cWf==NULL) { printf("out of memory\n"); exit(1); }
    FourierPhase(N,cWf);                   /* W_{2N}^k, k < N   */

    std::vector<cvector> work((size_t)N), res((size_t)N), work2((size_t)N);

    /* even bins */
    for(int j=0;j<N;j++) work[(size_t)j]=x[(size_t)j];
    FFT(lim,N,cW2,work.data(),res.data());
    for(int j=0;j<N;j++) work2[(size_t)j]=res[(size_t)j];

    /* odd bins */
    for(int j=0;j<N;j++)
    {   cvector t=x[(size_t)j];
        Cmul(cWf[j],t.x,work[(size_t)j].x);
        Cmul(cWf[j],t.y,work[(size_t)j].y);
        Cmul(cWf[j],t.z,work[(size_t)j].z);
    }
    FFT(lim,N,cW2,work.data(),res.data());

    out.assign((size_t)(2*N),cvector{});
    for(int l=0;l<N;l++)
    {   out[(size_t)(2*l)]  =work2[(size_t)l];
        out[(size_t)(2*l+1)]=res[(size_t)l];
    }

    free(cW2); free(cWf);
}

/****************************************************************************/
/******************************** Modes *************************************/
/****************************************************************************/

static int mode_delta(int N)
{   int lim=bits_for(N);
    cplx * cW=make_kernel_table(N);

    std::vector<cvector> in((size_t)N), out((size_t)N);
    std::vector<dcmplx> ref((size_t)N), refin((size_t)N);

    double worst=0.;
    int worst_pos=-1, worst_k=-1;

    printf("delta test, N=%d, kernel table cW=FourierPhase(%d) -> W_%d^k\n",
           N,N/2,N);

    for(int pos=0;pos<N;pos++)
    {   for(int j=0;j<N;j++)
        {   in[(size_t)j].x.r=0.f; in[(size_t)j].x.i=0.f;
            in[(size_t)j].y.r=0.f; in[(size_t)j].y.i=0.f;
            in[(size_t)j].z.r=0.f; in[(size_t)j].z.i=0.f;
            refin[(size_t)j]=dcmplx(0.,0.);
        }
        in[(size_t)pos].x.r=1.f;
        refin[(size_t)pos]=dcmplx(1.,0.);

        FFT(lim,N,cW,in.data(),out.data());
        ref_dft(N,refin,ref);

        for(int k=0;k<N;k++)
        {   dcmplx got(out[(size_t)k].x.r,out[(size_t)k].x.i);
            double e=rel_err(ref[(size_t)k],got);
            if(e>worst) { worst=e; worst_pos=pos; worst_k=k; }
        }

        if(pos==1)
        {   printf("  delta at 1, first bins:");
            for(int k=0;k<3&&k<N;k++)
                printf("  %+.4f%+.4fi",out[(size_t)k].x.r,out[(size_t)k].x.i);
            printf("\n");
        }
    }

    printf("  max relative error = %.3e (pos %d, bin %d)\n",worst,worst_pos,worst_k);
    printf("  conclusion: kernel is %s DFT_%d\n",
           worst<1e-4?"a natural-order":"NOT a natural-order",N);

    free(cW);
    return worst<1e-4?0:1;
}

static int mode_half(int N)
{   std::vector<cvector> x((size_t)N), out;
    std::vector<dcmplx> refin((size_t)(2*N)), ref;

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> u(-1.,1.);

    for(int j=0;j<N;j++)
    {   x[(size_t)j].x.r=(float)u(rng); x[(size_t)j].x.i=(float)u(rng);
        x[(size_t)j].y.r=(float)u(rng); x[(size_t)j].y.i=(float)u(rng);
        x[(size_t)j].z.r=(float)u(rng); x[(size_t)j].z.i=(float)u(rng);
    }

    half_composite(N,x,out);

    for(int j=0;j<2*N;j++) refin[(size_t)j]=dcmplx(0.,0.);
    for(int j=0;j<N;j++)   refin[(size_t)j]=dcmplx(x[(size_t)j].x.r,x[(size_t)j].x.i);
    ref_dft(2*N,refin,ref);

    double worst=0.;
    for(int k=0;k<2*N;k++)
    {   dcmplx got(out[(size_t)k].x.r,out[(size_t)k].x.i);
        double e=rel_err(ref[(size_t)k],got);
        if(e>worst) worst=e;
    }

    printf("half composite, N=%d -> padded %d: max rel err vs DFT_%d = %.3e\n",
           N,2*N,2*N,worst);
    printf("  conclusion: composite is %s DFT_%d of the zero-padded line\n",
           worst<1e-4?"a full":"NOT a",2*N);

    return worst<1e-4?0:1;
}

static int mode_compare(int N, int reps)
{   int lim=bits_for(N);
    cplx * cW=make_kernel_table(N);
    const ddafft::plan * pl=ddafft::get_plan(N);
    if(pl==NULL) { printf("no plan for N=%d\n",N); return 1; }

    std::vector<cvector> pristine((size_t)N), work((size_t)N);
    std::vector<cvector> oldout((size_t)N), newout((size_t)N);

    std::mt19937 rng(2024);
    std::uniform_real_distribution<double> u(-1.,1.);
    for(int j=0;j<N;j++)
    {   pristine[(size_t)j].x.r=(float)u(rng); pristine[(size_t)j].x.i=(float)u(rng);
        pristine[(size_t)j].y.r=(float)u(rng); pristine[(size_t)j].y.i=(float)u(rng);
        pristine[(size_t)j].z.r=(float)u(rng); pristine[(size_t)j].z.i=(float)u(rng);
    }

    work=pristine;
    FFT(lim,N,cW,work.data(),oldout.data());
    ddafft::fft_line(pl,newout.data(),pristine.data());

    double num=0., den=0.;
    for(int j=0;j<N;j++)
    {   const cvector &a=oldout[(size_t)j], &b=newout[(size_t)j];
        dcmplx oa(a.x.r,a.x.i), ob(b.x.r,b.x.i);
        num+=std::norm(oa-ob); den+=std::norm(oa);
        oa=dcmplx(a.y.r,a.y.i); ob=dcmplx(b.y.r,b.y.i);
        num+=std::norm(oa-ob); den+=std::norm(oa);
        oa=dcmplx(a.z.r,a.z.i); ob=dcmplx(b.z.r,b.z.i);
        num+=std::norm(oa-ob); den+=std::norm(oa);
    }
    double err=std::sqrt(num/den);

    /* timing: identical harness for both (refill + transform) */
    for(int r=0;r<10;r++)
    {   memcpy(work.data(),pristine.data(),(size_t)N*sizeof(cvector));
        FFT(lim,N,cW,work.data(),oldout.data());
        ddafft::fft_line(pl,newout.data(),work.data());
    }

    volatile float sink=0.f;

    auto t0=std::chrono::steady_clock::now();
    for(int r=0;r<reps;r++)
    {   memcpy(work.data(),pristine.data(),(size_t)N*sizeof(cvector));
        FFT(lim,N,cW,work.data(),oldout.data());
        sink+=oldout[(size_t)(r%N)].x.r;
    }
    auto t1=std::chrono::steady_clock::now();

    auto t2=std::chrono::steady_clock::now();
    for(int r=0;r<reps;r++)
    {   memcpy(work.data(),pristine.data(),(size_t)N*sizeof(cvector));
        ddafft::fft_line(pl,newout.data(),work.data());
        sink+=newout[(size_t)(r%N)].x.r;
    }
    auto t3=std::chrono::steady_clock::now();

    double ns_old=std::chrono::duration<double,std::nano>(t1-t0).count()/(double)reps;
    double ns_new=std::chrono::duration<double,std::nano>(t3-t2).count()/(double)reps;

    printf("compare N=%d: rel L2 old-vs-new = %.3e | old %.1f ns, new %.1f ns, speedup %.2fx\n",
           N,err,ns_old,ns_new,ns_old/ns_new);

    free(cW);
    return err<1e-5?0:1;
}

static int mode_newdelta(int N)
{   const ddafft::plan * pl=ddafft::get_plan(N);
    if(pl==NULL) { printf("no plan for N=%d\n",N); return 1; }

    std::vector<cvector> in((size_t)N), out((size_t)N);
    std::vector<dcmplx> ref((size_t)N), refin((size_t)N);

    double worst=0.;
    int wp=-1,wk=-1;

    for(int pos=0;pos<N;pos++)
    {   for(int j=0;j<N;j++)
        {   in[(size_t)j]=cvector{};
            refin[(size_t)j]=dcmplx(0.,0.);
        }
        in[(size_t)pos].x.r=1.f;
        refin[(size_t)pos]=dcmplx(1.,0.);

        ddafft::fft_line(pl,out.data(),in.data());
        ref_dft(N,refin,ref);

        for(int k=0;k<N;k++)
        {   dcmplx got(out[(size_t)k].x.r,out[(size_t)k].x.i);
            double e=rel_err(ref[(size_t)k],got);
            if(e>worst) { worst=e; wp=pos; wk=k; }
        }
    }

    printf("new delta: N=%d max rel err = %.3e (pos %d bin %d) -> %s DFT_%d\n",
           N,worst,wp,wk,worst<1e-4?"natural-order":"NOT a natural-order",N);
    return worst<1e-4?0:1;
}

static int mode_accuracy(int N)
{   int lim=bits_for(N);
    cplx * cW=make_kernel_table(N);
    const ddafft::plan * pl=ddafft::get_plan(N);
    if(pl==NULL) { printf("no plan for N=%d\n",N); return 1; }

    std::vector<cvector> in((size_t)N), work((size_t)N), oldout((size_t)N), newout((size_t)N);
    std::vector<dcmplx> refin((size_t)N), ref((size_t)N);

    std::mt19937 rng(31415);
    std::uniform_real_distribution<double> u(-1.,1.);
    for(int j=0;j<N;j++)
    {   in[(size_t)j].x.r=(float)u(rng); in[(size_t)j].x.i=(float)u(rng);
        in[(size_t)j].y.r=(float)u(rng); in[(size_t)j].y.i=(float)u(rng);
        in[(size_t)j].z.r=(float)u(rng); in[(size_t)j].z.i=(float)u(rng);
        refin[(size_t)j]=dcmplx(in[(size_t)j].x.r,in[(size_t)j].x.i);
    }

    work=in;
    FFT(lim,N,cW,work.data(),oldout.data());
    ddafft::fft_line(pl,newout.data(),in.data());
    ref_dft(N,refin,ref);

    double eo=0., en=0., den=0.;
    for(int k=0;k<N;k++)
    {   dcmplx r=ref[(size_t)k];
        eo += std::norm(dcmplx(oldout[(size_t)k].x.r,oldout[(size_t)k].x.i)-r);
        en += std::norm(dcmplx(newout[(size_t)k].x.r,newout[(size_t)k].x.i)-r);
        den+= std::norm(r);
    }

    printf("accuracy N=%d vs double DFT: old %.3e, new %.3e (relative L2)\n",
           N,std::sqrt(eo/den),std::sqrt(en/den));

    free(cW);
    return 0;
}

static int mode_speed(int N, int reps)
{   int lim=bits_for(N);
    cplx * cW=make_kernel_table(N);

    std::vector<cvector> pristine((size_t)N), work((size_t)N), out((size_t)N);

    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(-1.,1.);
    for(int j=0;j<N;j++)
    {   pristine[(size_t)j].x.r=(float)u(rng); pristine[(size_t)j].x.i=(float)u(rng);
        pristine[(size_t)j].y.r=(float)u(rng); pristine[(size_t)j].y.i=(float)u(rng);
        pristine[(size_t)j].z.r=(float)u(rng); pristine[(size_t)j].z.i=(float)u(rng);
    }

    /* warm-up */
    for(int r=0;r<10;r++)
    {   memcpy(work.data(),pristine.data(),(size_t)N*sizeof(cvector));
        FFT(lim,N,cW,work.data(),out.data());
    }

    volatile float sink=out[1].x.r;

    auto t0=std::chrono::steady_clock::now();
    for(int r=0;r<reps;r++)
    {   memcpy(work.data(),pristine.data(),(size_t)N*sizeof(cvector));
        FFT(lim,N,cW,work.data(),out.data());
        sink+=out[(size_t)(r%N)].y.i;
    }
    auto t1=std::chrono::steady_clock::now();

    double ns=std::chrono::duration<double,std::nano>(t1-t0).count();

    /* price of the memcpy alone, so it can be subtracted */
    auto s0=std::chrono::steady_clock::now();
    for(int r=0;r<reps;r++)
    {   memcpy(work.data(),pristine.data(),(size_t)N*sizeof(cvector));
        sink+=work[(size_t)(r%N)].z.i;
    }
    auto s1=std::chrono::steady_clock::now();
    double ns_copy=std::chrono::duration<double,std::nano>(s1-s0).count();

    printf("speed old kernel: N=%d reps=%d  %.1f ns/line (memcpy %.1f) -> %.1f ns/line net\n",
           N,reps,ns/(double)reps,ns_copy/(double)reps,
           (ns-ns_copy)/(double)reps);

    free(cW);
    return 0;
}

static int mode_roundtrip(int N)
{   int lim=bits_for(N);
    cplx * cW=make_kernel_table(N);

    std::vector<cvector> orig((size_t)N), in((size_t)N), out((size_t)N), back((size_t)N);

    std::mt19937 rng(99);
    std::uniform_real_distribution<double> u(-1.,1.);
    for(int j=0;j<N;j++)
    {   orig[(size_t)j].x.r=(float)u(rng); orig[(size_t)j].x.i=(float)u(rng);
        orig[(size_t)j].y.r=(float)u(rng); orig[(size_t)j].y.i=(float)u(rng);
        orig[(size_t)j].z.r=(float)u(rng); orig[(size_t)j].z.i=(float)u(rng);
    }

    /* the old kernel clobbers its input buffer, so hand it a copy */
    in=orig;
    FFT(lim,N,cW,in.data(),out.data());
    iFFT(lim,N,cW,out.data(),back.data());

    double num=0., den=0.;
    for(int j=0;j<N;j++)
    {   /* iFFT here is conj/FFT/conj without 1/N, so back = N * orig */
        dcmplx want((double)orig[(size_t)j].x.r*(double)N,(double)orig[(size_t)j].x.i*(double)N);
        num += std::norm(dcmplx((double)(back[(size_t)j].x.r),(double)(back[(size_t)j].x.i)) - want);
        want=dcmplx((double)orig[(size_t)j].y.r*(double)N,(double)orig[(size_t)j].y.i*(double)N);
        num += std::norm(dcmplx((double)(back[(size_t)j].y.r),(double)(back[(size_t)j].y.i)) - want);
        want=dcmplx((double)orig[(size_t)j].z.r*(double)N,(double)orig[(size_t)j].z.i*(double)N);
        num += std::norm(dcmplx((double)(back[(size_t)j].z.r),(double)(back[(size_t)j].z.i)) - want);
        den += std::norm(dcmplx((double)orig[(size_t)j].x.r*(double)N,(double)orig[(size_t)j].x.i*(double)N));
        den += std::norm(dcmplx((double)orig[(size_t)j].y.r*(double)N,(double)orig[(size_t)j].y.i*(double)N));
        den += std::norm(dcmplx((double)orig[(size_t)j].z.r*(double)N,(double)orig[(size_t)j].z.i*(double)N));
    }

    printf("roundtrip N=%d: relative L2 error vs N*x = %.3e\n",N,std::sqrt(num/den));
    free(cW);
    return 0;
}

/****************************************************************************/

int main(int argc, char ** argv)
{   if(argc<3)
    {   printf("usage: %s delta|newdelta|accuracy|compare|half|speed|roundtrip N [reps]\n",argv[0]);
        return 2;
    }

    const char * mode=argv[1];
    int N=atoi(argv[2]);

    if(N<2 || (N&(N-1))!=0)
    {   printf("N must be a power of two and at least 2\n");
        return 2;
    }

    /* The old kernel is driven by FourierPhase(N/2).  For N = 2 that is
       FourierPhase(1), which writes cW[1] past the end of a one-element
       table -- a pre-existing bug in fft_func.cpp, reachable only for a
       lattice dimension of 2.  The new kernel is unaffected, so only the
       new-kernel mode is allowed below N = 4. */
    if(N<4 && strcmp(mode,"newdelta")!=0)
    {   printf("modes that exercise the old kernel need N >= 4 "
               "(FourierPhase(1) overflows); use newdelta for N=2\n");
        return 2;
    }

    if(!strcmp(mode,"delta"))     return mode_delta(N);
    if(!strcmp(mode,"newdelta"))  return mode_newdelta(N);
    if(!strcmp(mode,"accuracy"))  return mode_accuracy(N);
    if(!strcmp(mode,"compare"))
    {   int reps=argc>3?atoi(argv[3]):20000;
        return mode_compare(N,reps);
    }
    if(!strcmp(mode,"half"))      return mode_half(N);
    if(!strcmp(mode,"roundtrip")) return mode_roundtrip(N);
    if(!strcmp(mode,"speed"))
    {   int reps=argc>3?atoi(argv[3]):20000;
        return mode_speed(N,reps);
    }

    printf("unknown mode %s\n",mode);
    return 2;
}
