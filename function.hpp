
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>

#include "cplx.h"
#include "cvector.hpp"
#include "mmuller.hpp"
#include "randomize.hpp"

typedef
    struct
    {   real r;     // radius of cylinder
        real l;     // length of cylinder
        vector b;   // begin point of cylinder
        vector e;   // end point of cylinder
        vector n;   // axis of cylinder |n|=1
    }
    cylinder;

/* Converting spherical coordinate system to Cartesian coordinate system */
/* angles are measured in degrees */
#define SphCart_d( r, phi, theta, result)   (result).x=r*cos(phi*PI/180.)*sin(theta*PI/180.),\
                                            (result).y=r*sin(phi*PI/180.)*sin(theta*PI/180.),\
                                            (result).z=r*cos(theta*PI/180.)

/* angles are measured in radians */
#define SphCart_r( r, phi, theta, result)   (result).x=r*cos(phi)*sin(theta),\
                                            (result).y=r*sin(phi)*sin(theta),\
                                            (result).z=r*cos(theta)
#define ZERO 0.0000000001

/*********************************** FFT ************************************/

void Change_binary_number(int number_old, int number_of_bits, int * number_new);
void Get_number_bits(int number, int * number_of_bits);
void FourierPhase(int i_max, cplx * cW);
void FourierPhase(int i_max, real ** Im_cW, cplx ** cW2);


void FFT(int lim_NoB, int i_max, cplx * cW, cmatrix * cmDATA_ini, cmatrix * cmDATA);
void FFT(int lim_NoB, int i_max, cplx * cW, cvector * cvDATA_ini, cvector * cvDATA);
void iFFT(int lim_NoB, int i_max, cplx * cW, cmatrix * cmDATA_ini, cmatrix * cmDATA);
void iFFT(int lim_NoB, int i_max, cplx * cW, cvector * cvDATA_ini, cvector * cvDATA);

void FFT(int i_max, real ** Im_cW, cplx ** cW2, cmatrix * cmDATA_ini,
         cmatrix * cmDATA, cmatrix * cA, cmatrix * cB, cmatrix * cC);
void FFT(int i_max, real ** Im_cW, cplx ** cW2, cvector * cvDATA_ini,
         cvector * cvDATA, cvector * cA, cvector * cB, cvector * cC);
void iFFT(int i_max, real ** Im_cW, cplx ** cW2, cvector * cvDATA_ini,
          cvector * cvDATA, cvector * cA, cvector * cB, cvector * cC);


void FFT_3D(int i_max, int j_max, int l_max, cplx * cW_i, cplx * cW_j,
            cplx * cW_l, cmatrix *** cmDATA);
void FFT_3D(int i_max, int j_max, int l_max, cplx * cW_2_i, cplx * cW_2_j,
            cplx * cW_2_l, cplx * cW_i, cplx * cW_j, cplx * cW_l,
            cvector *** cvDATA);
void iFFT_3D(int i_max, int j_max, int l_max, cplx * cW_i, cplx * cW_j,
             cplx * cW_l, cvector *** cvDATA);

void FFT_3D(int i_max, int j_max, int l_max, real ** Im_cW_i, real ** Im_cW_j,
            real ** Im_cW_l, cplx ** cW2_i, cplx ** cW2_j, cplx ** cW2_l,
            cmatrix *** cmDATA);
void FFT_3D(int i_max, int j_max, int l_max, real ** Im_cW_i, real ** Im_cW_j,
            real ** Im_cW_l, cplx ** cW2_i, cplx ** cW2_j, cplx ** cW2_l,
            cvector *** cvDATA);
void iFFT_3D(int i_max, int j_max, int l_max, real ** Im_cW_i, real ** Im_cW_j,
             real ** Im_cW_l, cplx ** cW2_i, cplx ** cW2_j, cplx ** cW2_l,
             cvector *** cvDATA);

void Convolution(int i_max, int j_max, int l_max, cplx * cW_2_i, cplx * cW_2_j,
                 cplx * cW_2_l, cplx * cW_i, cplx * cW_j, cplx * cW_l,
                 cmatrix *** cmDATA, cvector *** cvDATA);
void Convolution_conj(int i_max, int j_max, int l_max, cplx * cW_2_i,
                      cplx * cW_2_j, cplx * cW_2_l, cplx * cW_i,
                      cplx * cW_j, cplx * cW_l, cmatrix *** cmDATA,
                      cvector *** cvDATA);

void Convolution(int i_max, int j_max, int l_max, real ** Im_cW_i, real ** Im_cW_j,
                 real ** Im_cW_l, cplx ** cW2_i, cplx ** cW2_j, cplx ** cW2_l,
                 cmatrix *** cmDATA, cvector *** cvDATA);
void Convolution_conj(int i_max, int j_max, int l_max, real ** Im_cW_i,
                      real ** Im_cW_j, real ** Im_cW_l, cplx ** cW2_i,
                      cplx ** cW2_j, cplx ** cW2_l, cmatrix *** cmDATA,
                      cvector *** cvDATA);

/*********************************** DDA ************************************/

/* Scratch arrays used by the iterative solver. They are allocated once and
   reused for every illumination and every realization, so the program does
   not allocate and free hundreds of megabytes on each call. */

typedef
    struct
    {   int i_max;
        int j_max;
        int l_max;
        cvector *** cvP;        /* dipole moments, current iterate            */
        cvector *** cvP_min;    /* iterate with the smallest discrepancy      */
        cvector *** cvWork;     /* padded (2i x 2j x 2l) convolution buffer  */
        cvector *** cvR;        /* residual                                    */
        cvector *** cvS;        /* search direction                            */
        cvector *** cvG;        /* BiCGStab only, NULL otherwise               */
        cvector *** cvH;        /* BiCGStab only, NULL otherwise               */
    }
    dda_workspace;

void dda_workspace_create(dda_workspace * ws, int i_max, int j_max, int l_max,
                          int iMode_SolutionMethod);
void dda_workspace_free(dda_workspace * ws);

vector Rotation_x(real cos_angle, real sin_angle, vector coordinate_ini);
vector Rotation_y(real cos_angle, real sin_angle, vector coordinate_ini);
vector Rotation_z(real cos_angle, real sin_angle, vector coordinate_ini);
cvector Rotation_x(real cos_angle, real sin_angle, cvector coordinate_ini);
cvector Rotation_y(real cos_angle, real sin_angle, cvector coordinate_ini);
cvector Rotation_z(real cos_angle, real sin_angle, cvector coordinate_ini);

vector vrand(void);
void CreateAMatrix(int i_max, int j_max, int l_max, real d0, cmatrix *** cmDATA);

void FirstApproximation(int i_max, int j_max, int l_max, int * iListMV,
                        cplx * cListMaterials, char *** cAlpha, real d0,
                        vector kinc, cvector Einc, cvector *** cvB);

void Copying(int i_max, int j_max, int l_max, cvector *** cvA, cvector *** cvB);
void CopyingBack(int i_max, int j_max, int l_max, cvector *** cvA, cvector *** cvB);

ld_cplx ldcScalMul(int i_max, int j_max, int l_max, int * iListMV,
                   char *** cAlpha, cvector *** cv1, cvector *** cv2);
long double ldScalMulConj(int i_max, int j_max, int l_max, int * iListMV,
                          char *** cAlpha, cvector *** cv1, cvector *** cv2);
long double ldError(int i_max, int j_max, int l_max, int * iListMV,
                    char *** cAlpha, cvector *** cvR);

void LinearComb(int i_max, int j_max, int l_max, cplx cCoefficient, int * iListMV,
                char *** cAlpha, cvector *** cv1, cvector *** cv2, cvector *** cv3);
void LinearComb(int i_max, int j_max, int l_max, cplx cCoefficient1,
                cplx cCoefficient2, int * iListMV, char *** cAlpha,
                cvector *** cv1, cvector *** cv2, cvector *** cv3,
                cvector *** cv4);
void LinearComb(int i_max, int j_max, int l_max, cplx cCoefficient, int * iListMV,
                char *** cAlpha, cvector *** cv1, cvector *** cv2, int iMode);
void LinearComb(int i_max, int j_max, int l_max, real cCoefficient, int * iListMV,
                char *** cAlpha, cvector *** cv1, cvector *** cv2, int iMode);

int DDA_FFT(int iNoIllumination, int iMode_SolutionMethod, real rTau, int i_max,
            int j_max, int l_max, cplx * cListMaterials, int * iListMV,
            char *** cAlpha, int Max_Step, real Error_max, real d0, vector kinc,
            cvector Einc, cplx * cW_2_i, cplx * cW_2_j, cplx * cW_2_l, cplx * cW_i,
            cplx * cW_j, cplx * cW_l, real ** Im_cW_i, real ** Im_cW_j,
            real ** Im_cW_l, cplx ** cW2_i, cplx ** cW2_j, cplx ** cW2_l,
            int iMode_FFT, cmatrix *** cmA, int iNoScPl, int iNoA, vector ** ksc,
            cvector ** Esc, real * rCabs, real * rCext, dda_workspace * ws);

cvector FarField(int i_max, int j_max, int l_max, int * iListMV,
                 char *** cAlpha, real d0, cvector *** cvP, vector ksc_ort);

cplx Polarizability_CM(real d0, cplx cRefrIndex);
cplx Polarizability_CMRR(cplx alpha_CM);
cplx Polarizability_CMNR(real d0, cplx cRefrIndex);
cplx Polarizability_ILDR(real d0, cplx cRefrIndex, cplx alpha_CM);
cplx Polarizability_LDR(real d0, cplx cRefrIndex, cplx alpha_CM,
                        vector kinc_ort, cvector Einc_ort);

void CabsCext(real d0, vector kinc_ort, cvector Einc, int i_max, int j_max,
              int l_max, cplx * cListMaterials, char *** cAlpha,
              cvector *** cvP, real * Cabs, real * Cext);

/***************************** Generate particle ****************************/

int iGenerateEllipsoid(int i_max, int j_max, int l_max, real b_a, real c_a,
                       char *** cAlpha);

int iGenerateBisphere(int i_max, int j_max, int l_max, real R1_R2, char *** cAlpha);

int iGenerateParallelepiped(int i_max, int j_max, int l_max, real b_a, real c_a,
                            char *** cAlpha);

int iGenerateRandomEllipsoid(int i_max, int j_max, int l_max, real b_a, real c_a,
                             int iSorts, int * iSeeds, char *** cAlpha);

/********************** Generate anthrax-like particles *********************/

int iGenerateAnthrax(int i_max, int j_max, int l_max, real radius_worm,
                     real length_worm, int iWormsNumber, cylinder * worms,
                     char *** cAlpha);

int iGenerateHedgehog(int i_max, int j_max, int l_max, real radius_worm,
                      real length_worm, int iWormsNumber, cylinder * worms,
                      char *** cAlpha);

int iGenerateDeadAnthrax(int i_max, int j_max, int l_max, real radius_worm,
                         real length_worm, int iWormsNumber, cylinder * worms,
                         char *** cAlpha);

int iGenerateRandomParallelepiped(int i_max, int j_max, int l_max, real b_a,
                                  real c_a, int iSorts, int * iSeeds,
                                  char *** cAlpha);

int iGenerateNibbleEllipsoid(int i_max, int j_max, int l_max, real b_a, real c_a,
                             int iSorts, int * iSeeds, int iBites, real rDepth,
                             int * iListMV, char *** cAlpha);

int iGenerateHierarchicalEllipsoid(int iHierarchicalLevels, int i_max, int j_max, int l_max,
                                   real b_a, real c_a, int iSorts, int * iSeeds,
                                   char *** cAlpha);

int iGenerateSphereSpec(int i_max, int j_max, int l_max, real Radius, real rho,
                       char *** cAlpha);


void Particle_Boundary(int iNoD, int i_max_ini, int j_max_ini, int l_max_ini,
                       char *** cAlpha, cplx * cListMaterials, int * X_min,
                       int * Y_min, int * Z_min, int * X_max, int * Y_max,
                       int * Z_max);
