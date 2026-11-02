
#include "function.hpp"

/****************************************************************************/
/******************************** Functions *********************************/
/****************************************************************************/

vector Rotation_x(real cos_angle, real sin_angle, vector coordinate_ini)
{   vector vTmp;

    vTmp.x=coordinate_ini.x;
    vTmp.y=cos_angle*coordinate_ini.y+sin_angle*coordinate_ini.z;
    vTmp.z=-sin_angle*coordinate_ini.y+cos_angle*coordinate_ini.z;

    return vTmp;
}

vector Rotation_y(real cos_angle, real sin_angle, vector coordinate_ini)
{   vector vTmp;

    vTmp.x=cos_angle*coordinate_ini.x-sin_angle*coordinate_ini.z;
    vTmp.y=coordinate_ini.y;
    vTmp.z=sin_angle*coordinate_ini.x+cos_angle*coordinate_ini.z;

    return vTmp;
}

vector Rotation_z(real cos_angle, real sin_angle, vector coordinate_ini)
{   vector vTmp;

    vTmp.x=cos_angle*coordinate_ini.x+sin_angle*coordinate_ini.y;
    vTmp.y=-sin_angle*coordinate_ini.x+cos_angle*coordinate_ini.y;
    vTmp.z=coordinate_ini.z;

    return vTmp;
}

cvector Rotation_x(real cos_angle, real sin_angle, cvector coordinate_ini)
{   cvector cvTmp;

    cvTmp.x=coordinate_ini.x;
    cvTmp.y=cos_angle*coordinate_ini.y+sin_angle*coordinate_ini.z;
    cvTmp.z=-sin_angle*coordinate_ini.y+cos_angle*coordinate_ini.z;

    return cvTmp;
}

cvector Rotation_y(real cos_angle, real sin_angle, cvector coordinate_ini)
{   cvector cvTmp;

    cvTmp.x=cos_angle*coordinate_ini.x-sin_angle*coordinate_ini.z;
    cvTmp.y=coordinate_ini.y;
    cvTmp.z=sin_angle*coordinate_ini.x+cos_angle*coordinate_ini.z;

    return cvTmp;
}

cvector Rotation_z(real cos_angle, real sin_angle, cvector coordinate_ini)
{   cvector cvTmp;

    cvTmp.x=cos_angle*coordinate_ini.x+sin_angle*coordinate_ini.y;
    cvTmp.y=-sin_angle*coordinate_ini.x+cos_angle*coordinate_ini.y;
    cvTmp.z=coordinate_ini.z;

    return cvTmp;
}

vector vrand(void)
{   real alpha,beta,cos_beta,sin_beta;
    vector vTmp;

    alpha=(real)rand()*2*PI/RAND_MAX;
    cos_beta=2.*(real)rand()/(real)RAND_MAX-1.;
    beta=acos(cos_beta);
    sin_beta=sin(beta);

    // forming of the normalized wave vector
    vTmp.x=cos(alpha)*sin_beta;
    vTmp.y=sin(alpha)*sin_beta;
    vTmp.z=cos_beta;

    return vTmp;
}

void CreateAMatrix(int i_max, int j_max, int l_max, real d0, cmatrix *** cmDATA)
{   int i,j,l;
    real r1,r2,r3;
    cplx c0,c1,c2,c3,c4,c5;
    vector r;
    matrix rr;

    for(i=0;i<2*i_max;i++)
    {   if(i<i_max)
            r.x=(real)(i)*d0;
        else
            r.x=(real)(i-2*i_max)*d0;

        for(j=0;j<2*j_max;j++)
        {   if(j<j_max)
                r.y=(real)(j)*d0;
            else
                r.y=(real)(j-2*j_max)*d0;

            for(l=0;l<2*l_max;l++)
            {   if(l<l_max)
                    r.z=(real)(l)*d0;
                else
                    r.z=(real)(l-2*l_max)*d0;

                rr=dyad(r);
                r2=rr.e00+rr.e11+rr.e22;
                if(r2<ZERO)
                {   cmDATA[i][j][l].e00.r=0.;   cmDATA[i][j][l].e00.i=0.;
                    cmDATA[i][j][l].e01.r=0.;   cmDATA[i][j][l].e01.i=0.;
                    cmDATA[i][j][l].e02.r=0.;   cmDATA[i][j][l].e02.i=0.;
                    cmDATA[i][j][l].e11.r=0.;   cmDATA[i][j][l].e11.i=0.;
                    cmDATA[i][j][l].e12.r=0.;   cmDATA[i][j][l].e12.i=0.;
                    cmDATA[i][j][l].e22.r=0.;   cmDATA[i][j][l].e22.i=0.;
                }
                else
                {   r1=sqrt(r2);

                    c0.r=r2-1.; c0.i=r1;
                    c1.r=3.-r2; c1.i=-3.*r1;
                    r3=r1*r2;
                    Cexp_im(r1,c3);
                    CDdiv(c3,r3,c2);

                    Cmul(c0,c2,c3);
                    Cmul(c1,c2,c4);
                    CDdiv(c4,r2,c5);
                    DCmul(rr.e00,c5,cmDATA[i][j][l].e00);
                    DCmul(rr.e01,c5,cmDATA[i][j][l].e01);
                    DCmul(rr.e02,c5,cmDATA[i][j][l].e02);
                    DCmul(rr.e11,c5,cmDATA[i][j][l].e11);
                    DCmul(rr.e12,c5,cmDATA[i][j][l].e12);
                    DCmul(rr.e22,c5,cmDATA[i][j][l].e22);

                    Cadd(cmDATA[i][j][l].e00,c3,cmDATA[i][j][l].e00);
                    Cadd(cmDATA[i][j][l].e11,c3,cmDATA[i][j][l].e11);
                    Cadd(cmDATA[i][j][l].e22,c3,cmDATA[i][j][l].e22);
                }
            }
        }
    }

    return;
}

void Convolution(int i_max, int j_max, int l_max, cplx * cW_2_i, cplx * cW_2_j,
                 cplx * cW_2_l, cplx * cW_i, cplx * cW_j, cplx * cW_l,
                 cmatrix *** cmDATA, cvector *** cvDATA)
{   int i,j,l;
    cplx c0,c1,c2,c3,c4,c5,c6;
    FFT_3D(i_max,j_max,l_max,cW_2_i,cW_2_j,cW_2_l,cW_i,cW_j,cW_l,cvDATA);

    for(i=0;i<i_max;i++)
        for(j=0;j<j_max;j++)
            for(l=0;l<l_max;l++)
            {   if(cvDATA[i][j][l].x.r!=0.||cvDATA[i][j][l].x.i!=0.||
                   cvDATA[i][j][l].y.r!=0.||cvDATA[i][j][l].y.i!=0.||
                   cvDATA[i][j][l].z.r!=0.||cvDATA[i][j][l].z.i!=0.)
                {   Cmul(cmDATA[i][j][l].e00,cvDATA[i][j][l].x,c0);
                    Cmul(cmDATA[i][j][l].e01,cvDATA[i][j][l].y,c1);
                    Cmul(cmDATA[i][j][l].e02,cvDATA[i][j][l].z,c2);
                    Cadd(c0,c1,c3); Cadd(c2,c3,c4);

                    Cmul(cmDATA[i][j][l].e01,cvDATA[i][j][l].x,c0);
                    Cmul(cmDATA[i][j][l].e11,cvDATA[i][j][l].y,c1);
                    Cmul(cmDATA[i][j][l].e12,cvDATA[i][j][l].z,c2);
                    Cadd(c0,c1,c3); Cadd(c2,c3,c5);

                    Cmul(cmDATA[i][j][l].e02,cvDATA[i][j][l].x,c0);
                    Cmul(cmDATA[i][j][l].e12,cvDATA[i][j][l].y,c1);
                    Cmul(cmDATA[i][j][l].e22,cvDATA[i][j][l].z,c2);
                    Cadd(c0,c1,c3); Cadd(c2,c3,c6);

                    Ccopy(c4,cvDATA[i][j][l].x);
                    Ccopy(c5,cvDATA[i][j][l].y);
                    Ccopy(c6,cvDATA[i][j][l].z);
                }
            }

    iFFT_3D(i_max,j_max,l_max,cW_i,cW_j,cW_l,cvDATA);

    return;
}

void Convolution(int i_max, int j_max, int l_max, real ** Im_cW_i, real ** Im_cW_j,
                 real ** Im_cW_l, cplx ** cW2_i, cplx ** cW2_j, cplx ** cW2_l,
                 cmatrix *** cmDATA, cvector *** cvDATA)
{   int i,j,l;
    cplx c0,c1,c2,c3,c4,c5,c6;
    FFT_3D(i_max,j_max,l_max,Im_cW_i,Im_cW_j,Im_cW_l,cW2_i,cW2_j,cW2_l,cvDATA);

    for(i=0;i<i_max;i++)
        for(j=0;j<j_max;j++)
            for(l=0;l<l_max;l++)
            {   if(cvDATA[i][j][l].x.r!=0.||cvDATA[i][j][l].x.i!=0.||
                   cvDATA[i][j][l].y.r!=0.||cvDATA[i][j][l].y.i!=0.||
                   cvDATA[i][j][l].z.r!=0.||cvDATA[i][j][l].z.i!=0.)
                {   Cmul(cmDATA[i][j][l].e00,cvDATA[i][j][l].x,c0);
                    Cmul(cmDATA[i][j][l].e01,cvDATA[i][j][l].y,c1);
                    Cmul(cmDATA[i][j][l].e02,cvDATA[i][j][l].z,c2);
                    Cadd(c0,c1,c3); Cadd(c2,c3,c4);

                    Cmul(cmDATA[i][j][l].e01,cvDATA[i][j][l].x,c0);
                    Cmul(cmDATA[i][j][l].e11,cvDATA[i][j][l].y,c1);
                    Cmul(cmDATA[i][j][l].e12,cvDATA[i][j][l].z,c2);
                    Cadd(c0,c1,c3); Cadd(c2,c3,c5);

                    Cmul(cmDATA[i][j][l].e02,cvDATA[i][j][l].x,c0);
                    Cmul(cmDATA[i][j][l].e12,cvDATA[i][j][l].y,c1);
                    Cmul(cmDATA[i][j][l].e22,cvDATA[i][j][l].z,c2);
                    Cadd(c0,c1,c3); Cadd(c2,c3,c6);

                    Ccopy(c4,cvDATA[i][j][l].x);
                    Ccopy(c5,cvDATA[i][j][l].y);
                    Ccopy(c6,cvDATA[i][j][l].z);
                }
            }

    iFFT_3D(i_max,j_max,l_max,Im_cW_i,Im_cW_j,Im_cW_l,cW2_i,cW2_j,cW2_l,cvDATA);
    return;
}

void Convolution_conj(int i_max, int j_max, int l_max, cplx * cW_2_i,
                      cplx * cW_2_j, cplx * cW_2_l, cplx * cW_i,
                      cplx * cW_j, cplx * cW_l, cmatrix *** cmDATA,
                      cvector *** cvDATA)
{   int i,j,l;
    cplx c0,c1,c2,c3,c4,c5,c6;
    cmatrix cm0;

    FFT_3D(i_max,j_max,l_max,cW_2_i,cW_2_j,cW_2_l,cW_i,cW_j,cW_l,cvDATA);

    cm0=conj(cmDATA[0][0][0]);

    // single 0,0,0 element
    Cmul(cm0.e00,cvDATA[0][0][0].x,c0);
    Cmul(cm0.e01,cvDATA[0][0][0].y,c1);
    Cmul(cm0.e02,cvDATA[0][0][0].z,c2);
    Cadd(c0,c1,c3); Cadd(c2,c3,c4);

    Cmul(cm0.e01,cvDATA[0][0][0].x,c0);
    Cmul(cm0.e11,cvDATA[0][0][0].y,c1);
    Cmul(cm0.e12,cvDATA[0][0][0].z,c2);
    Cadd(c0,c1,c3); Cadd(c2,c3,c5);

    Cmul(cm0.e02,cvDATA[0][0][0].x,c0);
    Cmul(cm0.e12,cvDATA[0][0][0].y,c1);
    Cmul(cm0.e22,cvDATA[0][0][0].z,c2);
    Cadd(c0,c1,c3); Cadd(c2,c3,c6);

    Ccopy(c4,cvDATA[0][0][0].x);
    Ccopy(c5,cvDATA[0][0][0].y);
    Ccopy(c6,cvDATA[0][0][0].z);

    // three edges beginnig from poin 0,0,0
    j=0;    l=0;
    for(i=1;i<i_max;i++)
    {   cm0=conj(cmDATA[i_max-i][j][l]);

        Cmul(cm0.e00,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e01,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e02,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c4);

        Cmul(cm0.e01,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e11,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e12,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c5);

        Cmul(cm0.e02,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e12,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e22,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c6);

        Ccopy(c4,cvDATA[i][j][l].x);
        Ccopy(c5,cvDATA[i][j][l].y);
        Ccopy(c6,cvDATA[i][j][l].z);
    }

    i=0;    l=0;
    for(j=1;j<j_max;j++)
    {   cm0=conj(cmDATA[i][j_max-j][l]);

        Cmul(cm0.e00,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e01,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e02,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c4);

        Cmul(cm0.e01,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e11,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e12,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c5);

        Cmul(cm0.e02,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e12,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e22,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c6);

        Ccopy(c4,cvDATA[i][j][l].x);
        Ccopy(c5,cvDATA[i][j][l].y);
        Ccopy(c6,cvDATA[i][j][l].z);
    }

    i=0;    j=0;
    for(l=1;l<l_max;l++)
    {   cm0=conj(cmDATA[i][j][l_max-l]);

        Cmul(cm0.e00,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e01,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e02,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c4);

        Cmul(cm0.e01,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e11,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e12,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c5);

        Cmul(cm0.e02,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e12,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e22,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c6);

        Ccopy(c4,cvDATA[i][j][l].x);
        Ccopy(c5,cvDATA[i][j][l].y);
        Ccopy(c6,cvDATA[i][j][l].z);
    }

    // three sides beginnig from poin 0,0,0
    l=0;
    for(i=1;i<i_max;i++)
        for(j=1;j<j_max;j++)
        {   cm0=conj(cmDATA[i_max-i][j_max-j][l]);

            Cmul(cm0.e00,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e01,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e02,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c4);

            Cmul(cm0.e01,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e11,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e12,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c5);

            Cmul(cm0.e02,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e12,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e22,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c6);

            Ccopy(c4,cvDATA[i][j][l].x);
            Ccopy(c5,cvDATA[i][j][l].y);
            Ccopy(c6,cvDATA[i][j][l].z);
        }

    j=0;
    for(i=1;i<i_max;i++)
        for(l=1;l<l_max;l++)
        {   cm0=conj(cmDATA[i_max-i][j][l_max-l]);

            Cmul(cm0.e00,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e01,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e02,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c4);

            Cmul(cm0.e01,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e11,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e12,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c5);

            Cmul(cm0.e02,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e12,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e22,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c6);

            Ccopy(c4,cvDATA[i][j][l].x);
            Ccopy(c5,cvDATA[i][j][l].y);
            Ccopy(c6,cvDATA[i][j][l].z);
        }

    i=0;
    for(j=1;j<j_max;j++)
        for(l=1;l<l_max;l++)
        {   cm0=conj(cmDATA[i][j_max-j][l_max-l]);

            Cmul(cm0.e00,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e01,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e02,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c4);

            Cmul(cm0.e01,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e11,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e12,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c5);

            Cmul(cm0.e02,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e12,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e22,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c6);

            Ccopy(c4,cvDATA[i][j][l].x);
            Ccopy(c5,cvDATA[i][j][l].y);
            Ccopy(c6,cvDATA[i][j][l].z);
        }

    // others elements
    for(i=1;i<i_max;i++)
        for(j=1;j<j_max;j++)
            for(l=1;l<l_max;l++)
            {   cm0=conj(cmDATA[i_max-i][j_max-j][l_max-l]);

                Cmul(cm0.e00,cvDATA[i][j][l].x,c0);
                Cmul(cm0.e01,cvDATA[i][j][l].y,c1);
                Cmul(cm0.e02,cvDATA[i][j][l].z,c2);
                Cadd(c0,c1,c3); Cadd(c2,c3,c4);

                Cmul(cm0.e01,cvDATA[i][j][l].x,c0);
                Cmul(cm0.e11,cvDATA[i][j][l].y,c1);
                Cmul(cm0.e12,cvDATA[i][j][l].z,c2);
                Cadd(c0,c1,c3); Cadd(c2,c3,c5);

                Cmul(cm0.e02,cvDATA[i][j][l].x,c0);
                Cmul(cm0.e12,cvDATA[i][j][l].y,c1);
                Cmul(cm0.e22,cvDATA[i][j][l].z,c2);
                Cadd(c0,c1,c3); Cadd(c2,c3,c6);

                Ccopy(c4,cvDATA[i][j][l].x);
                Ccopy(c5,cvDATA[i][j][l].y);
                Ccopy(c6,cvDATA[i][j][l].z);
            }

    iFFT_3D(i_max,j_max,l_max,cW_i,cW_j,cW_l,cvDATA);

    return;
}

void Convolution_conj(int i_max, int j_max, int l_max, real ** Im_cW_i,
                      real ** Im_cW_j, real ** Im_cW_l, cplx ** cW2_i,
                      cplx ** cW2_j, cplx ** cW2_l, cmatrix *** cmDATA,
                      cvector *** cvDATA)
{   int i,j,l;
    cplx c0,c1,c2,c3,c4,c5,c6;
    cmatrix cm0;

    FFT_3D(i_max,j_max,l_max,Im_cW_i,Im_cW_j,Im_cW_l,cW2_i,cW2_j,cW2_l,cvDATA);

    cm0=conj(cmDATA[0][0][0]);

    // single 0,0,0 element
    Cmul(cm0.e00,cvDATA[0][0][0].x,c0);
    Cmul(cm0.e01,cvDATA[0][0][0].y,c1);
    Cmul(cm0.e02,cvDATA[0][0][0].z,c2);
    Cadd(c0,c1,c3); Cadd(c2,c3,c4);

    Cmul(cm0.e01,cvDATA[0][0][0].x,c0);
    Cmul(cm0.e11,cvDATA[0][0][0].y,c1);
    Cmul(cm0.e12,cvDATA[0][0][0].z,c2);
    Cadd(c0,c1,c3); Cadd(c2,c3,c5);

    Cmul(cm0.e02,cvDATA[0][0][0].x,c0);
    Cmul(cm0.e12,cvDATA[0][0][0].y,c1);
    Cmul(cm0.e22,cvDATA[0][0][0].z,c2);
    Cadd(c0,c1,c3); Cadd(c2,c3,c6);

    Ccopy(c4,cvDATA[0][0][0].x);
    Ccopy(c5,cvDATA[0][0][0].y);
    Ccopy(c6,cvDATA[0][0][0].z);

    // three edges beginnig from poin 0,0,0
    j=0;    l=0;
    for(i=1;i<i_max;i++)
    {   cm0=conj(cmDATA[i_max-i][j][l]);

        Cmul(cm0.e00,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e01,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e02,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c4);

        Cmul(cm0.e01,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e11,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e12,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c5);

        Cmul(cm0.e02,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e12,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e22,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c6);

        Ccopy(c4,cvDATA[i][j][l].x);
        Ccopy(c5,cvDATA[i][j][l].y);
        Ccopy(c6,cvDATA[i][j][l].z);
    }

    i=0;    l=0;
    for(j=1;j<j_max;j++)
    {   cm0=conj(cmDATA[i][j_max-j][l]);

        Cmul(cm0.e00,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e01,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e02,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c4);

        Cmul(cm0.e01,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e11,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e12,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c5);

        Cmul(cm0.e02,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e12,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e22,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c6);

        Ccopy(c4,cvDATA[i][j][l].x);
        Ccopy(c5,cvDATA[i][j][l].y);
        Ccopy(c6,cvDATA[i][j][l].z);
    }

    i=0;    j=0;
    for(l=1;l<l_max;l++)
    {   cm0=conj(cmDATA[i][j][l_max-l]);

        Cmul(cm0.e00,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e01,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e02,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c4);

        Cmul(cm0.e01,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e11,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e12,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c5);

        Cmul(cm0.e02,cvDATA[i][j][l].x,c0);
        Cmul(cm0.e12,cvDATA[i][j][l].y,c1);
        Cmul(cm0.e22,cvDATA[i][j][l].z,c2);
        Cadd(c0,c1,c3); Cadd(c2,c3,c6);

        Ccopy(c4,cvDATA[i][j][l].x);
        Ccopy(c5,cvDATA[i][j][l].y);
        Ccopy(c6,cvDATA[i][j][l].z);
    }

    // three sides beginnig from poin 0,0,0
    l=0;
    for(i=1;i<i_max;i++)
        for(j=1;j<j_max;j++)
        {   cm0=conj(cmDATA[i_max-i][j_max-j][l]);

            Cmul(cm0.e00,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e01,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e02,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c4);

            Cmul(cm0.e01,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e11,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e12,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c5);

            Cmul(cm0.e02,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e12,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e22,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c6);

            Ccopy(c4,cvDATA[i][j][l].x);
            Ccopy(c5,cvDATA[i][j][l].y);
            Ccopy(c6,cvDATA[i][j][l].z);
        }

    j=0;
    for(i=1;i<i_max;i++)
        for(l=1;l<l_max;l++)
        {   cm0=conj(cmDATA[i_max-i][j][l_max-l]);

            Cmul(cm0.e00,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e01,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e02,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c4);

            Cmul(cm0.e01,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e11,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e12,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c5);

            Cmul(cm0.e02,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e12,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e22,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c6);

            Ccopy(c4,cvDATA[i][j][l].x);
            Ccopy(c5,cvDATA[i][j][l].y);
            Ccopy(c6,cvDATA[i][j][l].z);
        }

    i=0;
    for(j=1;j<j_max;j++)
        for(l=1;l<l_max;l++)
        {   cm0=conj(cmDATA[i][j_max-j][l_max-l]);

            Cmul(cm0.e00,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e01,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e02,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c4);

            Cmul(cm0.e01,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e11,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e12,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c5);

            Cmul(cm0.e02,cvDATA[i][j][l].x,c0);
            Cmul(cm0.e12,cvDATA[i][j][l].y,c1);
            Cmul(cm0.e22,cvDATA[i][j][l].z,c2);
            Cadd(c0,c1,c3); Cadd(c2,c3,c6);

            Ccopy(c4,cvDATA[i][j][l].x);
            Ccopy(c5,cvDATA[i][j][l].y);
            Ccopy(c6,cvDATA[i][j][l].z);
        }

    // others elements
    for(i=1;i<i_max;i++)
        for(j=1;j<j_max;j++)
            for(l=1;l<l_max;l++)
            {   cm0=conj(cmDATA[i_max-i][j_max-j][l_max-l]);

                Cmul(cm0.e00,cvDATA[i][j][l].x,c0);
                Cmul(cm0.e01,cvDATA[i][j][l].y,c1);
                Cmul(cm0.e02,cvDATA[i][j][l].z,c2);
                Cadd(c0,c1,c3); Cadd(c2,c3,c4);

                Cmul(cm0.e01,cvDATA[i][j][l].x,c0);
                Cmul(cm0.e11,cvDATA[i][j][l].y,c1);
                Cmul(cm0.e12,cvDATA[i][j][l].z,c2);
                Cadd(c0,c1,c3); Cadd(c2,c3,c5);

                Cmul(cm0.e02,cvDATA[i][j][l].x,c0);
                Cmul(cm0.e12,cvDATA[i][j][l].y,c1);
                Cmul(cm0.e22,cvDATA[i][j][l].z,c2);
                Cadd(c0,c1,c3); Cadd(c2,c3,c6);

                Ccopy(c4,cvDATA[i][j][l].x);
                Ccopy(c5,cvDATA[i][j][l].y);
                Ccopy(c6,cvDATA[i][j][l].z);
            }

    iFFT_3D(i_max,j_max,l_max,Im_cW_i,Im_cW_j,Im_cW_l,cW2_i,cW2_j,cW2_l,cvDATA);
    return;
}

/****************************************************************************/
/****************************** First step **********************************/
/****************************************************************************/

void FirstApproximation(int i_max, int j_max, int l_max, int * iListMV,
                        cplx * cListMaterials, char *** cAlpha, real d0,
                        vector kinc, cvector Einc, cvector *** cvB)
{   int i,j,l,m;
    cplx c0,c1;
    vector r,vCenter={(real)((i_max-1)*d0/2.),(real)((j_max-1)*d0/2.),(real)((l_max-1)*d0/2.)};

    for(i=0;i<i_max;i++)
    {   r.x=(real)i*d0;
        for(j=0;j<j_max;j++)
        {   r.y=(real)j*d0;
            for(l=0;l<l_max;l++)
            {   r.z=(real)l*d0;

                //cvB[i][j][l]=(cListMaterials[(int)(cAlpha[i][j][l])]
                            //*exp_i((kinc*(r-vCenter))))*Einc;

                m=(int)cAlpha[i][j][l];
                if(iListMV[m]!=0)
                {   c0=exp_i((kinc*(r-vCenter)));
                    Cmul(c0,cListMaterials[m],c1);
                    Cmul(c1,Einc.x,cvB[i][j][l].x);
                    Cmul(c1,Einc.y,cvB[i][j][l].y);
                    Cmul(c1,Einc.z,cvB[i][j][l].z);
                }
                else
                {   cvB[i][j][l].x.r=0.;    cvB[i][j][l].x.i=0.;
                    cvB[i][j][l].y.r=0.;    cvB[i][j][l].y.i=0.;
                    cvB[i][j][l].z.r=0.;    cvB[i][j][l].z.i=0.;
                }
            }
        }
    }

    return;
}

void Copying(int i_max, int j_max, int l_max, cvector *** cvA, cvector *** cvB)
{   int i,j,l;

    for(i=0;i<i_max;i++)
        for(j=0;j<j_max;j++)
            for(l=0;l<l_max;l++)
                cvB[i][j][l]=cvA[i][j][l];

    return;
}

void CopyingBack(int i_max, int j_max, int l_max, cvector *** cvA,
                 cvector *** cvB)
{   int i,j,l;

    for(i=0;i<2*i_max;i++)
        for(j=0;j<2*j_max;j++)
            for(l=0;l<2*l_max;l++)
                if(i>=i_max||j>=j_max||l>=l_max)
                {   cvB[i][j][l].x=zero();
                    cvB[i][j][l].y=zero();
                    cvB[i][j][l].z=zero();
                }
                else
                    cvB[i][j][l]=cvA[i][j][l];

    return;
}

/****************************************************************************/
/***************************** Others steps *********************************/
/****************************************************************************/

ld_cplx ldcScalMul(int i_max, int j_max, int l_max, int * iListMV, char *** cAlpha,
                 cvector *** cv1, cvector *** cv2)
{   int i,j,l,m;
    ld_cplx cTmp={0.,0.},c0,c1,c2;

    for(i=0;i<i_max;i++)
        for(j=0;j<j_max;j++)
            for(l=0;l<l_max;l++)
            {   m=(int)cAlpha[i][j][l];
                if(iListMV[m]!=0)
                {   c0.r=(long double)(cv1[i][j][l].x.r*cv2[i][j][l].x.r
                        -cv1[i][j][l].x.i*cv2[i][j][l].x.i);
                    c0.i=(long double)(cv1[i][j][l].x.r*cv2[i][j][l].x.i
                        +cv1[i][j][l].x.i*cv2[i][j][l].x.r);

                    c1.r=(long double)(cv1[i][j][l].y.r*cv2[i][j][l].y.r
                        -cv1[i][j][l].y.i*cv2[i][j][l].y.i);
                    c1.i=(long double)(cv1[i][j][l].y.r*cv2[i][j][l].y.i
                        +cv1[i][j][l].y.i*cv2[i][j][l].y.r);

                    c2.r=(long double)(cv1[i][j][l].z.r*cv2[i][j][l].z.r
                        -cv1[i][j][l].z.i*cv2[i][j][l].z.i);
                    c2.i=(long double)(cv1[i][j][l].z.r*cv2[i][j][l].z.i
                        +cv1[i][j][l].z.i*cv2[i][j][l].z.r);

                    cTmp.r+=(c0.r+c1.r+c2.r);
                    cTmp.i+=(c0.i+c1.i+c2.i);
                }
            }

    return cTmp;
}

long double ldScalMulConj(int i_max, int j_max, int l_max, int * iListMV,
                         char *** cAlpha, cvector *** cv1, cvector *** cv2)
{   int i,j,l,m;
    long double rTmp=0.,ld0,ld1,ld2;

    for(i=0;i<i_max;i++)
        for(j=0;j<j_max;j++)
            for(l=0;l<l_max;l++)
            {   m=(int)cAlpha[i][j][l];
                if(iListMV[m]!=0)
                {
                    ld0=(long double)(cv1[i][j][l].x.r*cv2[i][j][l].x.r
                        +cv1[i][j][l].x.i*cv2[i][j][l].x.i);

                    ld1=(long double)(cv1[i][j][l].y.r*cv2[i][j][l].y.r
                        +cv1[i][j][l].y.i*cv2[i][j][l].y.i);

                    ld2=(long double)(cv1[i][j][l].z.r*cv2[i][j][l].z.r
                        +cv1[i][j][l].z.i*cv2[i][j][l].z.i);

                    rTmp+=(ld0+ld1+ld2);
                }
            }

    return rTmp;
}

long double ldError(int i_max, int j_max, int l_max, int * iListMV, char *** cAlpha,
            cvector *** cvR)
{   int i,j,l,m;
    long int n=0;
    long double rTmp=0.,ld0,ld1,ld2;

    for(i=0;i<i_max;i++)
        for(j=0;j<j_max;j++)
            for(l=0;l<l_max;l++)
            {   m=(int)cAlpha[i][j][l];
                if(iListMV[m]!=0)
                {   n++;

                    ld0=(long double)(cvR[i][j][l].x.r*cvR[i][j][l].x.r
                        +cvR[i][j][l].x.i*cvR[i][j][l].x.i);

                    ld1=(long double)(cvR[i][j][l].y.r*cvR[i][j][l].y.r
                        +cvR[i][j][l].y.i*cvR[i][j][l].y.i);

                    ld2=(long double)(cvR[i][j][l].z.r*cvR[i][j][l].z.r
                        +cvR[i][j][l].z.i*cvR[i][j][l].z.i);

                    rTmp=rTmp+(long double)sqrt(ld0+ld1+ld2);
                }
            }

    return (rTmp/(long double)n);
}

void LinearComb(int i_max, int j_max, int l_max, cplx cCoefficient, int * iListMV,
                char *** cAlpha, cvector *** cv1, cvector *** cv2, cvector *** cv3)
{   int i,j,l,m;
    cvector cv0;

    for(i=0;i<i_max;i++)
        for(j=0;j<j_max;j++)
            for(l=0;l<l_max;l++)
            {   //cv3[i][j][l]=cv1[i][j][l]+cCoefficient*cv2[i][j][l];

                m=(int)cAlpha[i][j][l];
                if(iListMV[m]!=0)
                {   Cmul(cCoefficient,cv2[i][j][l].x,cv0.x);
                    Cmul(cCoefficient,cv2[i][j][l].y,cv0.y);
                    Cmul(cCoefficient,cv2[i][j][l].z,cv0.z);

                    Cadd(cv1[i][j][l].x,cv0.x,cv3[i][j][l].x);
                    Cadd(cv1[i][j][l].y,cv0.y,cv3[i][j][l].y);
                    Cadd(cv1[i][j][l].z,cv0.z,cv3[i][j][l].z);
                }
                else
                {   cv3[i][j][l].x.r=0.;    cv3[i][j][l].x.i=0.;
                    cv3[i][j][l].y.r=0.;    cv3[i][j][l].y.i=0.;
                    cv3[i][j][l].z.r=0.;    cv3[i][j][l].z.i=0.;
                }
            }

    return;
}

void LinearComb(int i_max, int j_max, int l_max, cplx cCoefficient1,
                cplx cCoefficient2, int * iListMV, char *** cAlpha,
                cvector *** cv1, cvector *** cv2, cvector *** cv3,
                cvector *** cv4)
{   int i,j,l,m;
    cvector cv00,cv01,cv02;

    for(i=0;i<i_max;i++)
        for(j=0;j<j_max;j++)
            for(l=0;l<l_max;l++)
            {   //cv4[i][j][l]=cv1[i][j][l]+cCoefficient1*cv2[i][j][l]
                              //+cCoefficient2*cv3[i][j][l];

                m=(int)cAlpha[i][j][l];
                if(iListMV[m]!=0)
                {   Cmul(cCoefficient1,cv2[i][j][l].x,cv00.x);
                    Cmul(cCoefficient1,cv2[i][j][l].y,cv00.y);
                    Cmul(cCoefficient1,cv2[i][j][l].z,cv00.z);

                    Cmul(cCoefficient2,cv3[i][j][l].x,cv01.x);
                    Cmul(cCoefficient2,cv3[i][j][l].y,cv01.y);
                    Cmul(cCoefficient2,cv3[i][j][l].z,cv01.z);

                    Cadd(cv1[i][j][l].x,cv00.x,cv02.x);
                    Cadd(cv1[i][j][l].y,cv00.y,cv02.y);
                    Cadd(cv1[i][j][l].z,cv00.z,cv02.z);

                    Cadd(cv01.x,cv02.x,cv4[i][j][l].x);
                    Cadd(cv01.y,cv02.y,cv4[i][j][l].y);
                    Cadd(cv01.z,cv02.z,cv4[i][j][l].z);
                }
                else
                {   cv4[i][j][l].x.r=0.;    cv4[i][j][l].x.i=0.;
                    cv4[i][j][l].y.r=0.;    cv4[i][j][l].y.i=0.;
                    cv4[i][j][l].z.r=0.;    cv4[i][j][l].z.i=0.;
                }
            }

    return;
}

void LinearComb(int i_max, int j_max, int l_max, cplx cCoefficient, int * iListMV,
                char *** cAlpha, cvector *** cv1, cvector *** cv2, int iMode)
{   int i,j,l,m;
    cplx c0,c1,c2;

    if(iMode==1)
    {   for(i=0;i<i_max;i++)
            for(j=0;j<j_max;j++)
                for(l=0;l<l_max;l++)
                {   //cv1[i][j][l]=cv1[i][j][l]+cCoefficient*cv2[i][j][l];

                    m=(int)cAlpha[i][j][l];
                    if(iListMV[m]!=0)
                    {   Cmul(cCoefficient,cv2[i][j][l].x,c0);
                        Cmul(cCoefficient,cv2[i][j][l].y,c1);
                        Cmul(cCoefficient,cv2[i][j][l].z,c2);

                        Cadd(cv1[i][j][l].x,c0,cv1[i][j][l].x);
                        Cadd(cv1[i][j][l].y,c1,cv1[i][j][l].y);
                        Cadd(cv1[i][j][l].z,c2,cv1[i][j][l].z);
                    }
                }
    }

    if(iMode==2)
    {   for(i=0;i<i_max;i++)
            for(j=0;j<j_max;j++)
                for(l=0;l<l_max;l++)
                {   //cv2[i][j][l]=cv1[i][j][l]+cCoefficient*cv2[i][j][l];

                    m=(int)cAlpha[i][j][l];
                    if(iListMV[m]!=0)
                    {   Cmul(cCoefficient,cv2[i][j][l].x,c0);
                        Cmul(cCoefficient,cv2[i][j][l].y,c1);
                        Cmul(cCoefficient,cv2[i][j][l].z,c2);

                        Cadd(cv1[i][j][l].x,c0,cv2[i][j][l].x);
                        Cadd(cv1[i][j][l].y,c1,cv2[i][j][l].y);
                        Cadd(cv1[i][j][l].z,c2,cv2[i][j][l].z);
                    }
                }
    }

    return;
}

void LinearComb(int i_max, int j_max, int l_max, real cCoefficient, int * iListMV,
                char *** cAlpha, cvector *** cv1, cvector *** cv2, int iMode)
{   int i,j,l,m;
    cplx c0,c1,c2;

    if(iMode==1)
    {   for(i=0;i<i_max;i++)
            for(j=0;j<j_max;j++)
                for(l=0;l<l_max;l++)
                {   //cv1[i][j][l]=cv1[i][j][l]+cCoefficient*cv2[i][j][l];

                    m=(int)cAlpha[i][j][l];
                    if(iListMV[m]!=0)
                    {   DCmul(cCoefficient,cv2[i][j][l].x,c0);
                        DCmul(cCoefficient,cv2[i][j][l].y,c1);
                        DCmul(cCoefficient,cv2[i][j][l].z,c2);

                        Cadd(cv1[i][j][l].x,c0,cv1[i][j][l].x);
                        Cadd(cv1[i][j][l].y,c1,cv1[i][j][l].y);
                        Cadd(cv1[i][j][l].z,c2,cv1[i][j][l].z);
                    }
                }
    }

    if(iMode==2)
    {   for(i=0;i<i_max;i++)
            for(j=0;j<j_max;j++)
                for(l=0;l<l_max;l++)
                {   //cv2[i][j][l]=cv1[i][j][l]+cCoefficient*cv2[i][j][l];

                    m=(int)cAlpha[i][j][l];
                    if(iListMV[m]!=0)
                    {   DCmul(cCoefficient,cv2[i][j][l].x,c0);
                        DCmul(cCoefficient,cv2[i][j][l].y,c1);
                        DCmul(cCoefficient,cv2[i][j][l].z,c2);

                        Cadd(cv1[i][j][l].x,c0,cv2[i][j][l].x);
                        Cadd(cv1[i][j][l].y,c1,cv2[i][j][l].y);
                        Cadd(cv1[i][j][l].z,c2,cv2[i][j][l].z);
                    }
                }
    }

    return;
}

/****************************************************************************/
/***************************** Memory for fields *****************************/
/****************************************************************************/

/* Allocates a contiguous ni x nj x nl array of cvector. The rows are carved
   out of a single block, so the number of heap allocations is 3 instead of
   ni*nj+ni+1 and the whole field is one contiguous, cache-friendly region. */

/* Rows are padded by FIELD_PAD elements. The lattice dimensions are powers of
   two, so an unpadded layout puts every row at an exact power-of-two offset and
   makes the three FFT sweeps hit the same cache sets over and over. A small pad
   breaks that aliasing at the cost of a tiny amount of memory. */

#ifndef FIELD_PAD
#define FIELD_PAD 1
#endif

static cvector *** alloc_field3(int ni, int nj, int nl)
{   size_t stride_l=(size_t)nl+FIELD_PAD;
    size_t plane=(size_t)nj*stride_l+FIELD_PAD;
    size_t ncell=(size_t)ni*plane;
    cvector *** p;
    cvector ** rows;
    cvector * blk;
    int i,j;

    p=(cvector ***)malloc((size_t)ni*sizeof(cvector **));
    rows=(cvector **)malloc((size_t)ni*(size_t)nj*sizeof(cvector *));
    blk=(cvector *)malloc(ncell*sizeof(cvector));
    if(p==NULL||rows==NULL||blk==NULL)
    {   printf("\n No enough memory for a %dx%dx%d field (%lu bytes)! ",
               ni,nj,nl,(unsigned long)(ncell*sizeof(cvector)));
        exit(1);
    }

    for(i=0;i<ni;i++)
    {   p[i]=rows+(size_t)i*nj;
        for(j=0;j<nj;j++)
            p[i][j]=blk+(size_t)i*plane+(size_t)j*stride_l;
    }

    return p;
}

static void free_field3(cvector *** p)
{   if(p==NULL)  return;
    free(p[0][0]);
    free(p[0]);
    free(p);
}

/****************************************************************************/
/************************ Workspace for the solver ***************************/
/****************************************************************************/

void dda_workspace_create(dda_workspace * ws, int i_max, int j_max, int l_max,
                          int iMode_SolutionMethod)
{   ws->i_max=i_max;    ws->j_max=j_max;    ws->l_max=l_max;

    ws->cvP=alloc_field3(i_max,j_max,l_max);
    ws->cvP_min=alloc_field3(i_max,j_max,l_max);
    ws->cvR=alloc_field3(i_max,j_max,l_max);
    ws->cvS=alloc_field3(i_max,j_max,l_max);

/* The convolution works on a lattice padded by a factor of two in each
   direction, hence the doubled dimensions here. */

    ws->cvWork=alloc_field3(2*i_max,2*j_max,2*l_max);

    if(iMode_SolutionMethod==2)
    {   ws->cvG=alloc_field3(i_max,j_max,l_max);
        ws->cvH=alloc_field3(i_max,j_max,l_max);
    }
    else
    {   ws->cvG=NULL;   ws->cvH=NULL;   }

    return;
}

void dda_workspace_free(dda_workspace * ws)
{   free_field3(ws->cvP);
    free_field3(ws->cvP_min);
    free_field3(ws->cvR);
    free_field3(ws->cvS);
    free_field3(ws->cvWork);
    free_field3(ws->cvG);
    free_field3(ws->cvH);

    ws->cvP=NULL;   ws->cvP_min=NULL;   ws->cvWork=NULL;
    ws->cvR=NULL;   ws->cvS=NULL;       ws->cvG=NULL;   ws->cvH=NULL;

    return;
}

/****************************************************************************/
/****************************** DDA Function ********************************/
/****************************************************************************/

int DDA_FFT(int iNoIllumination, int iMode_SolutionMethod, real rTau, int i_max,
            int j_max, int l_max, cplx * cListMaterials, int * iListMV,
            char *** cAlpha, int Max_Step, real Error_max, real d0, vector kinc,
            cvector Einc, cplx * cW_2_i, cplx * cW_2_j, cplx * cW_2_l, cplx * cW_i,
            cplx * cW_j, cplx * cW_l, real ** Im_cW_i, real ** Im_cW_j,
            real ** Im_cW_l, cplx ** cW2_i, cplx ** cW2_j, cplx ** cW2_l,
            int iMode_FFT, cmatrix *** cmA, int iNoScPl, int iNoA, vector ** ksc,
            cvector ** Esc, real * rCabs, real * rCext, dda_workspace * ws)
{   int i,j,l,o,m,step_min;
    real rA,rB,rrCabs,rrCext;
    long double Error_min,Error,RR,RR_1,SS;
    cplx cA,cB,cW,c0,c1,c2,c3;
    ld_cplx ldcA,ldcB,ldcW,ldcDenominator,ldcRho_1,ldcRho,ldc0,ldc1;
    cvector *** cvP,*** cvP_min,*** cvWork,*** cvR,*** cvS,*** cvG,*** cvH;
    clock_t time0,time1;

    if(ws==NULL||ws->i_max!=i_max||ws->j_max!=j_max||ws->l_max!=l_max)
    {   printf("\n The workspace does not match the lattice dimensions! ");
        exit(1);
    }

    cvP=ws->cvP;        cvP_min=ws->cvP_min;    cvWork=ws->cvWork;
    cvR=ws->cvR;        cvS=ws->cvS;
    cvG=ws->cvG;        cvH=ws->cvH;

    time0=time(NULL);

    // Zero iteration // cvP==cvP_inc
    FirstApproximation(i_max,j_max,l_max,iListMV,cListMaterials,cAlpha,d0,kinc,
                       Einc,cvP);
    CopyingBack(i_max,j_max,l_max,cvP,cvWork);

    if(iMode_FFT==0)
        Convolution((2*i_max),(2*j_max),(2*l_max),cW_2_i,cW_2_j,cW_2_l,
                    cW_i,cW_j,cW_l,cmA,cvWork);
    if(iMode_FFT==1)
        Convolution((2*i_max),(2*j_max),(2*l_max),Im_cW_i,Im_cW_j,Im_cW_l,
                    cW2_i,cW2_j,cW2_l,cmA,cvWork);

    for(i=0;i<i_max;i++)
        for(j=0;j<j_max;j++)
            for(l=0;l<l_max;l++)
            {   //cvR[i][j][l]=cvP_ind[i][j][l]-cListMaterials[(int)(cAlpha[i][j][l])]
                            //*cvWork[i][j][l]-cvP_inc[i][j][l];

                m=(int)cAlpha[i][j][l];
                if(iListMV[m]!=0)
                {   Cmul(cListMaterials[m],cvWork[i][j][l].x,c0);
                    Cmul(cListMaterials[m],cvWork[i][j][l].y,c1);
                    Cmul(cListMaterials[m],cvWork[i][j][l].z,c2);

                    cvR[i][j][l].x.r=-c0.r; cvR[i][j][l].x.i=-c0.i;
                    cvR[i][j][l].y.r=-c1.r; cvR[i][j][l].y.i=-c1.i;
                    cvR[i][j][l].z.r=-c2.r; cvR[i][j][l].z.i=-c2.i;
                }
                else
                {   cvR[i][j][l].x.r=0.;    cvR[i][j][l].x.i=0.;
                    cvR[i][j][l].y.r=0.;    cvR[i][j][l].y.i=0.;
                    cvR[i][j][l].z.r=0.;    cvR[i][j][l].z.i=0.;
                }
            }

    if(iMode_SolutionMethod==0)
        Copying(i_max,j_max,l_max,cvR,cvS);

    if(iMode_SolutionMethod==1)
    {   RR_1=ldScalMulConj(i_max,j_max,l_max,iListMV,cAlpha,cvR,cvR);

        CopyingBack(i_max,j_max,l_max,cvR,cvWork);

        if(iMode_FFT==0)
            Convolution_conj((2*i_max),(2*j_max),(2*l_max),cW_2_i,cW_2_j,
                            cW_2_l,cW_i,cW_j,cW_l,cmA,cvWork);
        if(iMode_FFT==1)
            Convolution_conj((2*i_max),(2*j_max),(2*l_max),Im_cW_i,Im_cW_j,
                            Im_cW_l,cW2_i,cW2_j,cW2_l,cmA,cvWork);

        for(i=0;i<i_max;i++)
            for(j=0;j<j_max;j++)
                for(l=0;l<l_max;l++)
                {   //cvS[i][j][l]=cvR[i][j][l]-cListMaterials[(int)(cAlpha[i][j][l])]
                                //*cvWork[i][j][l];

                    m=(int)cAlpha[i][j][l];
                    c3=conj(cListMaterials[m]);

                    if(iListMV[m]!=0)
                    {   Cmul(c3,cvWork[i][j][l].x,c0);
                        Cmul(c3,cvWork[i][j][l].y,c1);
                        Cmul(c3,cvWork[i][j][l].z,c2);

                        Csub(cvR[i][j][l].x,c0,cvS[i][j][l].x);
                        Csub(cvR[i][j][l].y,c1,cvS[i][j][l].y);
                        Csub(cvR[i][j][l].z,c2,cvS[i][j][l].z);
                    }
                    else
                    {   cvS[i][j][l].x.r=0.;    cvS[i][j][l].x.i=0.;
                        cvS[i][j][l].y.r=0.;    cvS[i][j][l].y.i=0.;
                        cvS[i][j][l].z.r=0.;    cvS[i][j][l].z.i=0.;
                    }
                }
    }

    if(iMode_SolutionMethod==2)
    {   Copying(i_max,j_max,l_max,cvR,cvG);
        Copying(i_max,j_max,l_max,cvR,cvS);

        ldcRho_1=ldcScalMul(i_max,j_max,l_max,iListMV,cAlpha,cvR,cvG);

        CopyingBack(i_max,j_max,l_max,cvS,cvWork);

        if(iMode_FFT==0)
            Convolution((2*i_max),(2*j_max),(2*l_max),cW_2_i,cW_2_j,
                        cW_2_l,cW_i,cW_j,cW_l,cmA,cvWork);
        if(iMode_FFT==1)
            Convolution((2*i_max),(2*j_max),(2*l_max),Im_cW_i,Im_cW_j,
                        Im_cW_l,cW2_i,cW2_j,cW2_l,cmA,cvWork);

        for(i=0;i<i_max;i++)
            for(j=0;j<j_max;j++)
                for(l=0;l<l_max;l++)
                {   //cvH[i][j][l]=cvS[i][j][l]-cListMaterials[(int)(cAlpha[i][j][l])]
                                    //*cvWork[i][j][l];

                    m=(int)cAlpha[i][j][l];
                    if(iListMV[m]!=0)
                    {   Cmul(cListMaterials[m],cvWork[i][j][l].x,c0);
                        Cmul(cListMaterials[m],cvWork[i][j][l].y,c1);
                        Cmul(cListMaterials[m],cvWork[i][j][l].z,c2);

                        Csub(cvS[i][j][l].x,c0,cvH[i][j][l].x);
                        Csub(cvS[i][j][l].y,c1,cvH[i][j][l].y);
                        Csub(cvS[i][j][l].z,c2,cvH[i][j][l].z);
                    }
                    else
                    {   cvH[i][j][l].x.r=0.;    cvH[i][j][l].x.i=0.;
                        cvH[i][j][l].y.r=0.;    cvH[i][j][l].y.i=0.;
                        cvH[i][j][l].z.r=0.;    cvH[i][j][l].z.i=0.;
                    }
                }

        ldc0=ldcScalMul(i_max,j_max,l_max,iListMV,cAlpha,cvG,cvH);
        ldcA=ldcRho_1/ldc0;

        cA.r=(real)ldcA.r;  cA.i=(real)ldcA.i;

        LinearComb(i_max,j_max,l_max,(-cA),iListMV,cAlpha,cvR,cvH,1);

        CopyingBack(i_max,j_max,l_max,cvR,cvWork);

        if(iMode_FFT==0)
            Convolution((2*i_max),(2*j_max),(2*l_max),cW_2_i,cW_2_j,
                        cW_2_l,cW_i,cW_j,cW_l,cmA,cvWork);
        if(iMode_FFT==1)
            Convolution((2*i_max),(2*j_max),(2*l_max),Im_cW_i,Im_cW_j,
                        Im_cW_l,cW2_i,cW2_j,cW2_l,cmA,cvWork);

        for(i=0;i<i_max;i++)
            for(j=0;j<j_max;j++)
                for(l=0;l<l_max;l++)
                {   //cvWork[i][j][l]=cvR[i][j][l]-cListMaterials[(int)(cAlpha[i][j][l])]
                                    //*cvWork[i][j][l];

                    m=(int)cAlpha[i][j][l];
                    if(iListMV[m]!=0)
                    {   Cmul(cListMaterials[m],cvWork[i][j][l].x,c0);
                        Cmul(cListMaterials[m],cvWork[i][j][l].y,c1);
                        Cmul(cListMaterials[m],cvWork[i][j][l].z,c2);

                        Csub(cvR[i][j][l].x,c0,cvWork[i][j][l].x);
                        Csub(cvR[i][j][l].y,c1,cvWork[i][j][l].y);
                        Csub(cvR[i][j][l].z,c2,cvWork[i][j][l].z);
                    }
                    else
                    {   cvWork[i][j][l].x.r=0.;     cvWork[i][j][l].x.i=0.;
                        cvWork[i][j][l].y.r=0.;     cvWork[i][j][l].y.i=0.;
                        cvWork[i][j][l].z.r=0.;     cvWork[i][j][l].z.i=0.;
                    }
                }


        ldc0=ldcScalMul(i_max,j_max,l_max,iListMV,cAlpha,cvWork,cvR);
        ldc1=ldcScalMul(i_max,j_max,l_max,iListMV,cAlpha,cvWork,cvWork);
        ldcW=ldc0/ldc1;

        cW.r=(real)ldcW.r;  cW.i=(real)ldcW.i;

        LinearComb(i_max,j_max,l_max,(-cA),(-cW),iListMV,cAlpha,cvP,cvS,cvR,cvP);
        LinearComb(i_max,j_max,l_max,(-cW),iListMV,cAlpha,cvR,cvWork,1);
    }

    if(iMode_SolutionMethod==3)
    {   Copying(i_max,j_max,l_max,cvP,cvS);
        
        for(i=0;i<i_max;i++)
            for(j=0;j<j_max;j++)
                for(l=0;l<l_max;l++)
                {   //cvP[i][j][l]=cvS[i][j][l]-cvR[i][j][l];

                    Csub(cvS[i][j][l].x,cvR[i][j][l].x,cvP[i][j][l].x);
                    Csub(cvS[i][j][l].y,cvR[i][j][l].y,cvP[i][j][l].y);
                    Csub(cvS[i][j][l].z,cvR[i][j][l].z,cvP[i][j][l].z);
                }
    }

    Error=ldError(i_max,j_max,l_max,iListMV,cAlpha,cvR);
    printf("\n Step    0    discrepancy=%9.7Lf",Error);
    Error_min=Error;    step_min=0;

/* The iterate of step 0 is the best one until some iteration improves on it,
   so it has to be stored right away. Otherwise, if the iterative process never
   improves (e.g. the simple iteration method with a too large weight factor
   diverges), cvP_min would still be uninitialised when the code below falls
   back to it. */

    for(i=0;i<i_max;i++)
        for(j=0;j<j_max;j++)
            for(l=0;l<l_max;l++)
                cvP_min[i][j][l]=cvP[i][j][l];

    time1=time(NULL);
    printf("\ttime=%ld sec.",(long)(time1-time0));
    fflush(stdout);

    for(o=1;o<=Max_Step;o++)
    {   if(Error<=Error_max)    break;
        time0=time1;

        if(Error<Error_min)
        {   Error_min=Error;
            step_min=o-1;

            for(i=0;i<i_max;i++)
                for(j=0;j<j_max;j++)
                    for(l=0;l<l_max;l++)
                        cvP_min[i][j][l]=cvP[i][j][l];
        }

        if(iMode_SolutionMethod==2)
        {   ldcRho=ldcScalMul(i_max,j_max,l_max,iListMV,cAlpha,cvR,cvG);
            ldcB=(ldcRho*ldcA)/(ldcRho_1*ldcW);

            cB.r=(real)ldcB.r;  cB.i=(real)ldcB.i;

            LinearComb(i_max,j_max,l_max,(cB),(-(cB*cW)),iListMV,cAlpha,
                       cvR,cvS,cvH,cvS);
        }

        if(iMode_SolutionMethod==0||iMode_SolutionMethod==1||iMode_SolutionMethod==2)
        {   CopyingBack(i_max,j_max,l_max,cvS,cvWork);

            if(iMode_FFT==0)
                Convolution((2*i_max),(2*j_max),(2*l_max),cW_2_i,cW_2_j,
                            cW_2_l,cW_i,cW_j,cW_l,cmA,cvWork);
            if(iMode_FFT==1)
                Convolution((2*i_max),(2*j_max),(2*l_max),Im_cW_i,Im_cW_j,
                            Im_cW_l,cW2_i,cW2_j,cW2_l,cmA,cvWork);

            for(i=0;i<i_max;i++)
                for(j=0;j<j_max;j++)
                    for(l=0;l<l_max;l++)
                    {   //cvWork[i][j][l]=cvS[i][j][l]-cListMaterials[(int)(cAlpha[i][j][l])]
                                        //*cvWork[i][j][l];

                        m=(int)cAlpha[i][j][l];
                        if(iListMV[m]!=0)
                        {   Cmul(cListMaterials[m],cvWork[i][j][l].x,c0);
                            Cmul(cListMaterials[m],cvWork[i][j][l].y,c1);
                            Cmul(cListMaterials[m],cvWork[i][j][l].z,c2);

                            Csub(cvS[i][j][l].x,c0,cvWork[i][j][l].x);
                            Csub(cvS[i][j][l].y,c1,cvWork[i][j][l].y);
                            Csub(cvS[i][j][l].z,c2,cvWork[i][j][l].z);
                        }
                        else
                        {   cvWork[i][j][l].x.r=0.;     cvWork[i][j][l].x.i=0.;
                            cvWork[i][j][l].y.r=0.;     cvWork[i][j][l].y.i=0.;
                            cvWork[i][j][l].z.r=0.;     cvWork[i][j][l].z.i=0.;
                        }
                    }
        }

        if(iMode_SolutionMethod==0)
        {   ldcDenominator=ldcScalMul(i_max,j_max,l_max,iListMV,cAlpha,cvS,cvWork);
            ldcA=ldcScalMul(i_max,j_max,l_max,iListMV,cAlpha,cvR,cvS)/ldcDenominator;
            cA.r=(real)ldcA.r;  cA.i=(real)ldcA.i;

            LinearComb(i_max,j_max,l_max,(-cA),iListMV,cAlpha,cvP,cvS,1);
            LinearComb(i_max,j_max,l_max,(-cA),iListMV,cAlpha,cvR,cvWork,1);

            ldcB=ldcScalMul(i_max,j_max,l_max,iListMV,cAlpha,cvR,cvWork)/ldcDenominator;
            cB.r=(real)ldcB.r;  cB.i=(real)ldcB.i;

            LinearComb(i_max,j_max,l_max,(-cB),iListMV,cAlpha,cvR,cvS,2);
        }

        if(iMode_SolutionMethod==1)
        {   SS=ldScalMulConj(i_max,j_max,l_max,iListMV,cAlpha,cvS,cvS);
            rA=(real)(-RR_1/SS);

            LinearComb(i_max,j_max,l_max,rA,iListMV,cAlpha,cvR,cvWork,1);
            LinearComb(i_max,j_max,l_max,rA,iListMV,cAlpha,cvP,cvS,1);

            RR=ldScalMulConj(i_max,j_max,l_max,iListMV,cAlpha,cvR,cvR);
            rB=(real)(RR/RR_1); RR_1=RR;

            CopyingBack(i_max,j_max,l_max,cvR,cvWork);

            if(iMode_FFT==0)
                Convolution_conj((2*i_max),(2*j_max),(2*l_max),cW_2_i,cW_2_j,
                                cW_2_l,cW_i,cW_j,cW_l,cmA,cvWork);
            if(iMode_FFT==1)
                Convolution_conj((2*i_max),(2*j_max),(2*l_max),Im_cW_i,Im_cW_j,
                                Im_cW_l,cW2_i,cW2_j,cW2_l,cmA,cvWork);

            for(i=0;i<i_max;i++)
                for(j=0;j<j_max;j++)
                    for(l=0;l<l_max;l++)
                    {   //cvWork[i][j][l]=cvR[i][j][l]-cListMaterials[(int)(cAlpha[i][j][l])]
                                        //*cvWork[i][j][l];

                        m=(int)cAlpha[i][j][l];
                        c3=conj(cListMaterials[m]);

                        if(iListMV[m]!=0)
                        {   Cmul(c3,cvWork[i][j][l].x,c0);
                            Cmul(c3,cvWork[i][j][l].y,c1);
                            Cmul(c3,cvWork[i][j][l].z,c2);

                            Csub(cvR[i][j][l].x,c0,cvWork[i][j][l].x);
                            Csub(cvR[i][j][l].y,c1,cvWork[i][j][l].y);
                            Csub(cvR[i][j][l].z,c2,cvWork[i][j][l].z);
                        }
                        else
                        {   cvWork[i][j][l].x.r=0.;     cvWork[i][j][l].x.i=0.;
                            cvWork[i][j][l].y.r=0.;     cvWork[i][j][l].y.i=0.;
                            cvWork[i][j][l].z.r=0.;     cvWork[i][j][l].z.i=0.;
                        }
                    }

            LinearComb(i_max,j_max,l_max,rB,iListMV,cAlpha,cvWork,cvS,2);
        }

        if(iMode_SolutionMethod==2)
        {   Copying(i_max,j_max,l_max,cvWork,cvH);
            ldc0=ldcScalMul(i_max,j_max,l_max,iListMV,cAlpha,cvG,cvH);
            ldcA=ldcRho/ldc0;
            cA.r=(real)ldcA.r;  cA.i=(real)ldcA.i;

            LinearComb(i_max,j_max,l_max,(-cA),iListMV,cAlpha,cvR,cvH,1);

            CopyingBack(i_max,j_max,l_max,cvR,cvWork);

            if(iMode_FFT==0)
                Convolution((2*i_max),(2*j_max),(2*l_max),cW_2_i,cW_2_j,
                            cW_2_l,cW_i,cW_j,cW_l,cmA,cvWork);
            if(iMode_FFT==1)
                Convolution((2*i_max),(2*j_max),(2*l_max),Im_cW_i,Im_cW_j,
                            Im_cW_l,cW2_i,cW2_j,cW2_l,cmA,cvWork);

            for(i=0;i<i_max;i++)
                for(j=0;j<j_max;j++)
                    for(l=0;l<l_max;l++)
                    {   //cvWork[i][j][l]=cvR[i][j][l]-cListMaterials[(int)(cAlpha[i][j][l])]
                                        //*cvWork[i][j][l];

                        m=(int)cAlpha[i][j][l];
                        if(iListMV[m]!=0)
                        {   Cmul(cListMaterials[m],cvWork[i][j][l].x,c0);
                            Cmul(cListMaterials[m],cvWork[i][j][l].y,c1);
                            Cmul(cListMaterials[m],cvWork[i][j][l].z,c2);

                            Csub(cvR[i][j][l].x,c0,cvWork[i][j][l].x);
                            Csub(cvR[i][j][l].y,c1,cvWork[i][j][l].y);
                            Csub(cvR[i][j][l].z,c2,cvWork[i][j][l].z);
                        }
                        else
                        {   cvWork[i][j][l].x.r=0.;     cvWork[i][j][l].x.i=0.;
                            cvWork[i][j][l].y.r=0.;     cvWork[i][j][l].y.i=0.;
                            cvWork[i][j][l].z.r=0.;     cvWork[i][j][l].z.i=0.;
                        }
                    }


            ldc0=ldcScalMul(i_max,j_max,l_max,iListMV,cAlpha,cvWork,cvR);
            ldc1=ldcScalMul(i_max,j_max,l_max,iListMV,cAlpha,cvWork,cvWork);
            ldcW=ldc0/ldc1;
            cW.r=(real)ldcW.r;  cW.i=(real)ldcW.i;

            LinearComb(i_max,j_max,l_max,(-cA),(-cW),iListMV,cAlpha,cvP,cvS,cvR,cvP);
            LinearComb(i_max,j_max,l_max,(-cW),iListMV,cAlpha,cvR,cvWork,1);
            ldcRho_1=ldcRho;
        }

        if(iMode_SolutionMethod==3)
        {   CopyingBack(i_max,j_max,l_max,cvP,cvWork);

            if(iMode_FFT==0)
                Convolution((2*i_max),(2*j_max),(2*l_max),cW_2_i,cW_2_j,
                            cW_2_l,cW_i,cW_j,cW_l,cmA,cvWork);
            if(iMode_FFT==1)
                Convolution((2*i_max),(2*j_max),(2*l_max),Im_cW_i,Im_cW_j,
                            Im_cW_l,cW2_i,cW2_j,cW2_l,cmA,cvWork);

            for(i=0;i<i_max;i++)
                for(j=0;j<j_max;j++)
                    for(l=0;l<l_max;l++)
                    {   m=(int)cAlpha[i][j][l];

                        if(iListMV[m]!=0)
                        {   //cvWork[i][j][l]=cvP[i][j][l]-cListMaterials[(int)(cAlpha[i][j][l])]
                                            //*cvWork[i][j][l];
                            Cmul(cListMaterials[m],cvWork[i][j][l].x,c0);
                            Cmul(cListMaterials[m],cvWork[i][j][l].y,c1);
                            Cmul(cListMaterials[m],cvWork[i][j][l].z,c2);

                            Csub(cvP[i][j][l].x,c0,cvWork[i][j][l].x);
                            Csub(cvP[i][j][l].y,c1,cvWork[i][j][l].y);
                            Csub(cvP[i][j][l].z,c2,cvWork[i][j][l].z);

                            //cvR[i][j][l]=cvWork[i][j][l]-cvS[i][j][l];
                            Csub(cvWork[i][j][l].x,cvS[i][j][l].x,cvR[i][j][l].x);
                            Csub(cvWork[i][j][l].y,cvS[i][j][l].y,cvR[i][j][l].y);
                            Csub(cvWork[i][j][l].z,cvS[i][j][l].z,cvR[i][j][l].z);

                            //cvP_next[i][j][l]=cvP[i][j][l]-rTau*cvR[i][j][l];
                            DCmul(rTau,cvR[i][j][l].x,c0);
                            DCmul(rTau,cvR[i][j][l].y,c1);
                            DCmul(rTau,cvR[i][j][l].z,c2);

                            Csub(cvP[i][j][l].x,c0,cvWork[i][j][l].x);
                            Csub(cvP[i][j][l].y,c1,cvWork[i][j][l].y);
                            Csub(cvP[i][j][l].z,c2,cvWork[i][j][l].z);

                            cvP[i][j][l]=cvWork[i][j][l];
                        }
                        else
                        {   cvP[i][j][l].x.r=0.;    cvP[i][j][l].x.i=0.;
                            cvP[i][j][l].y.r=0.;    cvP[i][j][l].y.i=0.;
                            cvP[i][j][l].z.r=0.;    cvP[i][j][l].z.i=0.;

                            cvR[i][j][l]=cvP[i][j][l];
                        }
                    }
        }

        Error=ldError(i_max,j_max,l_max,iListMV,cAlpha,cvR);
        printf("\n Step %4d    discrepancy=%9.7Lf",o,Error);
        time1=time(NULL);
        printf("\ttime=%ld sec.",(long)(time1-time0));
        if(div(o,10).rem==0)
            printf("  irradiation #%i",(iNoIllumination+1));

        fflush(stdout);
    }

    if(Error>Error_min)
    {   for(i=0;i<i_max;i++)
            for(j=0;j<j_max;j++)
                for(l=0;l<l_max;l++)
                    cvP[i][j][l]=cvP_min[i][j][l];

        printf("\n Result of the last iteration has been replaced with the most optimal!");
        printf("\n Step %3d    discrepancy=%9.7Lf",step_min,Error_min);
        fflush(stdout);
    }

    time0=time(NULL);

    CabsCext(d0,kinc,Einc,i_max,j_max,l_max,cListMaterials,cAlpha,cvP,&rrCabs,
            &rrCext);

    (* rCabs)=rrCabs;   (* rCext)=rrCext;

    for(i=0;i<iNoScPl;i++)
        for(j=0;j<iNoA;j++)
            Esc[i][j]=FarField(i_max,j_max,l_max,iListMV,cAlpha,d0,cvP,ksc[i][j]);

    time1=time(NULL);
    printf("\n Calculations in %d scattering planes have been done for ",iNoScPl);
    printf("%ld sec.",(long)(time1-time0));
    fflush(stdout);


    return 0;
}

/****************************************************************************/
/********************* Calculation field in the far zone ********************/
/****************************************************************************/

cvector FarField(int i_max, int j_max, int l_max, int * iListMV,
                 char *** cAlpha, real d0, cvector *** cvP, vector ksc_ort)
{   int i,j,l,m;
    real r0;
    cplx c0,c1,c2,c3,c4,c5,c6,c7;
    vector r,vCenter={(real)((i_max-1)*d0/2.),(real)((j_max-1)*d0/2.),(real)((l_max-1)*d0/2.)};
    cvector cvTmp={{0.,0.},{0.,0.},{0.,0.}};

    for(i=0;i<i_max;i++)
    {   r.x=(real)i*d0-vCenter.x;
        for(j=0;j<j_max;j++)
        {   r.y=(real)j*d0-vCenter.y;
            for(l=0;l<l_max;l++)
            {
                m=(int)cAlpha[i][j][l];
                if(iListMV[m]!=0)
                {   r.z=(real)l*d0-vCenter.z;

                    r0=ksc_ort.x*r.x+ksc_ort.y*r.y+ksc_ort.z*r.z;
                    Cexp_im(r0,c0);

                    //cvTmp=cvTmp+(cvP[i][j][l]-ksc_ort*(ksc_ort*cvP[i][j][l]))*c0;

                    DCmul(ksc_ort.x,cvP[i][j][l].x,c1);
                    DCmul(ksc_ort.y,cvP[i][j][l].y,c2);
                    DCmul(ksc_ort.z,cvP[i][j][l].z,c3);
                    Cadd(c1,c2,c4); Cadd(c3,c4,c1);

                    DCmul(ksc_ort.x,c1,c2);
                    DCmul(ksc_ort.y,c1,c3);
                    DCmul(ksc_ort.z,c1,c4);

                    Csub(cvP[i][j][l].x,c2,c5);
                    Csub(cvP[i][j][l].y,c3,c6);
                    Csub(cvP[i][j][l].z,c4,c7);

                    Cmul(c0,c5,c2);
                    Cmul(c0,c6,c3);
                    Cmul(c0,c7,c4);

                    Cadd(cvTmp.x,c2,cvTmp.x);
                    Cadd(cvTmp.y,c3,cvTmp.y);
                    Cadd(cvTmp.z,c4,cvTmp.z);
                }
            }
        }
    }

    return cvTmp;
}

/****************************************************************************/
/*************************** Calculation of alpha ***************************/
/****************************************************************************/

cplx Polarizability_CM(real d0, cplx cRefrIndex)
{   cplx cTmp;

    cTmp=(d0*d0*d0)*(0.75/PI)*(cRefrIndex*cRefrIndex-1)/(cRefrIndex*cRefrIndex+2);
    return cTmp;
}

cplx Polarizability_CMRR(cplx alpha_CM)
{   cplx c0={0.,1.};
    cplx cTmp;

    cTmp=alpha_CM/(1.-((2./3.)*c0*alpha_CM));
    return cTmp;
}

cplx Polarizability_CMNR(real d0, cplx alpha_CM)
{   real b=-sqrt(4.*PI/3.)/d0;
    cplx c0={0.,1.},c1;
    cplx cTmp;

    c1=alpha_CM/(1.+alpha_CM*b);
    cTmp=c1/(1.-((2./3.)*c0*c1));

    return cTmp;
}

cplx Polarizability_ILDR(real d0, cplx cRefrIndex, cplx alpha_CM)
{   real b1=-1.8915316,b2=0.1648469,b3=-1.7700004,S=.2;
    cplx epsilon=cRefrIndex*cRefrIndex,c0={0.,1.};
    cplx cTmp;

    cTmp=alpha_CM/(1+(alpha_CM/d0)*(b1+b2*epsilon+S*b3*epsilon-(2./3.)*d0*c0));
    return cTmp;
}

cplx Polarizability_LDR(real d0, cplx cRefrIndex, cplx alpha_CM,
                        vector kinc_ort, cvector Einc_ort)
{   real b1=-1.8915316,b2=0.1648469,b3=-1.7700004,S;
    real rx=kinc_ort.x*Einc_ort.x.r,ry=kinc_ort.y*Einc_ort.y.r,rz=kinc_ort.z*Einc_ort.z.r;
    cplx epsilon=cRefrIndex*cRefrIndex,c0={0.,1.};
    cplx cTmp;

    S=rx*rx+ry*ry+rz*rz;

    cTmp=alpha_CM/(1+(alpha_CM/d0)*(b1+b2*epsilon+S*b3*epsilon-(2./3.)*d0*c0));
    return cTmp;
}

void CabsCext(real d0, vector kinc_ort, cvector Einc, int i_max, int j_max,
              int l_max, cplx * cListMaterials, char *** cAlpha,
              cvector *** cvP, real * Cabs, real * Cext)
{   int i,j,l,m;
    cplx c0,c1,c2,c3,c4,c5;
    real r0=((4.*PI)/re(Einc*conj(Einc))),r1=0.,r2=0.,r3,r4=2./3.;
    vector r,vCenter={(real)(i_max*d0/2.),(real)(j_max*d0/2.),(real)(l_max*d0/2.)};
    cvector cv0=conj(Einc);

    for(i=0;i<i_max;i++)
    {   r.x=(real)i*d0;
        for(j=0;j<j_max;j++)
        {   r.y=(real)j*d0;
            for(l=0;l<l_max;l++)
            {   r.z=(real)l*d0;

                m=(int)(cAlpha[i][j][l]);
                if(cListMaterials[m].r==0.&&cListMaterials[m].i==0.);
                else
                {   r3=-(kinc_ort.x*(r.x-vCenter.x)
                        +kinc_ort.y*(r.y-vCenter.y)
                        +kinc_ort.z*(r.z-vCenter.z));
                    c0=exp_i(r3);

                    Cmul(cvP[i][j][l].x,cv0.x,c1);
                    Cmul(cvP[i][j][l].y,cv0.y,c2);
                    Cmul(cvP[i][j][l].z,cv0.z,c3);
                    Cadd(c1,c2,c4); Cadd(c3,c4,c5);

                    Cmul(c0,c5,c1);
                    r1+=c1.i;

                    Cmuls(cvP[i][j][l].x,cvP[i][j][l].x,c1);
                    Cmuls(cvP[i][j][l].y,cvP[i][j][l].y,c2);
                    Cmuls(cvP[i][j][l].z,cvP[i][j][l].z,c3);
                    Cadd(c1,c2,c4); Cadd(c3,c4,c0);

                    c1.r=cListMaterials[m].r;   c1.i=-cListMaterials[m].i;
                    Cdiv(c0,c1,c2);

                    r2+=(c2.i-r4*c0.r);
                }
            }
        }
    }

    (* Cext)=r0*r1; (* Cabs)=r0*r2;
    return;
}
