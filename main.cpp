
//Time dependence is exp(-iwt)
//Refractive index is n=n'+in''

#include "function.hpp"
#include "randomize.hpp"
#include "options.hpp"

#include <thread>
#include <vector>
#include <mutex>
#include <atomic>
#include <unistd.h>


/****************************************************************************/
/****************************** Parallel workers ****************************/
/****************************************************************************/

/* Amount of physical memory in bytes, 0 if it cannot be determined. */

static double physical_memory(void)
{   long pages=sysconf(_SC_PHYS_PAGES);
    long size=sysconf(_SC_PAGE_SIZE);

    if(pages<=0||size<=0)   return 0.;
    return (double)pages*(double)size;
}

/* Bytes that one worker needs on top of the shared A-matrix: the four field
   arrays, the padded convolution buffer, the two BiCGStab arrays and the
   occupancy map. */

/* The count requested by --threads, or DDA_THREADS from the environment if
   the option is absent, one worker per core otherwise, but never more workers
   than the private workspaces allow.

   One worker needs the six i*j*l field arrays (cvP, cvP_min, cvR, cvS, cvG,
   cvH), the padded 2i x 2j x 2l convolution buffer and the occupancy map;
   the A-matrix is shared by all of them. */

static int worker_count(int i_max, int j_max, int l_max, int iRequested)
{   const char * s=getenv("DDA_THREADS");
    double shared,footprint,budget;
    int n,fits;

    n=(int)sysconf(_SC_NPROCESSORS_ONLN);
    if(n<1) n=1;

    if(iRequested>0)    n=iRequested;
    else if(s!=NULL)
    {   int v=atoi(s);
        if(v>0) n=v;
    }

    {   double cells=(double)i_max*(double)j_max*(double)l_max;

        footprint=6.*cells*sizeof(cvector)+8.*cells*sizeof(cvector)+cells;
        shared=8.*cells*sizeof(cmatrix);
    }

    budget=0.7*physical_memory()-shared;
    if(budget<=0.||footprint<=0.)   return n;

    fits=(int)(budget/footprint);
    if(fits<1)  fits=1;
    if(fits<n)
    {   printf("\n %d worker(s) fit into memory, %.0f MB each; using %d.",
               fits,footprint/1048576.,fits);
        printf("\n Use --threads to change this. ");
        fflush(stdout);
        n=fits;
    }

    return n;
}

/* Contiguous storage for the occupancy map of one worker. */

static char *** alloc_char3(int ni, int nj, int nl)
{   size_t ncell=(size_t)ni*(size_t)nj*(size_t)nl;
    char *** p=(char ***)malloc((size_t)ni*sizeof(char **));
    char ** rows=(char **)malloc((size_t)ni*(size_t)nj*sizeof(char *));
    char * blk=(char *)malloc(ncell*sizeof(char));
    int i,j;

    if(p==NULL||rows==NULL||blk==NULL)
    {   printf("\n No enough memory for the occupancy map! ");
        exit(1);
    }

    for(i=0;i<ni;i++)
    {   p[i]=rows+(size_t)i*nj;
        for(j=0;j<nj;j++)
            p[i][j]=blk+((size_t)i*nj+j)*nl;
    }

    return p;
}

static void free_char3(char *** p)
{   if(p==NULL)  return;
    free(p[0][0]);
    free(p[0]);
    free(p);
}

/* Everything that a worker owns privately. The A-matrix, the FFT tables and
   the reference geometry are shared between the workers and only read. */

typedef
    struct
    {   char *** cAlpha;             /* occupancy map of the current sample    */
        cplx * cListMaterials;       /* polarizabilities, may vary by sample   */
        vector * vNormal;            /* normals of the scattering planes       */
        vector ** ksc;               /* scattering directions                  */
        vector kinc;                 /* incident wave vector                   */
        cvector Einc[2];             /* incident field, two polarizations      */
        cvector ** Esc[2];           /* scattered field, two polarizations    */
        cylinder * worms;            /* geometry of the worm-like samples      */
        dda_workspace ws;            /* scratch arrays of the solver           */
        mmuller * mmAcc;             /* Muller matrix of this sample           */
        mmuller * mm0;               /* scratch for the averaging              */
    }
    dda_worker;


int main(int argc, char ** argv)
{   time_t time_begin,time_end;

    int iRealizationNumber,iRealizationCounter,iOrientationMode,iWormsNumber=-1;
    int iMode_FFT,iMode_SolutionMethod,iHistoryCondition=0;
    int iHierarchicalLevels,iChoiceStructure,iChoicePolarizability,iNoScPl;
    int i_max,j_max,l_max,X_min,X_max,Y_min,Y_max,Z_min,Z_max,Max_Step;

    real x,Error_max,psi,iota,kappa,theta_inf,theta_up,theta_step,b_a,c_a,rVolume;
    real rho,radius_worm,length_worm,lr_worm,Cabs,Cext,rCabs,rCext;
    cylinder * worms;
    cplx * cRefrIndex, * cListMaterials;

/******************************** Other variable ****************************/

    int i,j,l,m,iNoA,iNoD,iSorts,iBites,X,Y,Z,T;
    int * iSeeds,* iListMV;
    char * filename1,* filename2;
    real kd,theta,cos_psi,sin_psi,cos_iota,sin_iota,cos_kappa,sin_kappa,rDepth;
    real r0,r1,R1,R2,R1_R2,rTau;
    real VectTrans[3][3];
    cplx alpha;
    real ** Im_cW_i,** Im_cW_j,** Im_cW_l;
    cplx * cW_2_i,* cW_2_j,* cW_2_l,* cW_i,* cW_j,* cW_l;
    cplx ** cW2_i,** cW2_j,** cW2_l;
    vector kinc_ini={0.,0.,1.},kinc;
    cvector Einc_ini[2]={{{1.,0.},{0.,0.},{0.,0.}},{{0.,0.},{1.,0.},{0.,0.}}},Einc[2];
    mjones  mjScat;
    mmuller * mmAccumulator,* mm0;

    cvector ** Esc[2];
    vector * vNormal_ini,* vNormal;
    vector ** ksc_ini,** ksc;

    char *** cAlpha;
    cmatrix *** cmA;
    dda_worker * vWorkers;
    int iWorkers,iFirstRealization,iTasks;
    real psiFixed,iotaFixed,kappaFixed;
    FILE * output_data,* storage_data;

    printf("\nComputation of light scattering by arbitrary target:");
    printf("\nDiscrete Dipole Approximation (DDA) with Fast Fourier Transform (FFT)");
    printf("\nVersion 6.3, Copyright (C) 2001-2012 E.S.Zubko");
    printf("\nVersion 7.2, Modification and optimization: K.A.Shmirko (kshmirko@dvo.ru) 2026\n");
    fflush(stdout);

    filename1=(char *)malloc(1024*sizeof(char));
    filename2=(char *)malloc(1024*sizeof(char));

/****************************************************************************/
/******************************** Input data ********************************/
/****************************************************************************/

    sim_options opt;

    {   int parse_rc=parse_options(argc,argv,&opt);

        if(parse_rc==PARSE_EXIT_OK) return 0;
        if(parse_rc!=PARSE_OK)      return 2;
    }

    if(opt.has_seed)    random_init_seed(opt.seed);
    else                random_init();
    time_begin=time(NULL);

    theta_inf=(real)opt.theta_min;
    theta_up=(real)opt.theta_max;
    theta_step=(real)opt.theta_step;

    iNoA=rounding(((theta_up-theta_inf)/theta_step)+1.);

    mmAccumulator=(mmuller *)malloc(iNoA*sizeof(mmuller));
    if(mmAccumulator==NULL)
    {   printf("\n No memory for array mmAccumulator! ");
        exit(1);
    }

    mm0=(mmuller *)malloc(iNoA*sizeof(mmuller));
    if(mm0==NULL)
    {   printf("\n No memory for array mm0! ");
        exit(1);
    }

    strcpy(filename1,opt.output_file.c_str());

    // begin of reading of the dst-file
    if(opt.restart)     remove(opt.state_file.c_str());

    storage_data=fopen(opt.state_file.c_str(),"rb+");
    if(storage_data==NULL)          // file doesn't exist; it needs to be created
    {   storage_data=fopen(opt.state_file.c_str(),"wb+");
        if(storage_data==NULL)
        {   printf("\n Cannot create the dst-file! ");
            exit(2);
        }

        iRealizationCounter=0;
        l=fwrite(&iRealizationCounter,sizeof(iRealizationCounter),1,storage_data);
        Cabs=0.;
        l=fwrite(&Cabs,sizeof(Cabs),1,storage_data);
        Cext=0.;
        l=fwrite(&Cext,sizeof(Cext),1,storage_data);

        for(i=0;i<iNoA;i++) mmAccumulator[i]=mmNull();

        for(i=0;i<iNoA;i++)
            l=fwrite(&mmAccumulator[i],sizeof(mmAccumulator[i]),1,storage_data);

        fflush(storage_data);
    }
    fclose(storage_data);   // closing the dst-file

    storage_data=fopen(opt.state_file.c_str(),"rb+");
    if(storage_data!=NULL)
    {   l=fread(&iRealizationCounter,sizeof(iRealizationCounter),1,storage_data);
        if(l!=1)
        {   printf("\n Cannot read data from the dst-file !");
            exit(4);
        }

        l=fread(&Cabs,sizeof(Cabs),1,storage_data);
        if(l!=1)
        {   printf("\n Cannot read data from the dst-file! ");
            exit(4);
        }

        l=fread(&Cext,sizeof(Cext),1,storage_data);
        if(l!=1)
        {   printf("\n Cannot read data from the dst-file! ");
            exit(4);
        }

        for(i=0;i<iNoA;i++)
        {   l=fread(&mmAccumulator[i],sizeof(mmAccumulator[i]),1,storage_data);
            if(l!=1)
            {   printf("\n Cannot read data from the dst-file! ");
                exit(4);
            }
        }
    }
    else
    {   printf("\n Some trouble has happened with dst-file.");
        printf("\n The program will be terminated. Please, find out what is going on.");
        exit(4);
    }
    // end of reading of the dst-file

    iRealizationNumber=opt.samples;
    iOrientationMode=opt.orientation;

    if(iOrientationMode==0)
    {   psi=(real)opt.psi*PI/180.;
        iota=(real)opt.iota*PI/180.;
        kappa=(real)opt.kappa*PI/180.;
    }

    iNoScPl=opt.planes;

    iChoiceStructure=opt.target;

    iSorts=(int)opt.materials.size();
    cRefrIndex=opt.materials.empty()? NULL : &opt.materials[0];
    iSeeds=opt.seeds.empty()? NULL : &opt.seeds[0];

    b_a=(real)opt.aspect_b;
    c_a=(real)opt.aspect_c;
    R1_R2=(real)opt.radius_ratio;
    iHierarchicalLevels=opt.hierarchical_levels;
    iBites=opt.bites;
    rDepth=(real)opt.depth;
    rho=(real)opt.packing;
    radius_worm=(real)opt.monomer_radius;
    lr_worm=(real)opt.monomer_length_ratio;

    if(iChoiceStructure==11)    strcpy(filename2,opt.shape_file.c_str());

    i_max=opt.grid[0];
    j_max=opt.grid[1];
    l_max=opt.grid[2];

    x=(real)opt.size_parameter;
    iChoicePolarizability=opt.polarizability;
    Error_max=(real)opt.tolerance;
    Max_Step=opt.max_steps;
    iMode_FFT=opt.fft_mode;
    iMode_SolutionMethod=opt.solver;
    rTau=(real)opt.tau;


    kd=2.*x/(real)i_max;

    if(iRealizationCounter>=iRealizationNumber)
        goto NO_CALC;

    printf("\n Total number of cells in the cubic lattice = %d",i_max*j_max*l_max);
    fflush(stdout);

/****************************************************************************/
/**************************** Memory allocation *****************************/
/****************************************************************************/

    vNormal_ini=(vector *)malloc(iNoScPl*sizeof(vector));
    vNormal=(vector *)malloc(iNoScPl*sizeof(vector));
    if(vNormal_ini==NULL||vNormal==NULL)
    {   printf("\n No memory for array vNormal_ini or vNormal!");
        exit(1);
    }
    vNormal_ini[0].x=0.;    vNormal_ini[0].y=1.;    vNormal_ini[0].z=0.;

    ksc_ini=(vector **)malloc(iNoScPl*sizeof(vector *));
    ksc=(vector **)malloc(iNoScPl*sizeof(vector *));
    if(ksc_ini==NULL||ksc==NULL)
    {   printf("\n No memory for array ksc_ini or ksc!");
        exit(1);
    }
    for(m=0;m<iNoScPl;m++)
    {   ksc_ini[m]=(vector *)malloc(iNoA*sizeof(vector));
        ksc[m]=(vector *)malloc(iNoA*sizeof(vector));

        if(ksc_ini[m]==NULL||ksc[m]==NULL)
        {   printf("\n No memory for array ksc_ini[%d] or ksc[%d]!",m,m);
            exit(1);
        }
    }

    for(m=0;m<iNoA;m++)
    {   theta=theta_inf+theta_step*m;
        SphCart_d(1.,0.,theta,ksc_ini[0][m]);
    }

    for(i=0;i<2;i++)
    {   Esc[i]=(cvector **)malloc(iNoScPl*sizeof(cvector *));
        if(Esc[i]==NULL)
        {   printf("\n No memory for array Esc[%d]!",i);
            exit(1);
        }

        for(m=0;m<iNoScPl;m++)
        {   Esc[i][m]=(cvector *)malloc(iNoA*sizeof(cvector));

            if(Esc[i][m]==NULL)
            {   printf("\n No memory for array Esc[%d][%d]!",i,m);
                exit(1);
            }
        }
    }

    // Begin of definition of initial collection of scattering plane
    r0=180./((real)iNoScPl);    r1=0.;
    for(m=1;m<iNoScPl;m++)
    {   r1+=r0;
        cos_psi=cos(r1*(PI/180.));  sin_psi=sin(r1*(PI/180.));

        vNormal_ini[m]=Rotation_z(cos_psi,sin_psi,vNormal_ini[0]);

        for(i=0;i<iNoA;i++)
            ksc_ini[m][i]=Rotation_z(cos_psi,sin_psi,ksc_ini[0][i]);

    }
    // End of definition of initial collection of scattering plane

    cmA=(cmatrix ***)malloc(2*i_max*sizeof(cmatrix **));
    if(cmA==NULL)
    {   printf("\n No enough memory for array cmA! ");
        exit(1);
    }
    for(i=0;i<2*i_max;i++)
    {   cmA[i]=(cmatrix **)malloc(2*j_max*sizeof(cmatrix *));
        if(cmA[i]==NULL)
        {   printf("\n No enough memory for array cmA[%d]! ",i);
            exit(1);
        }
    }
    for(i=0;i<2*i_max;i++)
        for(j=0;j<2*j_max;j++)
        {   cmA[i][j]=(cmatrix *)malloc(2*l_max*sizeof(cmatrix));
            if(cmA[i][j]==NULL)
            {   printf("\n No enough memory for array cmA[%d][%d]! ",i,j);
                exit(1);
            }
        }

    if(iMode_FFT==0)
    {   cW_2_i=(cplx *)malloc((i_max/2)*sizeof(cplx));
        cW_2_j=(cplx *)malloc((j_max/2)*sizeof(cplx));
        cW_2_l=(cplx *)malloc((l_max/2)*sizeof(cplx));
        cW_i=(cplx *)malloc((i_max)*sizeof(cplx));
        cW_j=(cplx *)malloc((j_max)*sizeof(cplx));
        cW_l=(cplx *)malloc((l_max)*sizeof(cplx));
        if(cW_2_i==NULL||cW_2_j==NULL||cW_2_l==NULL||cW_i==NULL||cW_j==NULL||cW_l==NULL)
        {   printf("\n No enough memory for one of cW arrays! ");
            exit(1);
        }
    }

    if(iMode_FFT==1)
    {   cW2_i=(cplx **)malloc((i_max)*sizeof(cplx *));
        cW2_j=(cplx **)malloc((j_max)*sizeof(cplx *));
        cW2_l=(cplx **)malloc((l_max)*sizeof(cplx *));
        Im_cW_i=(real **)malloc((i_max)*sizeof(real *));
        Im_cW_j=(real **)malloc((j_max)*sizeof(real *));
        Im_cW_l=(real **)malloc((l_max)*sizeof(real *));

        if(cW2_i==NULL||cW2_j==NULL||cW2_l==NULL||
           Im_cW_i==NULL||Im_cW_j==NULL||Im_cW_l==NULL)
        {   printf("\n No enough memory for one of cW arrays! ");
            exit(1);
        }

        for(i=0;i<i_max;i++)
        {   cW2_i[i]=(cplx *)malloc((i_max)*sizeof(cplx));
            Im_cW_i[i]=(real *)malloc((i_max)*sizeof(real));
            if(cW2_i[i]==NULL||Im_cW_i[i]==NULL)
            {   printf("\n No enough memory for cW_i or Im_cW_i array!");
                exit(1);
            }
        }

        for(j=0;j<j_max;j++)
        {   cW2_j[j]=(cplx *)malloc((j_max)*sizeof(cplx));
            Im_cW_j[j]=(real *)malloc((j_max)*sizeof(real));
            if(cW2_j[j]==NULL||Im_cW_j[j]==NULL)
            {   printf("\n No enough memory for cW_j or Im_cW_j array!");
                exit(1);
            }
        }

        for(l=0;l<l_max;l++)
        {   cW2_l[l]=(cplx *)malloc((l_max)*sizeof(cplx));
            Im_cW_l[l]=(real *)malloc((l_max)*sizeof(real));
            if(cW2_l[l]==NULL||Im_cW_l[l]==NULL)
            {   printf("\n No enough memory for cW_l or Im_cW_l array!");
                exit(1);
            }
        }
    }

    CreateAMatrix(i_max,j_max,l_max,kd,cmA);

    if(iMode_FFT==0)
    {   FourierPhase((i_max/2),cW_2_i);
        FourierPhase((j_max/2),cW_2_j);
        FourierPhase((l_max/2),cW_2_l);
        FourierPhase(i_max,cW_i);
        FourierPhase(j_max,cW_j);
        FourierPhase(l_max,cW_l);

        FFT_3D((2*i_max),(2*j_max),(2*l_max),cW_i,cW_j,cW_l,cmA);
    }

    if(iMode_FFT==1)
    {   FourierPhase(2*i_max,Im_cW_i,cW2_i);
        FourierPhase(2*j_max,Im_cW_j,cW2_j);
        FourierPhase(2*l_max,Im_cW_l,cW2_l);

        FFT_3D((2*i_max),(2*j_max),(2*l_max),Im_cW_i,Im_cW_j,Im_cW_l,
                cW2_i,cW2_j,cW2_l,cmA);
    }


    cAlpha=(char ***)malloc(i_max*sizeof(char **));
    if(cAlpha==NULL)
    {   printf("\n No enough memory for array cAlpha! ");
        exit(1);
    }
    for(i=0;i<i_max;i++)
    {   cAlpha[i]=(char **)malloc(j_max*sizeof(char *));
        if(cAlpha[i]==NULL)
        {   printf("\n No enough memory for array cAlpha[%d]! ",i);
            exit(1);
        }
    }
    for(i=0;i<i_max;i++)
        for(j=0;j<j_max;j++)
        {   cAlpha[i][j]=(char *)malloc(l_max*sizeof(char));
            if(cAlpha[i][j]==NULL)
            {   printf("\n No enough memory for array cAlpha[%d][%d]! ",i,j);
                exit(1);
            }
        }



/****************************************************************************/
/******************************* Calculation ********************************/
/****************************************************************************/

    if(iOrientationMode==0)
    {   cos_psi=cos(-psi);      sin_psi=sin(-psi);
        cos_iota=cos(-iota);    sin_iota=sin(-iota);
        cos_kappa=cos(-kappa);  sin_kappa=sin(-kappa);

        for(m=0;m<iNoScPl;m++)
        {   vNormal[m]=Rotation_z(cos_psi,sin_psi,vNormal_ini[m]);
            vNormal[m]=Rotation_y(cos_iota,sin_iota,vNormal[m]);
            vNormal[m]=Rotation_x(cos_kappa,sin_kappa,vNormal[m]);

            for(i=0;i<iNoA;i++)
            {   ksc[m][i]=Rotation_z(cos_psi,sin_psi,ksc_ini[m][i]);
                ksc[m][i]=Rotation_y(cos_iota,sin_iota,ksc[m][i]);
                ksc[m][i]=Rotation_x(cos_kappa,sin_kappa,ksc[m][i]);
            }
        }

        kinc=Rotation_z(cos_psi,sin_psi,kinc_ini);
        kinc=Rotation_y(cos_iota,sin_iota,kinc);
        kinc=Rotation_x(cos_kappa,sin_kappa,kinc);

        for(m=0;m<2;m++)
        {   Einc[m]=Rotation_z(cos_psi,sin_psi,Einc_ini[m]);
            Einc[m]=Rotation_y(cos_iota,sin_iota,Einc[m]);
            Einc[m]=Rotation_x(cos_kappa,sin_kappa,Einc[m]);
        }
    }

    // Volume calculations
    if(iChoiceStructure==1||iChoiceStructure==4||iChoiceStructure==6
       ||iChoiceStructure==7||iChoiceStructure==8||iChoiceStructure==14)
    {   iNoD=iGenerateEllipsoid(i_max,j_max,l_max,1.,1.,cAlpha);
        rVolume=(PI/6.)*(i_max*j_max*l_max);
    }

    if(iChoiceStructure==2||iChoiceStructure==5||iChoiceStructure==10||
        iChoiceStructure==13)
    {   iNoD=iGenerateEllipsoid(i_max,j_max,l_max,b_a,c_a,cAlpha);
        rVolume=(PI/6.)*(i_max*(i_max*b_a)*(i_max*c_a));
    }

    if(iChoiceStructure==3||iChoiceStructure==9)
    {   iNoD=iGenerateParallelepiped(i_max,j_max,l_max,b_a,c_a,cAlpha);
        rVolume=i_max*j_max*l_max;
    }

    if(iChoiceStructure==11)
    {   iNoD=i_max*j_max*l_max; rVolume=(real)iNoD;
    }

    if(iChoiceStructure==12)
    {   iNoD=iGenerateBisphere(i_max,j_max,l_max,R1_R2,cAlpha);

        if(R1_R2>1.)
        {   R1=(real)j_max/2.;
            R2=R1/R1_R2;
        }
        else
        {   R2=(real)j_max/2.;
            R1=R2*R1_R2;
        }
        rVolume=4*PI/3.*(R1*R1*R1+R2*R2*R2);
    }

    printf("\n Number of cells in the initial matrix      = %d",iNoD);
    printf("\n kd = %.6f",(float)kd);
    printf("\n x  = %.6f",(float)x);

    if(iChoiceStructure==6||iChoiceStructure==7||iChoiceStructure==8)
    {   iWormsNumber=(int)(rho*(real)(i_max*j_max*l_max)/(8.*(radius_worm*radius_worm*radius_worm)*(0.75*lr_worm+1)))+1;
        length_worm=radius_worm*lr_worm;

        worms=(cylinder *)malloc(sizeof(cylinder)*iWormsNumber);
        if(worms==NULL)
        {   printf("\n No enough memory for array worm! ");
            exit(1);
        }

        printf("\n Number of worms = %d",iWormsNumber);
        printf("\n Volume density of particles = %4.2lf",rho);
    }
    fflush(stdout);

    cListMaterials=(cplx *)malloc((iSorts+1)*sizeof(cplx));
    iListMV=(int *)malloc((iSorts+1)*sizeof(int));
    if(cListMaterials==NULL||iListMV==NULL)
    {   printf("\n No enough memory for arrays cListMaterials or iListMV! ");
        exit(1);
    }

    cListMaterials[0]=zero();   iListMV[0]=0;
    for(i=0;i<iSorts;i++)
    {   alpha=Polarizability_CM(kd,cRefrIndex[i]);

        if(cRefrIndex[i].r==1.&&cRefrIndex[i].i==0.)
            iListMV[(i+1)]=0;
        else
            iListMV[(i+1)]=1;

        if(iChoicePolarizability==2)
            alpha=Polarizability_CMRR(alpha);

        if(iChoicePolarizability==3)
            alpha=Polarizability_CMNR(kd,alpha);

        if(iChoicePolarizability==4||iChoicePolarizability==5)
            alpha=Polarizability_ILDR(kd,cRefrIndex[i],alpha);

        cListMaterials[(i+1)]=alpha;

        if(iChoiceStructure==4||iChoiceStructure==5||iChoiceStructure==9
           ||iChoiceStructure==10||iChoiceStructure==11||iChoiceStructure==12
           ||iChoiceStructure==13)
        {   if(iChoiceStructure==11||iChoiceStructure==12)
            {   printf("\n material %d: m = %.6f %.6fi",(i+1),cRefrIndex[i].r,cRefrIndex[i].i);
                printf("\n          alpha = %.6f %.6fi",alpha.r,alpha.i);
            }
            else
            {   printf("\n material %d(%d): m = %.6f %.6fi",(i+1),iSeeds[i],cRefrIndex[i].r,cRefrIndex[i].i);
                printf("\n              alpha = %.6f %.6fi",alpha.r,alpha.i);
            }
        }
        else
        {   printf("\n material: m = %.6f %.6fi",cRefrIndex[i].r,cRefrIndex[i].i);
            printf("\n       alpha = %.6f %.6fi",alpha.r,alpha.i);
        }
    }
    printf("\n");   fflush(stdout);

    if(iChoiceStructure==11)
    {   output_data=fopen(filename2,"r");   // Open initial data file for reading
        if(output_data==NULL)
        {   printf("Cannot read the initial data file!\n");
            exit(2);
        }
                                            // Reading initial data
        fseek(output_data,0,SEEK_SET);
        fscanf(output_data,"%s",filename2);
        iNoD=atoi(filename2);

        for(X=0;X<i_max;X++)
            for(Y=0;Y<j_max;Y++)
                for(Z=0;Z<l_max;Z++)
                    cAlpha[X][Y][Z]=0;

        for(i=0;i<iNoD;i++)
        {   fscanf(output_data,"%s",filename2); X=atoi(filename2);
            fscanf(output_data,"%s",filename2); Y=atoi(filename2);
            fscanf(output_data,"%s",filename2); Z=atoi(filename2);
            fscanf(output_data,"%s",filename2); T=atoi(filename2);

            if(X>=i_max||Y>=j_max||Z>=l_max)
            {   printf("\n Sample particle in the file is out of the expected boundaries.");
                printf("\n The program will be terminated. Please, check the file with particle.");
                fflush(stdout);
                exit(3);
            }

            cAlpha[X][Y][Z]=T;
        }
        fclose(output_data);
    }

/****************************************************************************/
/********************** Parallel solution over the samples ******************/
/****************************************************************************/

/* The samples are independent of one another: each one generates its own
   particle, solves the scattering problem for it and contributes to the
   averages. They are therefore distributed over several workers, each with
   its own occupancy map and its own set of solver arrays. The A-matrix and
   the FFT tables stay shared and are only read. */

    iFirstRealization=iRealizationCounter;
    iTasks=iRealizationNumber-iFirstRealization;
    iWorkers=worker_count(i_max,j_max,l_max,opt.threads);

/* No point in allocating more workspaces than there is work to do. */

    if(iWorkers>iTasks) iWorkers=iTasks;
    if(iWorkers<1)      iWorkers=1;

    printf("\n Samples to process: %d",iTasks);
    printf("\n Worker threads: %d",iWorkers);
    fflush(stdout);

    psiFixed=psi;   iotaFixed=iota;   kappaFixed=kappa;

    vWorkers=(dda_worker *)malloc(iWorkers*sizeof(dda_worker));
    if(vWorkers==NULL)
    {   printf("\n No enough memory for the worker descriptions! ");
        exit(1);
    }

    for(T=0;T<iWorkers;T++)
    {   dda_worker * w=&vWorkers[T];

        w->cAlpha=alloc_char3(i_max,j_max,l_max);
        for(i=0;i<i_max;i++)
            for(j=0;j<j_max;j++)
                for(l=0;l<l_max;l++)
                    w->cAlpha[i][j][l]=cAlpha[i][j][l];

        w->cListMaterials=(cplx *)malloc((iSorts+1)*sizeof(cplx));
        w->vNormal=(vector *)malloc(iNoScPl*sizeof(vector));
        w->ksc=(vector **)malloc(iNoScPl*sizeof(vector *));
        w->mmAcc=(mmuller *)malloc(iNoA*sizeof(mmuller));
        w->mm0=(mmuller *)malloc(iNoA*sizeof(mmuller));
        if(w->cListMaterials==NULL||w->vNormal==NULL||w->ksc==NULL||
           w->mmAcc==NULL||w->mm0==NULL)
        {   printf("\n No enough memory for the worker state! ");
            exit(1);
        }

        for(i=0;i<=iSorts;i++)  w->cListMaterials[i]=cListMaterials[i];
        for(i=0;i<iNoA;i++)     w->mmAcc[i]=mmNull();

        for(i=0;i<iNoScPl;i++)
        {   w->ksc[i]=(vector *)malloc(iNoA*sizeof(vector));
            if(w->ksc[i]==NULL)
            {   printf("\n No enough memory for the scattering directions! ");
                exit(1);
            }
            for(j=0;j<iNoA;j++) w->ksc[i][j]=ksc[i][j];
            w->vNormal[i]=vNormal[i];
        }

        w->kinc=kinc;

        for(i=0;i<2;i++)
        {   w->Einc[i]=Einc[i];
            w->Esc[i]=(cvector **)malloc(iNoScPl*sizeof(cvector *));
            if(w->Esc[i]==NULL)
            {   printf("\n No enough memory for the scattered field! ");
                exit(1);
            }
            for(j=0;j<iNoScPl;j++)
            {   w->Esc[i][j]=(cvector *)malloc(iNoA*sizeof(cvector));
                if(w->Esc[i][j]==NULL)
                {   printf("\n No enough memory for the scattered field! ");
                    exit(1);
                }
            }
        }

        w->worms=NULL;
        if(iChoiceStructure==6||iChoiceStructure==7||iChoiceStructure==8)
        {   w->worms=(cylinder *)malloc(sizeof(cylinder)*iWormsNumber);
            if(w->worms==NULL)
            {   printf("\n No enough memory for the worm geometry! ");
                exit(1);
            }
        }

        dda_workspace_create(&w->ws,i_max,j_max,l_max,iMode_SolutionMethod);
        w->ws.iVerbose=(iWorkers==1);
    }

    {   std::mutex outMutex;
        std::atomic<int> iNextTask(0);

        auto run_samples = [&](int iWorkerId)
        {   dda_worker * w=&vWorkers[iWorkerId];

/* The names below are the ones the sample loop used when it lived in main();
   the mutable ones are now private to this worker. */

            char *** cAlpha=w->cAlpha;
            cplx * cListMaterials=w->cListMaterials;
            vector * vNormal=w->vNormal;
            vector ** ksc=w->ksc;
            vector & kinc=w->kinc;
            cvector (&Einc)[2]=w->Einc;
            cvector ** (&Esc)[2]=w->Esc;
            cylinder * worms=w->worms;
            mmuller * mm0=w->mm0;

            int i,j,l,m;
            int iNoD=0,iSampleNumber;
            int X_min=0,X_max=0,Y_min=0,Y_max=0,Z_min=0,Z_max=0;
            real psiW=psiFixed,iotaW=iotaFixed,kappaW=kappaFixed;
            real cos_psiW=0.,sin_psiW=0.,cos_iotaW=0.,sin_iotaW=0.;
            real cos_kappaW=0.,sin_kappaW=0.;
            real VectTrans[3][3],rCabs,rCext,rCabsSum,rCextSum;
            cplx alphaW;
            mjones mjScat;

            dda_srand(random_base_seed()
                      +0x9E3779B97F4A7C15ULL*(unsigned long long)(iWorkerId+1));

            for(;;)
            {   int iTask=iNextTask.fetch_add(1);
                if(iTask>=iTasks)   break;

                iSampleNumber=iFirstRealization+iTask+1;
                rCabsSum=0.;    rCextSum=0.;

       if(iChoiceStructure==1)
                iNoD=iGenerateEllipsoid(i_max,j_max,l_max,1.,1.,cAlpha);

            if(iChoiceStructure==2)
                iNoD=iGenerateEllipsoid(i_max,j_max,l_max,b_a,c_a,cAlpha);

            if(iChoiceStructure==3)
                iNoD=iGenerateParallelepiped(i_max,j_max,l_max,b_a,c_a,cAlpha);

            if(iChoiceStructure==4)
                iNoD=iGenerateRandomEllipsoid(i_max,j_max,l_max,1.,1.,iSorts,iSeeds,
                                              cAlpha);
            if(iChoiceStructure==5)
                iNoD=iGenerateRandomEllipsoid(i_max,j_max,l_max,b_a,c_a,iSorts,
                                              iSeeds,cAlpha);

            if(iChoiceStructure==6) // the anthrax-like particle
                iNoD=iGenerateAnthrax(i_max,j_max,l_max,radius_worm,length_worm,
                                      iWormsNumber,worms,cAlpha);

            if(iChoiceStructure==7) // the hedgehog-like particle
                iNoD=iGenerateHedgehog(i_max,j_max,l_max,radius_worm,length_worm,
                                       iWormsNumber,worms,cAlpha);

            if(iChoiceStructure==8) // the dead anthrax-like particle
                iNoD=iGenerateDeadAnthrax(i_max,j_max,l_max,radius_worm,length_worm,
                                          iWormsNumber,worms,cAlpha);

            if(iChoiceStructure==9)
                iNoD=iGenerateRandomParallelepiped(i_max,j_max,l_max,b_a,c_a,
                                                   iSorts,iSeeds,cAlpha);

            if(iChoiceStructure==10)
                iNoD=iGenerateNibbleEllipsoid(i_max,j_max,l_max,b_a,c_a,iSorts,
                                              iSeeds,iBites,rDepth,iListMV,cAlpha);

            if(iChoiceStructure==12)
                iNoD=iGenerateBisphere(i_max,j_max,l_max,R1_R2,cAlpha);

            if(iChoiceStructure==13)
                iNoD=iGenerateHierarchicalEllipsoid(iHierarchicalLevels,i_max,
                                        j_max,l_max,b_a,c_a,iSorts,iSeeds,cAlpha);

            if(iChoiceStructure==14)
                iNoD=iGenerateSphereSpec(i_max,j_max,l_max,i_max/2.,rho,cAlpha);

            if(iWorkers==1)  printf("\n Sample #%d (%d dipoles)",iSampleNumber,iNoD);

            Particle_Boundary(iNoD,i_max,j_max,l_max,cAlpha,cListMaterials,
                              &X_min,&Y_min,&Z_min,&X_max,&Y_max,&Z_max);

            if(iWorkers==1)  printf("\n The sample occupies the initial matrix within boundaries:");
            if(iWorkers==1)  printf("\n axis X - from %3d to %3d",X_min,X_max);
            if(iWorkers==1)  printf("\n axis Y - from %3d to %3d",Y_min,Y_max);
            if(iWorkers==1)  printf("\n axis Z - from %3d to %3d",Z_min,Z_max);

            if(iOrientationMode==0)
            {   if(iWorkers==1)  printf("\n Orientation of the sample (via rotations around the coordinate axes):");
                if(iWorkers==1)  printf("\n axis Z = %.2f    axis Y = %.2f    axis X = %.2f",psiW*180./PI,iotaW*180./PI,kappaW*180./PI);
                fflush(stdout);
            }

            if(iOrientationMode==1)
            {   cos_psiW=(real)dda_rand()*2./(real)DDA_RAND_MAX-1.;
                psiW=acos(cos_psiW);      sin_psiW=sin(psiW);       // c1 and s1

                iotaW=(real)dda_rand()*2.*PI/(real)DDA_RAND_MAX;
                cos_iotaW=cos(iotaW);     sin_iotaW=sin(iotaW);     // c2 and s2

                kappaW=(real)dda_rand()*2.*PI/(real)DDA_RAND_MAX;
                cos_kappaW=cos(kappaW);   sin_kappaW=sin(kappaW);   // c3 and s3

                VectTrans[0][0]=cos_iotaW*cos_kappaW-cos_psiW*sin_iotaW*sin_kappaW;
                VectTrans[0][1]=sin_iotaW*cos_kappaW+cos_psiW*cos_iotaW*sin_kappaW;
                VectTrans[0][2]=sin_psiW*sin_kappaW;

                VectTrans[1][0]=-cos_iotaW*sin_kappaW-cos_psiW*sin_iotaW*cos_kappaW;
                VectTrans[1][1]=-sin_iotaW*sin_kappaW+cos_psiW*cos_iotaW*cos_kappaW;
                VectTrans[1][2]=sin_psiW*cos_kappaW;

                VectTrans[2][0]=sin_psiW*sin_iotaW;
                VectTrans[2][1]=-sin_psiW*cos_iotaW;
                VectTrans[2][2]=cos_psiW;

                for(m=0;m<iNoScPl;m++)
                {   for(i=0;i<iNoA;i++)
                    {   ksc[m][i].x=VectTrans[0][0]*ksc_ini[m][i].x+VectTrans[0][1]*ksc_ini[m][i].y+VectTrans[0][2]*ksc_ini[m][i].z;
                        ksc[m][i].y=VectTrans[1][0]*ksc_ini[m][i].x+VectTrans[1][1]*ksc_ini[m][i].y+VectTrans[1][2]*ksc_ini[m][i].z;
                        ksc[m][i].z=VectTrans[2][0]*ksc_ini[m][i].x+VectTrans[2][1]*ksc_ini[m][i].y+VectTrans[2][2]*ksc_ini[m][i].z;
                    }

                    vNormal[m].x=VectTrans[0][0]*vNormal_ini[m].x+VectTrans[0][1]*vNormal_ini[m].y+VectTrans[0][2]*vNormal_ini[m].z;
                    vNormal[m].y=VectTrans[1][0]*vNormal_ini[m].x+VectTrans[1][1]*vNormal_ini[m].y+VectTrans[1][2]*vNormal_ini[m].z;
                    vNormal[m].z=VectTrans[2][0]*vNormal_ini[m].x+VectTrans[2][1]*vNormal_ini[m].y+VectTrans[2][2]*vNormal_ini[m].z;
                }

                kinc.x=VectTrans[0][0]*kinc_ini.x+VectTrans[0][1]*kinc_ini.y+VectTrans[0][2]*kinc_ini.z;
                kinc.y=VectTrans[1][0]*kinc_ini.x+VectTrans[1][1]*kinc_ini.y+VectTrans[1][2]*kinc_ini.z;
                kinc.z=VectTrans[2][0]*kinc_ini.x+VectTrans[2][1]*kinc_ini.y+VectTrans[2][2]*kinc_ini.z;

                for(m=0;m<2;m++)
                {   Einc[m].x=VectTrans[0][0]*Einc_ini[m].x+VectTrans[0][1]*Einc_ini[m].y+VectTrans[0][2]*Einc_ini[m].z;
                    Einc[m].y=VectTrans[1][0]*Einc_ini[m].x+VectTrans[1][1]*Einc_ini[m].y+VectTrans[1][2]*Einc_ini[m].z;
                    Einc[m].z=VectTrans[2][0]*Einc_ini[m].x+VectTrans[2][1]*Einc_ini[m].y+VectTrans[2][2]*Einc_ini[m].z;
                }

                if(iWorkers==1)  printf("\n Orientation of the sample (via Eiler's angles):");
                if(iWorkers==1)  printf("\n first = %.2f    second = %.2f    third = %.2f",psiW*180./PI,iotaW*180./PI,kappaW*180./PI);
                fflush(stdout);
            }

                                            // Solve the scattering problem for two
            for(m=0;m<2;m++)                // polarization of incident field
            {   if(iChoicePolarizability==5)
                    for(i=0;i<iSorts;i++)
                    {   alphaW=Polarizability_CM(kd,cRefrIndex[i]);
                        alphaW=Polarizability_LDR(kd,cRefrIndex[i],alphaW,kinc,Einc[m]);
                                    cListMaterials[(i+1)]=alphaW;
                    }

                DDA_FFT(m,iMode_SolutionMethod,rTau,i_max,j_max,l_max,cListMaterials,
                        iListMV,cAlpha,Max_Step,Error_max,kd,kinc,Einc[m],
                        cW_2_i,cW_2_j,cW_2_l,cW_i,cW_j,cW_l,Im_cW_i,Im_cW_j,
                        Im_cW_l,cW2_i,cW2_j,cW2_l,iMode_FFT,cmA,iNoScPl,iNoA,ksc,
                        Esc[m],&rCabs,&rCext,&w->ws);

                rCabsSum+=.5*rCabs;    rCextSum+=.5*rCext;
            }

                                            // Transform to Muller matrix
            for(m=0;m<iNoA;m++)
            {   mm0[m]=mmNull();

                for(i=0;i<iNoScPl;i++)
                {   mjScat=BuildMatrixJones(kinc,ksc[i][m],vNormal[i],Einc[0],
                                            Einc[1],Esc[0][i][m],Esc[1][i][m]);

                    mm0[m]=mm0[m]+JonesToMuller(mjScat);
                }

                mm0[m]=mm0[m]/(real)iNoScPl;
            }




/* Accumulate the result of this sample and refresh the state file. */

                {   std::lock_guard<std::mutex> lock(outMutex);

                    Cabs+=rCabsSum;     Cext+=rCextSum;
                    for(m=0;m<iNoA;m++)
                        mmAccumulator[m]=mmAccumulator[m]+mm0[m];
                    iRealizationCounter++;

                    if(iWorkers>1)
                    {   printf("\n Sample #%d done: %d dipoles, %d plane(s), %d angle(s)",
                               iSampleNumber,iNoD,iNoScPl,iNoA);
                        fflush(stdout);
                    }

                    fseek(storage_data,0,SEEK_SET);
                    l=fwrite(&iRealizationCounter,sizeof(iRealizationCounter),1,
                             storage_data);
                    l=fwrite(&Cabs,sizeof(Cabs),1,storage_data);
                    l=fwrite(&Cext,sizeof(Cext),1,storage_data);
                    for(i=0;i<iNoA;i++)
                        l=fwrite(&mmAccumulator[i],sizeof(mmAccumulator[i]),1,
                                 storage_data);

                    fflush(storage_data);
                }
            }
        };

        {   std::vector<std::thread> vThreads;

            vThreads.reserve(iWorkers);
            for(T=0;T<iWorkers;T++)
                vThreads.emplace_back(run_samples,T);
            for(T=0;T<(int)vThreads.size();T++)
                vThreads[T].join();
        }
    }

    for(T=0;T<iWorkers;T++)
    {   dda_worker * w=&vWorkers[T];

        free_char3(w->cAlpha);
        free(w->cListMaterials);
        free(w->vNormal);
        free(w->mmAcc);
        free(w->mm0);
        for(i=0;i<iNoScPl;i++)
        {   free(w->ksc[i]);
            free(w->Esc[0][i]);
            free(w->Esc[1][i]);
        }
        free(w->ksc);   free(w->Esc[0]);    free(w->Esc[1]);
        free(w->worms);
        dda_workspace_free(&w->ws);
    }
    free(vWorkers);

    fclose(storage_data);

/****************************************************************************/
/******************************* Memory free ********************************/
/****************************************************************************/

    free(cListMaterials);   free(iListMV);
    free(vNormal);  free(vNormal_ini);

    for(m=0;m<iNoScPl;m++)
    {   free(ksc_ini[m]);
        free(ksc[m]);
        free(Esc[0][m]);
        free(Esc[1][m]);
    }
    free(ksc_ini);  free(ksc);  free(Esc[0]);   free(Esc[1]);


    for(i=0;i<i_max;i++)
        for(j=0;j<j_max;j++)
            free(cAlpha[i][j]);

    for(i=0;i<i_max;i++)
        free(cAlpha[i]);

    free(cAlpha);



    for(i=0;i<2*i_max;i++)
        for(j=0;j<2*j_max;j++)
            free(cmA[i][j]);

    for(i=0;i<2*i_max;i++)
        free(cmA[i]);

    free(cmA);


    if(iMode_FFT==0)
    {   free(cW_2_i);   free(cW_2_j);   free(cW_2_l);
        free(cW_i);     free(cW_j);     free(cW_l);
    }

    if(iMode_FFT==1)
    {   for(i=0;i<i_max;i++)
        {   free(cW2_i[i]); free(Im_cW_i[i]);
        }

        for(j=0;j<j_max;j++)
        {   free(cW2_j[j]); free(Im_cW_j[j]);
        }

        for(l=0;l<l_max;l++)
        {   free(cW2_l[l]); free(Im_cW_l[l]);
        }

        free(cW2_i);    free(cW2_j);    free(cW2_l);
        free(Im_cW_i);  free(Im_cW_j);  free(Im_cW_l);
    }

    if(iChoiceStructure==6||iChoiceStructure==7||iChoiceStructure==8)
        free(worms);

/****************************************************************************/
/************************** Write the result to file ************************/
/****************************************************************************/

    NO_CALC: printf("\n Averaging over %d samples is done.",iRealizationCounter);
    fflush(stdout);

    output_data=fopen(filename1,"w");
    if(output_data==NULL)
    {   printf ("Cannot read the output data file!\n");
        exit(2);
    }

    kd=1./(real)iRealizationCounter;
    fprintf(output_data,"Cabs=%.8f\t\tCext=%.8f",(float)(kd*Cabs),
                                                 (float)(kd*Cext));
    for(m=0;m<iNoA;m++)
    {   theta=theta_inf+theta_step*m;

        fprintf(output_data,"\n%7.3f\t",(float)theta);
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[0][0]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[0][1]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[0][2]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[0][3]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[1][0]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[1][1]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[1][2]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[1][3]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[2][0]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[2][1]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[2][2]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[2][3]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[3][0]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[3][1]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[3][2]));
        fprintf(output_data,"\t%+10.8e",(float)(kd*mmAccumulator[m].e[3][3]));
    }
    fclose(output_data);

    free(filename1);    free(filename2);    free(mmAccumulator);    free(mm0);

    time_end=time(NULL);
    printf("\n\nElapsed time = %ld sec.",(long)(time_end-time_begin));
    return 0;
}
