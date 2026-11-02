
//Time dependence is exp(-iwt)
//Refractive index is n=n'+in''

#include "function.hpp"
#include "randomize.hpp"


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

    int i,j,l,m,iNoA,iNoD,iFlag,iSorts,iBites,X,Y,Z,T;
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
    dda_workspace ws;
    FILE * output_data,* storage_data;

    printf("\nComputation of light scattering by arbitrary target:");
    printf("\nDiscrete Dipole Approximation (DDA) with Fast Fourier Transform (FFT)");
    printf("\nVersion 6.3, Copyright (C) 2001-2012 E.S.Zubko\n");
    fflush(stdout);

    filename1=(char *)malloc(1024*sizeof(char));
    filename2=(char *)malloc(1024*sizeof(char));
    random_init();
    //srand((unsigned)time(NULL));
    time_begin=time(NULL);

/****************************************************************************/
/******************************** Input data ********************************/
/****************************************************************************/

    if(argc>1)  theta_inf=atof(argv[1]);
    else
    {   printf("\n Input the bottom limit of range of phase angles ");
        scanf("%f",&theta_inf);
    }

    if(argc>2)  theta_up=atof(argv[2]);
    else
    {   printf("\n Input the upper limit of range of phase angles ");
        scanf("%f",&theta_up);
    }

    if(argc>3)  theta_step=atof(argv[3]);
    else
    {   printf("\n Input the increment of phase angles ");
        scanf("%f",&theta_step);
    }

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

    if(argc>4)  strcpy(filename1,argv[4]);
    else
    {   printf("\n Input name for the binary result file (dst) ");
        scanf("%s",filename1);
    }

    // begin of reading of the dst-file
    storage_data=fopen(filename1,"rb+");
    if(storage_data==NULL)          // file doesn't exist; it needs to be created
    {   storage_data=fopen(filename1,"wb+");
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

    storage_data=fopen(filename1,"rb+");
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

    if(argc>5)  strcpy(filename1,argv[5]);
    else
    {   printf("\n Input name for the text result file (dat) ");
        scanf("%s",filename1);
    }

    if(argc>6)  iRealizationNumber=atoi(argv[6]);
    else
    {   printf("\n Input number of samples/orientations ");
        scanf("%i",&iRealizationNumber);
    }

    iOrientationMode=-1;
    while(iOrientationMode!=0&&iOrientationMode!=1)
    if(argc>7)  iOrientationMode=atoi(argv[7]);
    else
    {   printf("\n Please, choose the orientation mode ");
        printf("\n 0 - no averaging ");
        printf("\n 1 - Monte-Carlo averaging (random orientations) ");
        scanf("%i",&iOrientationMode);
    }

    if(iOrientationMode==0)
    {   if(argc>8)  psi=atof(argv[8]);
        else
        {   printf("\n Input the angle psi for the first rotation around axis Z (in degree) ");
            scanf("%f",&psi);
        }
        psi=psi*PI/180.;

        if(argc>9)  iota=atof(argv[9]);
        else
        {   printf("\n Input the angle iota for the second rotation around axis Y (in degree) ");
            scanf("%f",&iota);
        }
        iota=iota*PI/180.;

        if(argc>10)  kappa=atof(argv[10]);
        else
        {   printf("\n Input the angle kappa for the third rotation around axis X (in degree) ");
            scanf("%f",&kappa);
        }
        kappa=kappa*PI/180.;

        if(argc>11) iNoScPl=atoi(argv[11]);
        else
        {   printf("\n Input number of the scattering planes ");
            scanf("%d",&iNoScPl);
        }

        iFlag=4;
    }

    if(iOrientationMode==1)
    {   if(argc>8)  iNoScPl=atoi(argv[8]);
        else
        {   printf("\n Input number of the scattering planes ");
            scanf("%d",&iNoScPl);
        }

        iFlag=1;
    }

    iChoiceStructure=0;
    while(iChoiceStructure!=1&&iChoiceStructure!=2&&iChoiceStructure!=3&&
          iChoiceStructure!=4&&iChoiceStructure!=5&&iChoiceStructure!=6&&
          iChoiceStructure!=7&&iChoiceStructure!=8&&iChoiceStructure!=9&&
          iChoiceStructure!=10&&iChoiceStructure!=11&&iChoiceStructure!=12&&
          iChoiceStructure!=13&&iChoiceStructure!=14)
    if(argc>(8+iFlag))  iChoiceStructure=atoi(argv[(8+iFlag)]);
    else
    {   printf("\n Please, choose the morphology type of target ");
        printf("\n Input  1 for sphere ");
        printf("\n Input  2 for ellipsoid ");
        printf("\n Input  3 for parallelepiped ");
        printf("\n Input  4 for random particle in spherical matrix ");
        printf("\n Input  5 for random particle in ellipsoidal matrix ");
        printf("\n Input  6 for random antrax-like particle ");
        printf("\n Input  7 for random hedgehog-like particle ");
        printf("\n Input  8 for random dead antrax-like particle ");
        printf("\n Input  9 for random particle in parallelepiped matrix ");
        printf("\n Input 10 for random nibble particle in ellipsoidal matrix ");
        printf("\n Input 11 for particle from file ");
        printf("\n Input 12 for two spheres in contact (along X-axis) ");
        printf("\n Input 13 for hierarchical random particle in ellipsoidal matrix ");
        printf("\n Input 14 for sphere with arbitrary radius\n");
        scanf("%i",&iChoiceStructure);
    }

    if(iChoiceStructure!=4&&iChoiceStructure!=5&&iChoiceStructure!=9
       &&iChoiceStructure!=10&&iChoiceStructure!=11&&iChoiceStructure!=12)
    {   iSeeds=(int *)malloc(sizeof(int));
        cRefrIndex=(cplx *)malloc(sizeof(cplx));
        if(cRefrIndex==NULL||iSeeds==NULL)
        {   printf("\n No enough memory for array cRefrIndex! ");
            exit(1);
        }
        iSorts=1;
    }

    if(iChoiceStructure==1)
    {   if(argc>(9+iFlag)) cRefrIndex[0].r=atof(argv[(9+iFlag)]);
        else
        {   printf("\n Input the real part of refractive index of the constituent material ");
            scanf("%f",&cRefrIndex[0].r);
        }

        if(argc>(10+iFlag)) cRefrIndex[0].i=atof(argv[(10+iFlag)]);
        else
        {   printf("\n Input the imaginary part of refractive index of the constituent material ");
            scanf("%f",&cRefrIndex[0].i);
        }

        iFlag+=2;
    }

    if(iChoiceStructure==2)
    {   if(argc>(9+iFlag))  b_a=atof(argv[(9+iFlag)]);
        else
        {   printf("\n Input the aspect ratio of the ellipsoid b/a ");
            scanf("%f",&b_a);
        }

        if(argc>(10+iFlag)) c_a=atof(argv[(10+iFlag)]);
        else
        {   printf("\n Input the aspect ratio of the ellipsoid c/a ");
            scanf("%f",&c_a);
        }

        if(argc>(11+iFlag)) cRefrIndex[0].r=atof(argv[(11+iFlag)]);
        else
        {   printf("\n Input the real part of refractive index of the constituent material ");
            scanf("%f",&cRefrIndex[0].r);
        }

        if(argc>(12+iFlag)) cRefrIndex[0].i=atof(argv[(12+iFlag)]);
        else
        {   printf("\n Input the imaginary part of refractive index of the constituent material ");
            scanf("%f",&cRefrIndex[0].i);
        }

        iFlag+=4;
    }

    if(iChoiceStructure==3)
    {   if(argc>(9+iFlag))  b_a=atof(argv[(9+iFlag)]);
        else
        {   printf("\n Input the aspect ratio of the parallelepiped b/a ");
            scanf("%f",&b_a);
        }

        if(argc>(10+iFlag)) c_a=atof(argv[(10+iFlag)]);
        else
        {   printf("\n Input the aspect ratio of the parallelepiped c/a ");
            scanf("%f",&c_a);
        }

        if(argc>(11+iFlag)) cRefrIndex[0].r=atof(argv[(11+iFlag)]);
        else
        {   printf("\n Input the real part of refractive index of the constituent material ");
            scanf("%f",&cRefrIndex[0].r);
        }

        if(argc>(12+iFlag)) cRefrIndex[0].i=atof(argv[(12+iFlag)]);
        else
        {   printf("\n Input the imaginary part of refractive index of the constituent material ");
            scanf("%f",&cRefrIndex[0].i);
        }

        iFlag+=4;
    }

    if(iChoiceStructure==4)
    {   if(argc>(9+iFlag))  iSorts=atoi(argv[(9+iFlag)]);
        else
        {   printf("\n Input number of the constituent materials ");
            scanf("%i",&iSorts);
        }

        iSeeds=(int *)malloc(iSorts*sizeof(int));
        cRefrIndex=(cplx *)malloc(iSorts*sizeof(cplx));
        if(iSeeds==NULL||cRefrIndex==NULL)
        {   printf("\n No enough memory for array iSeeds or cRefrIndex! ");
            exit(1);
        }

        for(i=0;i<iSorts;i++)
        {   if(argc>(10+iFlag+3*i)) iSeeds[i]=atoi(argv[(10+iFlag+3*i)]);
            else
            {   printf("\n Input number of %d material seeds ",(i+1));
                scanf("%i",&iSeeds[i]);
            }

            if(argc>(10+iFlag+3*i+1))   cRefrIndex[i].r=atof(argv[(10+iFlag+3*i+1)]);
            else
            {   printf("\n Input the real part of refractive index of %d material ",(i+1));
                scanf("%f",&cRefrIndex[i].r);
            }

            if(argc>(10+iFlag+3*i+2))   cRefrIndex[i].i=atof(argv[(10+iFlag+3*i+2)]);
            else
            {   printf("\n Input the imaginary part of refractive index of %d material ",(i+1));
                scanf("%f",&cRefrIndex[i].i);
            }
        }

        iFlag+=(1+3*iSorts);
    }

    if(iChoiceStructure==5)
    {   if(argc>(9+iFlag))  b_a=atof(argv[(9+iFlag)]);
        else
        {   printf("\n Input the aspect ratio of the ellipsoid b/a ");
            scanf("%f",&b_a);
        }

        if(argc>(10+iFlag))  c_a=atof(argv[(10+iFlag)]);
        else
        {   printf("\n Input the aspect ratio of the ellipsoid c/a ");
            scanf("%f",&c_a);
        }

        if(argc>(11+iFlag))  iSorts=atoi(argv[(11+iFlag)]);
        else
        {   printf("\n Input number of the constituent materials ");
            scanf("%i",&iSorts);
        }

        iSeeds=(int *)malloc(iSorts*sizeof(int));
        cRefrIndex=(cplx *)malloc(iSorts*sizeof(cplx));
        if(iSeeds==NULL||cRefrIndex==NULL)
        {   printf("\n No enough memory for array iSeeds or cRefrIndex! ");
            exit(1);
        }

        for(i=0;i<iSorts;i++)
        {   if(argc>(12+iFlag+3*i)) iSeeds[i]=atoi(argv[(12+iFlag+3*i)]);
            else
            {   printf("\n Input number of %d material seeds ",(i+1));
                scanf("%i",&iSeeds[i]);
            }

            if(argc>(12+iFlag+3*i+1))   cRefrIndex[i].r=atof(argv[(12+iFlag+3*i+1)]);
            else
            {   printf("\n Input the real part of refractive index of %d material ",(i+1));
                scanf("%f",&cRefrIndex[i].r);
            }

            if(argc>(12+iFlag+3*i+2))   cRefrIndex[i].i=atof(argv[(12+iFlag+3*i+2)]);
            else
            {   printf("\n Input the imaginary part of refractive index of %d material ",(i+1));
                scanf("%f",&cRefrIndex[i].i);
            }
        }

        iFlag+=3*(1+iSorts);
    }

    if(iChoiceStructure==6||iChoiceStructure==7||iChoiceStructure==8)
    {   if(argc>(9+iFlag))  rho=atof(argv[(9+iFlag)]);
        else
        {   printf("\n Input the volume density of the target particle ");
            scanf("%f",&rho);
        }

        if(argc>(10+iFlag)) radius_worm=atof(argv[(10+iFlag)]);
        else
        {   printf("\n Input radius of the monomer (in cells)");
            scanf("%f",&radius_worm);
        }

        if(argc>(11+iFlag)) lr_worm=atof(argv[(11+iFlag)]);
        else
        {   printf("\n Input the ratio length over radius for the monomer ");
            scanf("%f",&lr_worm);
        }

        if(argc>(12+iFlag)) cRefrIndex[0].r=atof(argv[(12+iFlag)]);
        else
        {   printf("\n Input the real part of refractive index of the constituent material ");
            scanf("%f",&cRefrIndex[0].r);
        }

        if(argc>(13+iFlag)) cRefrIndex[0].i=atof(argv[(13+iFlag)]);
        else
        {   printf("\n Input the imaginary part of refractive index of the constituent material ");
            scanf("%f",&cRefrIndex[0].i);
        }

        iFlag+=5;
    }

    if(iChoiceStructure==9)
    {   if(argc>(9+iFlag))  b_a=atof(argv[(9+iFlag)]);
        else
        {   printf("\n Input the aspect ratio of the ellipsoid b/a ");
            scanf("%f",&b_a);
        }

        if(argc>(10+iFlag))  c_a=atof(argv[(10+iFlag)]);
        else
        {   printf("\n Input the aspect ratio of the ellipsoid c/a ");
            scanf("%f",&c_a);
        }

        if(argc>(11+iFlag))  iSorts=atoi(argv[(11+iFlag)]);
        else
        {   printf("\n Input number of the constituent materials ");
            scanf("%i",&iSorts);
        }

        iSeeds=(int *)malloc(iSorts*sizeof(int));
        cRefrIndex=(cplx *)malloc(iSorts*sizeof(cplx));
        if(iSeeds==NULL||cRefrIndex==NULL)
        {   printf("\n No enough memory for array iSeeds or cRefrIndex! ");
            exit(1);
        }

        for(i=0;i<iSorts;i++)
        {   if(argc>(12+iFlag+3*i)) iSeeds[i]=atoi(argv[(12+iFlag+3*i)]);
            else
            {   printf("\n Input number of %d material seeds ",(i+1));
                scanf("%i",&iSeeds[i]);
            }

            if(argc>(12+iFlag+3*i+1))   cRefrIndex[i].r=atof(argv[(12+iFlag+3*i+1)]);
            else
            {   printf("\n Input the real part of refractive index of %d material ",(i+1));
                scanf("%f",&cRefrIndex[i].r);
            }

            if(argc>(12+iFlag+3*i+2))   cRefrIndex[i].i=atof(argv[(12+iFlag+3*i+2)]);
            else
            {   printf("\n Input the imaginary part of refractive index of %d material ",(i+1));
                scanf("%f",&cRefrIndex[i].i);
            }
        }

        iFlag+=3*(1+iSorts);
    }

    if(iChoiceStructure==10)
    {   if(argc>(9+iFlag))  b_a=atof(argv[(9+iFlag)]);
        else
        {   printf("\n Input the aspect ratio of the ellipsoid b/a ");
            scanf("%f",&b_a);
        }

        if(argc>(10+iFlag)) c_a=atof(argv[(10+iFlag)]);
        else
        {   printf("\n Input the aspect ratio of the ellipsoid c/a ");
            scanf("%f",&c_a);
        }

        if(argc>(11+iFlag)) iSorts=atoi(argv[(11+iFlag)]);
        else
        {   printf("\n Input number of the constituent materials ");
            scanf("%i",&iSorts);
        }

        if(argc>(12+iFlag)) iBites=atoi(argv[(12+iFlag)]);
        else
        {   printf("\n Input number of the outer ceed cells for empty space ");
            scanf("%i",&iBites);
        }

        if(argc>(13+iFlag)) rDepth=atof(argv[(13+iFlag)]);
        else
        {   printf("\n Input depth of the surface layer with these ceed cells (in radii of the initial matrix)" );
            scanf("%f",&rDepth);
        }

        iSeeds=(int *)malloc(iSorts*sizeof(int));
        cRefrIndex=(cplx *)malloc(iSorts*sizeof(cplx));
        if(iSeeds==NULL||cRefrIndex==NULL)
        {   printf("\n No enough memory for array iSeeds or cRefrIndex! ");
            exit(1);
        }

        for(i=0;i<iSorts;i++)
        {   if(argc>(14+iFlag+3*i)) iSeeds[i]=atoi(argv[(14+iFlag+3*i)]);
            else
            {   printf("\n Input number of %d material seeds ",(i+1));
                scanf("%i",&iSeeds[i]);
            }

            if(argc>(14+iFlag+3*i+1))   cRefrIndex[i].r=atof(argv[(14+iFlag+3*i+1)]);
            else
            {   printf("\n Input the real part of refractive index of %d material ",(i+1));
                scanf("%f",&cRefrIndex[i].r);
            }

            if(argc>(14+iFlag+3*i+2))   cRefrIndex[i].i=atof(argv[(14+iFlag+3*i+2)]);
            else
            {   printf("\n Input the imaginary part of refractive index of %d material ",(i+1));
                scanf("%f",&cRefrIndex[i].i);
            }
        }

        iFlag+=(5+3*iSorts);
    }

    if(iChoiceStructure==11)
    {   if(argc>(9+iFlag))  strcpy(filename2,argv[(9+iFlag)]);
        else
        {   printf("\n Input name of the file containing data on the target ");
            scanf("%s",filename2);
        }

        if(argc>(10+iFlag)) iSorts=atoi(argv[(10+iFlag)]);
        else
        {   printf("\n Input number of the constituent materials ");
            scanf("%i",&iSorts);
        }

        cRefrIndex=(cplx *)malloc(iSorts*sizeof(cplx));
        if(cRefrIndex==NULL)
        {   printf("\n No enough memory for array cRefrIndex! ");
            exit(1);
        }

        for(i=0;i<iSorts;i++)
        {   if(argc>(11+iFlag+2*i))     cRefrIndex[i].r=atof(argv[(11+iFlag+2*i)]);
            else
            {   printf("\n Input the real part of refractive index of %d material ",(i+1));
                scanf("%f",&cRefrIndex[i].r);
            }

            if(argc>(11+iFlag+2*i+1))   cRefrIndex[i].i=atof(argv[(11+iFlag+2*i+1)]);
            else
            {   printf("\n Input the imaginary part of refractive index of %d material ",(i+1));
                scanf("%f",&cRefrIndex[i].i);
            }
        }

        iFlag+=(2+2*iSorts);
    }

    if(iChoiceStructure==12)
    {   if(argc>(9+iFlag))  R1_R2=atof(argv[(9+iFlag)]);
        else
        {   printf("\n Input ratio of the first sphere radius over the second sphere radius ");
            scanf("%f",&R1_R2);
        }

        iSorts=2;
        iSeeds=(int *)malloc(iSorts*sizeof(int));
        cRefrIndex=(cplx *)malloc(iSorts*sizeof(cplx));
        if(iSeeds==NULL||cRefrIndex==NULL)
        {   printf("\n No enough memory for array iSeeds or cRefrIndex! ");
            exit(1);
        }

        for(i=0;i<iSorts;i++)
        {   if(argc>(10+iFlag+2*i))     cRefrIndex[i].r=atof(argv[(10+iFlag+2*i)]);
            else
            {   printf("\n Input the real part of refractive index of %d material ",(i+1));
                scanf("%f",&cRefrIndex[i].r);
            }

            if(argc>(10+iFlag+2*i+1))   cRefrIndex[i].i=atof(argv[(10+iFlag+2*i+1)]);
            else
            {   printf("\n Input the imaginary part of refractive index of %d material ",(i+1));
                scanf("%f",&cRefrIndex[i].i);
            }
        }

        iFlag+=5;
    }

    if(iChoiceStructure==13)
    {   if(argc>(9+iFlag))  b_a=atof(argv[(9+iFlag)]);
        else
        {   printf("\n Input the aspect ratio of the ellipsoid b/a ");
            scanf("%f",&b_a);
        }

        if(argc>(10+iFlag)) c_a=atof(argv[(10+iFlag)]);
        else
        {   printf("\n Input the aspect ratio of the ellipsoid c/a ");
            scanf("%f",&c_a);
        }

        if(argc>(11+iFlag)) iHierarchicalLevels=atoi(argv[(11+iFlag)]);
        else
        {   printf("\n Input number of hierarchical levels ");
            scanf("%i",&iHierarchicalLevels);
        }

        if(argc>(12+iFlag)) iSorts=atoi(argv[(12+iFlag)]);
        else
        {   printf("\n Input number of the constituent materials");
            printf("\n Note: in the current version of the program it should be 2 always;");
            printf("\n whereas, the first material should not be the 'free space'! ");
            scanf("%i",&iSorts);
        }

        iSeeds=(int *)malloc(iSorts*sizeof(int));
        cRefrIndex=(cplx *)malloc(iSorts*sizeof(cplx));
        if(iSeeds==NULL||cRefrIndex==NULL)
        {   printf("\n No enough memory for array iSeeds or cRefrIndex! ");
            exit(1);
        }

        for(i=0;i<iSorts;i++)
        {   if(argc>(13+iFlag+3*i))     iSeeds[i]=atoi(argv[(13+iFlag+3*i)]);
            else
            {   printf("\n Input number of %d material seeds ",(i+1));
                scanf("%i",&iSeeds[i]);
            }

            if(argc>(13+iFlag+3*i+1))   cRefrIndex[i].r=atof(argv[(13+iFlag+3*i+1)]);
            else
            {   printf("\n Input the real part of refractive index of %d material ",(i+1));
                scanf("%f",&cRefrIndex[i].r);
            }

            if(argc>(13+iFlag+3*i+2))   cRefrIndex[i].i=atof(argv[(13+iFlag+3*i+2)]);
            else
            {   printf("\n Input the imaginary part of refractive index of %d material ",(i+1));
                scanf("%f",&cRefrIndex[i].i);
            }
        }

        iFlag+=(4+3*iSorts);
    }

    if(iChoiceStructure==14)
    {   if(argc>(9+iFlag)) cRefrIndex[0].r=atof(argv[(9+iFlag)]);
        else
        {   printf("\n Input the real part of refractive index of the constituent material ");
            scanf("%f",&cRefrIndex[0].r);
        }

        if(argc>(10+iFlag)) cRefrIndex[0].i=atof(argv[(10+iFlag)]);
        else
        {   printf("\n Input the imaginary part of refractive index of the constituent material ");
            scanf("%f",&cRefrIndex[0].i);
        }

        if(argc>(11+iFlag)) rho=atof(argv[(11+iFlag)]);
        else
        {   printf("\n Input the density rho (the part of the initial spherical volume filled by dipoles) ");
            scanf("%f",&rho);
        }

        iFlag+=3;
    }

    if(argc>(9+iFlag))  i_max=atoi(argv[(9+iFlag)]);
    else
    {   printf("\n Input the maximum number of cells along the axis X (shoulf be 2^integer) ");
        scanf("%i",&i_max);
    }

    if(argc>(10+iFlag)) j_max=atoi(argv[(10+iFlag)]);
    else
    {   printf("\n Input the maximum number of cells along the axis Y (shoulf be 2^integer) ");
        scanf("%i",&j_max);
    }

    if(argc>(11+iFlag)) l_max=atoi(argv[(11+iFlag)]);
    else
    {   printf("\n Input the maximum number of cells along the axis Z (shoulf be 2^integer) ");
        scanf("%i",&l_max);
    }

    if(argc>(12+iFlag))  x=atof(argv[(12+iFlag)]);
    else
    {   printf("\n Input the size parameter x of the target (in respect to the side X) ");
        scanf("%f",&x);
    }

    iChoicePolarizability=0;
    while(iChoicePolarizability!=1&&iChoicePolarizability!=2&&iChoicePolarizability!=3
        &&iChoicePolarizability!=4&&iChoicePolarizability!=5)
    {   if(argc>(13+iFlag))  iChoicePolarizability=atoi(argv[(13+iFlag)]);
        else
        {   printf("\n Please, choose the relation between polarizability and refractive index.");
            printf("\n Input 1 for Clausius-Mossotti");
            printf("\n Input 2 for Clausius-Mossotti plus radiative reaction (by Draine)");
            printf("\n Input 3 for Clausius-Mossotti plus nonradiative reaction (by Goedecke&O'Brien)");
            printf("\n Input 4 for Isotropized Lattice Dispersion Relation (by Draine)");
            printf("\n Input 5 for Lattice Dispersion Relation (by Draine)");
            printf("\n Note: Do not choose 4 if the wave-vector of incident field is parallel");
            printf("\n to one of the coordinate axes, as it reduces the accuracy of computations!\n");
            scanf("%i",&iChoicePolarizability);
        }
    }

    if(argc>(14+iFlag)) Error_max=atof(argv[(14+iFlag)]);
    else
    {   printf("\n Input the parameter Error_max (it specifies the termination of the iterative process) ");
        scanf("%f",&Error_max);
    }

    if(argc>(15+iFlag)) Max_Step=atoi(argv[(15+iFlag)]);
    else
    {   printf("\n Input the maximum number of steps in the iterative process ");
        scanf("%i",&Max_Step);
    }

    iMode_FFT=-1;
    while(iMode_FFT!=0&&iMode_FFT!=1)
    if(argc>(16+iFlag)) iMode_FFT=atoi(argv[(16+iFlag)]);
    else
    {   printf("\n Please, choose the algorithm of the FFT ");
        printf("\n Input 0 for the Cooley-Tukey method ");
        printf("\n Input 1 for the Rader-Brenner method ");
        scanf("%i",&iMode_FFT);
    }

    iMode_SolutionMethod=-1;
    while(iMode_SolutionMethod!=0&&iMode_SolutionMethod!=1&&iMode_SolutionMethod!=2
          &&iMode_SolutionMethod!=3)
    if(argc>(17+iFlag)) iMode_SolutionMethod=atoi(argv[(17+iFlag)]);
    else
    {   printf("\n Please, choose the solution method of system of linear algebraic equations");
        printf("\n Input 0 for Conjugate Gradient Method (by Voevodin)");
        printf("\n Input 1 for AA* Minimal Iterations Method (by Voevodin)");
        printf("\n Input 2 for BiConjugate Gradient Stabilized Method");
        printf("\n Input 3 for Simple Iteration Method\n");
        scanf("%i",&iMode_SolutionMethod);
    }

    if(iMode_SolutionMethod==3)
    {   if(argc>(18+iFlag)) rTau=atof(argv[(18+iFlag)]);
        else
        {   printf("\n Input the weight factor for Simple Iteration Method ");
            scanf("%f",&rTau);
        }

        iFlag+=1;
    }


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

/* The scratch arrays of the solver are allocated once for the whole run. */

    dda_workspace_create(&ws,i_max,j_max,l_max,iMode_SolutionMethod);


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

    while(iRealizationCounter<iRealizationNumber)
    {   if(iChoiceStructure==1)
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

        iRealizationCounter++;

        printf("\n Sample #%d (%d dipoles)",iRealizationCounter,iNoD);

        Particle_Boundary(iNoD,i_max,j_max,l_max,cAlpha,cListMaterials,
                          &X_min,&Y_min,&Z_min,&X_max,&Y_max,&Z_max);

        printf("\n The sample occupies the initial matrix within boundaries:");
        printf("\n axis X - from %3d to %3d",X_min,X_max);
        printf("\n axis Y - from %3d to %3d",Y_min,Y_max);
        printf("\n axis Z - from %3d to %3d",Z_min,Z_max);

        if(iOrientationMode==0)
        {   printf("\n Orientation of the sample (via rotations around the coordinate axes):");
            printf("\n axis Z = %.2f    axis Y = %.2f    axis X = %.2f",psi*180./PI,iota*180./PI,kappa*180./PI);
            fflush(stdout);
        }

        if(iOrientationMode==1)
        {   cos_psi=(real)dda_rand()*2./(real)DDA_RAND_MAX-1.;
            psi=acos(cos_psi);      sin_psi=sin(psi);       // c1 and s1

            iota=(real)dda_rand()*2.*PI/(real)DDA_RAND_MAX;
            cos_iota=cos(iota);     sin_iota=sin(iota);     // c2 and s2

            kappa=(real)dda_rand()*2.*PI/(real)DDA_RAND_MAX;
            cos_kappa=cos(kappa);   sin_kappa=sin(kappa);   // c3 and s3

            VectTrans[0][0]=cos_iota*cos_kappa-cos_psi*sin_iota*sin_kappa;
            VectTrans[0][1]=sin_iota*cos_kappa+cos_psi*cos_iota*sin_kappa;
            VectTrans[0][2]=sin_psi*sin_kappa;

            VectTrans[1][0]=-cos_iota*sin_kappa-cos_psi*sin_iota*cos_kappa;
            VectTrans[1][1]=-sin_iota*sin_kappa+cos_psi*cos_iota*cos_kappa;
            VectTrans[1][2]=sin_psi*cos_kappa;

            VectTrans[2][0]=sin_psi*sin_iota;
            VectTrans[2][1]=-sin_psi*cos_iota;
            VectTrans[2][2]=cos_psi;

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

            printf("\n Orientation of the sample (via Eiler's angles):");
            printf("\n first = %.2f    second = %.2f    third = %.2f",psi*180./PI,iota*180./PI,kappa*180./PI);
            fflush(stdout);
        }

                                        // Solve the scattering problem for two
        for(m=0;m<2;m++)                // polarization of incident field
        {   if(iChoicePolarizability==5)
                for(i=0;i<iSorts;i++)
                {   alpha=Polarizability_CM(kd,cRefrIndex[i]);
                    alpha=Polarizability_LDR(kd,cRefrIndex[i],alpha,kinc,Einc[m]);
                    cListMaterials[(i+1)]=alpha;
                }

            DDA_FFT(m,iMode_SolutionMethod,rTau,i_max,j_max,l_max,cListMaterials,
                    iListMV,cAlpha,Max_Step,Error_max,kd,kinc,Einc[m],
                    cW_2_i,cW_2_j,cW_2_l,cW_i,cW_j,cW_l,Im_cW_i,Im_cW_j,
                    Im_cW_l,cW2_i,cW2_j,cW2_l,iMode_FFT,cmA,iNoScPl,iNoA,ksc,
                    Esc[m],&rCabs,&rCext,&ws);

            Cabs+=.5*rCabs;    Cext+=.5*rCext;
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

        for(m=0;m<iNoA;m++)
            mmAccumulator[m]=mmAccumulator[m]+mm0[m];

        fseek(storage_data,0,SEEK_SET);
        l=fwrite(&iRealizationCounter,sizeof(iRealizationCounter),1,storage_data);
        l=fwrite(&Cabs,sizeof(Cabs),1,storage_data);
        l=fwrite(&Cext,sizeof(Cext),1,storage_data);
        for(i=0;i<iNoA;i++)
            l=fwrite(&mmAccumulator[i],sizeof(mmAccumulator[i]),1,storage_data);

        fflush(storage_data);
    }
    fclose(storage_data);

/****************************************************************************/
/******************************* Memory free ********************************/
/****************************************************************************/

    free(iSeeds);   free(cRefrIndex);   free(cListMaterials);   free(iListMV);
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

    dda_workspace_free(&ws);
    

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
