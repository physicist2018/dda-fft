
#include "function.hpp"

/****************************************************************************/
/************** Binary inversion of given number (used in FFT) **************/
/****************************************************************************/

void Change_binary_number( int number_old, int number_of_bits, int * number_new )
{   int i,i0;
	short int * bits;
	double d0,d1;

	bits=(short int *)malloc(number_of_bits*sizeof(short int));
	if(bits==NULL)  printf("No enough memory for array bits");

/********************* Representation in binary sort ************************/

	d0=number_old;
	for(i=0;i<number_of_bits;i++)
	{
		d1=d0/2;    d0=(int)d1;

		if(d1-d0!=0)
			bits[i]=1;
		else
			bits[i]=0;
	}

/******************** Invesion of number in binary sort *********************/

	for(i=0;i<(int)(number_of_bits/2);i++)
	{
		i0=bits[i];
		bits[i]=bits[(number_of_bits-1)-i];
		bits[(number_of_bits-1)-i]=i0;
	}

/************ Definition of the decimal form of record of outcome ***********/

	* number_new=0;
	for(i=0;i<number_of_bits;i++)
		(* number_new)+=(int)pow(2.,(double)i)*bits[i];

	free(bits);
	return;
}

/****************************************************************************/
/************************* Definition number of bits ************************/
/****************************************************************************/

void Get_number_bits(int number, int * number_of_bits)
{   double d0;

	if(number==0)   d0=0;
	else            d0=log((double)number)/log(2.);

	(* number_of_bits)=(int)d0+1;

	return;
}

/****************************************************************************/
/******************* Calculation the exp(-i((2*PI*j)/(N))) ******************/
/****************************************************************************/

void FourierPhase(int i_max, cplx * cW)
{   int i;
	cplx cw=exp_i((real)(-2.*PI/(2*i_max)));

	cW[0].r=1.; cW[0].i=0.;
	cW[1]=cw;

	for(i=2;i<(i_max);i++)
		cW[i]=cW[i-1]*cw;

	return;
}

/****************************************************************************/
/********************************* FFT base *********************************/
/****************************************************************************/

void FFT(int lim_NoB, int i_max, cplx * cW, cmatrix * cmDATA_ini,
		 cmatrix * cmDATA)
{   int i,j,l,k,m,i_new,j_new,iNolFFT,iStep;
	cmatrix cm0;

/************************* The sorting of initial data **********************/

	for(i=0;i<i_max;i++)
	{   Change_binary_number(i,lim_NoB,&i_new);
		//cmDATA[i]=cmDATA_ini[i_new];
		Ccopy(cmDATA_ini[i_new].e00,cmDATA[i].e00);
		Ccopy(cmDATA_ini[i_new].e01,cmDATA[i].e01);
		Ccopy(cmDATA_ini[i_new].e02,cmDATA[i].e02);
		Ccopy(cmDATA_ini[i_new].e11,cmDATA[i].e11);
		Ccopy(cmDATA_ini[i_new].e12,cmDATA[i].e12);
		Ccopy(cmDATA_ini[i_new].e22,cmDATA[i].e22);
	}

/************************ Log2(N) point-to-point FT *************************/
/***************************** Initial step *********************************/

	iNolFFT=i_max/2;            // number of local FFT at the first step
	iStep=i_max/iNolFFT;
	for(i=0;i<iNolFFT;i++)
	{   i_new=iStep*i;
		//cmDATA_ini[i_new]  =cmDATA[i_new]+cmDATA[i_new+1];
		//cmDATA_ini[i_new+1]=cmDATA[i_new]-cmDATA[i_new+1];

		Cadd(cmDATA[i_new].e00,cmDATA[i_new+1].e00,cmDATA_ini[i_new].e00);
		Cadd(cmDATA[i_new].e01,cmDATA[i_new+1].e01,cmDATA_ini[i_new].e01);
		Cadd(cmDATA[i_new].e02,cmDATA[i_new+1].e02,cmDATA_ini[i_new].e02);
		Cadd(cmDATA[i_new].e11,cmDATA[i_new+1].e11,cmDATA_ini[i_new].e11);
		Cadd(cmDATA[i_new].e12,cmDATA[i_new+1].e12,cmDATA_ini[i_new].e12);
		Cadd(cmDATA[i_new].e22,cmDATA[i_new+1].e22,cmDATA_ini[i_new].e22);

		Csub(cmDATA[i_new].e00,cmDATA[i_new+1].e00,cmDATA_ini[i_new+1].e00);
		Csub(cmDATA[i_new].e01,cmDATA[i_new+1].e01,cmDATA_ini[i_new+1].e01);
		Csub(cmDATA[i_new].e02,cmDATA[i_new+1].e02,cmDATA_ini[i_new+1].e02);
		Csub(cmDATA[i_new].e11,cmDATA[i_new+1].e11,cmDATA_ini[i_new+1].e11);
		Csub(cmDATA[i_new].e12,cmDATA[i_new+1].e12,cmDATA_ini[i_new+1].e12);
		Csub(cmDATA[i_new].e22,cmDATA[i_new+1].e22,cmDATA_ini[i_new+1].e22);
	}

/***************************** Others steps *********************************/

	for(l=1;l<lim_NoB;l++)
	{   iNolFFT/=2;             // number of local FFT (At begin iNolFFT=i_max/2)
		iStep=i_max/iNolFFT;    // number of point in local FFT
		for(j=0;j<iNolFFT;j++)
		{   j_new=iStep*j;
			for(i=0;i<(iStep/2);i++)
			{   i_new=j_new+i;  k=i*i_max/iStep;  m=i_new+(iStep/2);

				//cm0=cW[k]*cmDATA_ini[m];
				//cmDATA[i_new]=cmDATA_ini[i_new]+cm0;
				//cmDATA[m]=cmDATA_ini[i_new]-cm0;

				Cmul(cW[k],cmDATA_ini[m].e00,cm0.e00);
				Cmul(cW[k],cmDATA_ini[m].e01,cm0.e01);
				Cmul(cW[k],cmDATA_ini[m].e02,cm0.e02);
				Cmul(cW[k],cmDATA_ini[m].e11,cm0.e11);
				Cmul(cW[k],cmDATA_ini[m].e12,cm0.e12);
				Cmul(cW[k],cmDATA_ini[m].e22,cm0.e22);

				Cadd(cm0.e00,cmDATA_ini[i_new].e00,cmDATA[i_new].e00);
				Cadd(cm0.e01,cmDATA_ini[i_new].e01,cmDATA[i_new].e01);
				Cadd(cm0.e02,cmDATA_ini[i_new].e02,cmDATA[i_new].e02);
				Cadd(cm0.e11,cmDATA_ini[i_new].e11,cmDATA[i_new].e11);
				Cadd(cm0.e12,cmDATA_ini[i_new].e12,cmDATA[i_new].e12);
				Cadd(cm0.e22,cmDATA_ini[i_new].e22,cmDATA[i_new].e22);

				Csub(cmDATA_ini[i_new].e00,cm0.e00,cmDATA[m].e00);
				Csub(cmDATA_ini[i_new].e01,cm0.e01,cmDATA[m].e01);
				Csub(cmDATA_ini[i_new].e02,cm0.e02,cmDATA[m].e02);
				Csub(cmDATA_ini[i_new].e11,cm0.e11,cmDATA[m].e11);
				Csub(cmDATA_ini[i_new].e12,cm0.e12,cmDATA[m].e12);
				Csub(cmDATA_ini[i_new].e22,cm0.e22,cmDATA[m].e22);
			}
		}

		for(i=0;i<i_max;i++)
		{   //cmDATA_ini[i]=cmDATA[i];

			Ccopy(cmDATA[i].e00,cmDATA_ini[i].e00);
			Ccopy(cmDATA[i].e01,cmDATA_ini[i].e01);
			Ccopy(cmDATA[i].e02,cmDATA_ini[i].e02);
			Ccopy(cmDATA[i].e11,cmDATA_ini[i].e11);
			Ccopy(cmDATA[i].e12,cmDATA_ini[i].e12);
			Ccopy(cmDATA[i].e22,cmDATA_ini[i].e22);
		}
	}
	return;
}

void FFT(int lim_NoB, int i_max, cplx * cW, cvector * cvDATA_ini,
		 cvector * cvDATA)
{   int i,j,l,k,m,i_new,j_new,iNolFFT,iStep;
	cvector cv0;

/************************* The sorting of initial data **********************/

	for(i=0;i<i_max;i++)
	{   Change_binary_number(i,lim_NoB,&i_new);
		//cvDATA[i]=cvDATA_ini[i_new];
		Ccopy(cvDATA_ini[i_new].x,cvDATA[i].x);
		Ccopy(cvDATA_ini[i_new].y,cvDATA[i].y);
		Ccopy(cvDATA_ini[i_new].z,cvDATA[i].z);
	}

/************************ Log2(N) point-to-point FT *************************/
/***************************** Initial step *********************************/

	iNolFFT=i_max/2;            // number of local FFT at the first step
	iStep=i_max/iNolFFT;
	for(i=0;i<iNolFFT;i++)
	{   i_new=iStep*i;

		//cvDATA_ini[i_new]  =cvDATA[i_new]+cvDATA[i_new+1];
		//cvDATA_ini[i_new+1]=cvDATA[i_new]-cvDATA[i_new+1];

		Cadd(cvDATA[i_new].x,cvDATA[i_new+1].x,cvDATA_ini[i_new].x);
		Cadd(cvDATA[i_new].y,cvDATA[i_new+1].y,cvDATA_ini[i_new].y);
		Cadd(cvDATA[i_new].z,cvDATA[i_new+1].z,cvDATA_ini[i_new].z);

		Csub(cvDATA[i_new].x,cvDATA[i_new+1].x,cvDATA_ini[i_new+1].x);
		Csub(cvDATA[i_new].y,cvDATA[i_new+1].y,cvDATA_ini[i_new+1].y);
		Csub(cvDATA[i_new].z,cvDATA[i_new+1].z,cvDATA_ini[i_new+1].z);
	}

/***************************** Others steps *********************************/

	for(l=1;l<lim_NoB;l++)
	{   iNolFFT/=2;             // number of local FFT (At begin iNolFFT=i_max/2)
		iStep=i_max/iNolFFT;    // number of point in local FFT
		for(j=0;j<iNolFFT;j++)
		{   j_new=iStep*j;
			for(i=0;i<(iStep/2);i++)
			{   i_new=j_new+i;  k=i*i_max/iStep;  m=i_new+(iStep/2);

				//cv0=cW[k]*cvDATA_ini[m];
				//cvDATA[i_new]=cvDATA_ini[i_new]+cv0;
				//cvDATA[m]=cvDATA_ini[i_new]-cv0;

				Cmul(cW[k],cvDATA_ini[m].x,cv0.x);
				Cmul(cW[k],cvDATA_ini[m].y,cv0.y);
				Cmul(cW[k],cvDATA_ini[m].z,cv0.z);

				Cadd(cvDATA_ini[i_new].x,cv0.x,cvDATA[i_new].x);
				Cadd(cvDATA_ini[i_new].y,cv0.y,cvDATA[i_new].y);
				Cadd(cvDATA_ini[i_new].z,cv0.z,cvDATA[i_new].z);

				Csub(cvDATA_ini[i_new].x,cv0.x,cvDATA[m].x);
				Csub(cvDATA_ini[i_new].y,cv0.y,cvDATA[m].y);
				Csub(cvDATA_ini[i_new].z,cv0.z,cvDATA[m].z);
			}
		}

		for(i=0;i<i_max;i++)
		{   //cvDATA_ini[i]=cvDATA[i];
			Ccopy(cvDATA[i].x,cvDATA_ini[i].x);
			Ccopy(cvDATA[i].y,cvDATA_ini[i].y);
			Ccopy(cvDATA[i].z,cvDATA_ini[i].z);
		}
	}
	return;
}

/****************************************************************************/
/**************************** Inverse FFT base ******************************/
/****************************************************************************/

void iFFT(int lim_NoB, int i_max, cplx * cW, cmatrix * cmDATA_ini,
		  cmatrix * cmDATA)
{   int i;

	for(i=0;i<i_max;i++)
	{   cmDATA_ini[i].e00=conj(cmDATA_ini[i].e00);
		cmDATA_ini[i].e01=conj(cmDATA_ini[i].e01);
		cmDATA_ini[i].e02=conj(cmDATA_ini[i].e02);
		cmDATA_ini[i].e11=conj(cmDATA_ini[i].e11);
		cmDATA_ini[i].e12=conj(cmDATA_ini[i].e12);
		cmDATA_ini[i].e22=conj(cmDATA_ini[i].e22);
	}

	FFT(lim_NoB,i_max,cW,cmDATA_ini,cmDATA);

	for(i=0;i<i_max;i++)
	{   cmDATA[i].e00=conj(cmDATA[i].e00);
		cmDATA[i].e01=conj(cmDATA[i].e01);
		cmDATA[i].e02=conj(cmDATA[i].e02);
		cmDATA[i].e11=conj(cmDATA[i].e11);
		cmDATA[i].e12=conj(cmDATA[i].e12);
		cmDATA[i].e22=conj(cmDATA[i].e22);
	}

	return;
}

void iFFT(int lim_NoB, int i_max, cplx * cW, cvector * cvDATA_ini,
		  cvector * cvDATA)
{   int i;

	for(i=0;i<i_max;i++)
	{   cvDATA_ini[i].x=conj(cvDATA_ini[i].x);
		cvDATA_ini[i].y=conj(cvDATA_ini[i].y);
		cvDATA_ini[i].z=conj(cvDATA_ini[i].z);
	}

	FFT(lim_NoB,i_max,cW,cvDATA_ini,cvDATA);

	for(i=0;i<i_max;i++)
	{   cvDATA[i].x=conj(cvDATA[i].x);
		cvDATA[i].y=conj(cvDATA[i].y);
		cvDATA[i].z=conj(cvDATA[i].z);
	}

	return;
}

/****************************************************************************/
/********************************** FFT 3D **********************************/
/****************************************************************************/

void FFT_3D(int i_max, int j_max, int l_max, cplx * cW_i, cplx * cW_j,
			cplx * cW_l, cmatrix *** cmDATA)
{   int i,j,l,lim_NoB;
	real rNorma=1./(i_max*j_max*l_max);
	cmatrix * cmWork,* cmResult;

/*********************************** z-FFT **********************************/

	cmWork=(cmatrix *)malloc(l_max*sizeof(cmatrix));
	if(cmWork==NULL)
	{   printf("\n No enough memory for array cmWork[l_max]! ");
		exit(1);
	}

	cmResult=(cmatrix *)malloc(l_max*sizeof(cmatrix));
	if(cmResult==NULL)
	{   printf("\n No enough memory for array cmResult[l_max]! ");
		exit(1);
	}

	Get_number_bits((l_max-1),&lim_NoB);

	for(i=0;i<i_max;i++)
		for(j=0;j<j_max;j++)
		{   for(l=0;l<l_max;l++)
				cmWork[l]=cmDATA[i][j][l];

			FFT(lim_NoB,l_max,cW_l,cmWork,cmResult);

			for(l=0;l<l_max;l++)
				cmDATA[i][j][l]=cmResult[l];
		}

	free(cmWork);   free(cmResult);

/*********************************** y-FFT **********************************/

	cmWork=(cmatrix *)malloc(j_max*sizeof(cmatrix));
	if(cmWork==NULL)
	{   printf("\n No enough memory for array cmWork[j_max]! ");
		exit(1);
	}

	cmResult=(cmatrix *)malloc(j_max*sizeof(cmatrix));
	if(cmResult==NULL)
	{   printf("\n No enough memory for array cmResult[j_max]! ");
		exit(1);
	}

	Get_number_bits((j_max-1),&lim_NoB);

	for(i=0;i<i_max;i++)
		for(l=0;l<l_max;l++)
		{   for(j=0;j<j_max;j++)
				cmWork[j]=cmDATA[i][j][l];

			FFT(lim_NoB,j_max,cW_j,cmWork,cmResult);

			for(j=0;j<j_max;j++)
				cmDATA[i][j][l]=cmResult[j];
		}

	free(cmWork);   free(cmResult);

/*********************************** x-FFT **********************************/

	cmWork=(cmatrix *)malloc(i_max*sizeof(cmatrix));
	if(cmWork==NULL)
	{   printf("\n No enough memory for array cmWork[i_max]! ");
		exit(1);
	}

	cmResult=(cmatrix *)malloc(i_max*sizeof(cmatrix));
	if(cmResult==NULL)
	{   printf("\n No enough memory for array cmResult[i_max]! ");
		exit(1);
	}

	Get_number_bits((i_max-1),&lim_NoB);

	for(j=0;j<j_max;j++)
		for(l=0;l<l_max;l++)
		{   for(i=0;i<i_max;i++)
				cmWork[i]=cmDATA[i][j][l];

			FFT(lim_NoB,i_max,cW_i,cmWork,cmResult);

			for(i=0;i<i_max;i++)
				cmDATA[i][j][l]=cmResult[i];
		}

	free(cmWork);   free(cmResult);

	for(i=0;i<i_max;i++)
		for(j=0;j<j_max;j++)
			for(l=0;l<l_max;l++)
				cmDATA[i][j][l]=rNorma*cmDATA[i][j][l];

	return;
}

void FFT_3D(int i_max, int j_max, int l_max, cplx * cW_2_i, cplx * cW_2_j,
			cplx * cW_2_l, cplx * cW_i, cplx * cW_j, cplx * cW_l,
			cvector *** cvDATA)
{   int i,j,l,lim_NoB,i_max_2=i_max/2,j_max_2=j_max/2,l_max_2=l_max/2;
	cvector * cvWork,* cvWork2,* cvResult;

/*********************************** z-FFT **********************************/

	cvWork=(cvector *)malloc(l_max_2*sizeof(cvector));
	if(cvWork==NULL)
	{   printf("\n No enough memory for array cvWork[l_max]! ");
		exit(1);
	}

	cvWork2=(cvector *)malloc(l_max_2*sizeof(cvector));
	if(cvWork2==NULL)
	{   printf("\n No enough memory for array cvWork2[l_max]! ");
		exit(1);
	}

	cvResult=(cvector *)malloc(l_max_2*sizeof(cvector));
	if(cvResult==NULL)
	{   printf("\n No enough memory for array cvResult[l_max]! ");
		exit(1);
	}

	Get_number_bits((l_max_2-1),&lim_NoB);

	for(i=0;i<i_max_2;i++)
		for(j=0;j<j_max_2;j++)
		{   for(l=0;l<l_max_2;l++)
				cvWork[l]=cvDATA[i][j][l];

			FFT(lim_NoB,l_max_2,cW_2_l,cvWork,cvResult);

			for(l=0;l<l_max_2;l++)
				cvWork2[l]=cvResult[l];

			for(l=0;l<l_max_2;l++)
			{   Cmul(cW_l[l],cvDATA[i][j][l].x,cvWork[l].x);
				Cmul(cW_l[l],cvDATA[i][j][l].y,cvWork[l].y);
				Cmul(cW_l[l],cvDATA[i][j][l].z,cvWork[l].z);
			}

			FFT(lim_NoB,l_max_2,cW_2_l,cvWork,cvResult);

			for(l=0;l<l_max_2;l++)
			{   cvDATA[i][j][2*l]=cvWork2[l];
				cvDATA[i][j][2*l+1]=cvResult[l];
			}
		}

	free(cvWork);   free(cvWork2);  free(cvResult);

/*********************************** y-FFT **********************************/

	cvWork=(cvector *)malloc(j_max_2*sizeof(cvector));
	if(cvWork==NULL)
	{   printf("\n No enough memory for array cvWork[j_max]! ");
		exit(1);
	}

	cvWork2=(cvector *)malloc(j_max_2*sizeof(cvector));
	if(cvWork2==NULL)
	{   printf("\n No enough memory for array cvWork2[j_max]! ");
		exit(1);
	}

	cvResult=(cvector *)malloc(j_max_2*sizeof(cvector));
	if(cvResult==NULL)
	{   printf("\n No enough memory for array cvResult[j_max]! ");
		exit(1);
	}

	Get_number_bits((j_max_2-1),&lim_NoB);

	for(i=0;i<i_max_2;i++)
		for(l=0;l<l_max;l++)
		{   for(j=0;j<j_max_2;j++)
				cvWork[j]=cvDATA[i][j][l];

			FFT(lim_NoB,j_max_2,cW_2_j,cvWork,cvResult);

			for(j=0;j<j_max_2;j++)
				cvWork2[j]=cvResult[j];

			for(j=0;j<j_max_2;j++)
			{   Cmul(cW_j[j],cvDATA[i][j][l].x,cvWork[j].x);
				Cmul(cW_j[j],cvDATA[i][j][l].y,cvWork[j].y);
				Cmul(cW_j[j],cvDATA[i][j][l].z,cvWork[j].z);
			}

			FFT(lim_NoB,j_max_2,cW_2_j,cvWork,cvResult);

			for(j=0;j<j_max_2;j++)
			{   cvDATA[i][2*j][l]=cvWork2[j];
				cvDATA[i][2*j+1][l]=cvResult[j];
			}
		}

	free(cvWork);   free(cvWork2);  free(cvResult);

/*********************************** x-FFT **********************************/

	cvWork=(cvector *)malloc(i_max_2*sizeof(cvector));
	if(cvWork==NULL)
	{   printf("\n No enough memory for array cvWork[i_max]! ");
		exit(1);
	}

	cvWork2=(cvector *)malloc(i_max_2*sizeof(cvector));
	if(cvWork2==NULL)
	{   printf("\n No enough memory for array cvWork2[i_max]! ");
		exit(1);
	}

	cvResult=(cvector *)malloc(i_max_2*sizeof(cvector));
	if(cvResult==NULL)
	{   printf("\n No enough memory for array cvResult[i_max]! ");
		exit(1);
	}

	Get_number_bits((i_max_2-1),&lim_NoB);

	for(j=0;j<j_max;j++)
		for(l=0;l<l_max;l++)
		{   for(i=0;i<i_max_2;i++)
				cvWork[i]=cvDATA[i][j][l];

			FFT(lim_NoB,i_max_2,cW_2_i,cvWork,cvResult);

			for(i=0;i<i_max_2;i++)
				cvWork2[i]=cvResult[i];

			for(i=0;i<i_max_2;i++)
			{   Cmul(cW_i[i],cvDATA[i][j][l].x,cvWork[i].x);
				Cmul(cW_i[i],cvDATA[i][j][l].y,cvWork[i].y);
				Cmul(cW_i[i],cvDATA[i][j][l].z,cvWork[i].z);
			}

			FFT(lim_NoB,i_max_2,cW_2_i,cvWork,cvResult);

			for(i=0;i<i_max_2;i++)
			{   cvDATA[2*i][j][l]=cvWork2[i];
				cvDATA[2*i+1][j][l]=cvResult[i];
			}
		}

	free(cvWork);   free(cvWork2);  free(cvResult);

	return;
}

/****************************************************************************/
/****************************** Inverse FFT 3D ******************************/
/****************************************************************************/

void iFFT_3D(int i_max, int j_max, int l_max, cplx * cW_i, cplx * cW_j,
			 cplx * cW_l, cvector *** cvDATA)
{   int i,j,l,lim_NoB,i_max_2=i_max/2,j_max_2=j_max/2;
	cvector * cvWork,* cvResult;

/*********************************** x-FFT **********************************/

	cvWork=(cvector *)malloc(i_max*sizeof(cvector));
	if(cvWork==NULL)
	{   printf("\n No enough memory for array cvWork[i_max]! ");
		exit(1);
	}

	cvResult=(cvector *)malloc(i_max*sizeof(cvector));
	if(cvResult==NULL)
	{   printf("\n No enough memory for array cvResult[i_max]! ");
		exit(1);
	}

	Get_number_bits((i_max-1),&lim_NoB);

	for(j=0;j<j_max;j++)
		for(l=0;l<l_max;l++)
		{   for(i=0;i<i_max;i++)
				cvWork[i]=cvDATA[i][j][l];

			iFFT(lim_NoB,i_max,cW_i,cvWork,cvResult);

			for(i=0;i<i_max;i++)
				cvDATA[i][j][l]=cvResult[i];
		}

	free(cvWork);   free(cvResult);

/*********************************** y-FFT **********************************/

	cvWork=(cvector *)malloc(j_max*sizeof(cvector));
	if(cvWork==NULL)
	{   printf("\n No enough memory for array cvWork[j_max]! ");
		exit(1);
	}

	cvResult=(cvector *)malloc(j_max*sizeof(cvector));
	if(cvResult==NULL)
	{   printf("\n No enough memory for array cvResult[j_max]! ");
		exit(1);
	}

	Get_number_bits((j_max-1),&lim_NoB);

	for(i=0;i<i_max_2;i++)
		for(l=0;l<l_max;l++)
		{   for(j=0;j<j_max;j++)
				cvWork[j]=cvDATA[i][j][l];

			iFFT(lim_NoB,j_max,cW_j,cvWork,cvResult);

			for(j=0;j<j_max;j++)
				cvDATA[i][j][l]=cvResult[j];
		}

	free(cvWork);   free(cvResult);

/*********************************** z-FFT **********************************/

	cvWork=(cvector *)malloc(l_max*sizeof(cvector));
	if(cvWork==NULL)
	{   printf("\n No enough memory for array cvWork[l_max]! ");
		exit(1);
	}

	cvResult=(cvector *)malloc(l_max*sizeof(cvector));
	if(cvResult==NULL)
	{   printf("\n No enough memory for array cvResult[l_max]! ");
		exit(1);
	}

	Get_number_bits((l_max-1),&lim_NoB);

	for(i=0;i<i_max_2;i++)
		for(j=0;j<j_max_2;j++)
		{   for(l=0;l<l_max;l++)
				cvWork[l]=cvDATA[i][j][l];

			iFFT(lim_NoB,l_max,cW_l,cvWork,cvResult);

			for(l=0;l<l_max;l++)
				cvDATA[i][j][l]=cvResult[l];
		}

	free(cvWork);   free(cvResult);

	return;
}

/***************** FFT functions for Rader-Brenner method *******************/

/****************************************************************************/
/******************* Calculation the exp(-i((2*PI*j)/(N))) ******************/
/****************************************************************************/

void FourierPhase(int i_max, real ** Im_cW, cplx ** cW2)
{   int i,j,N_2=i_max/2;
	real rF=2.*PI/i_max;
	cplx c0;

	for(j=0;j<N_2;j++)
	{   Im_cW[0][j]=0.;
		cW2[0][j].r=1.; cW2[0][j].i=0.;
	}

	for(i=1;i<N_2;i++)
		for(j=i;j<N_2;j++)
		{   Cexp_im((rF*(real)(i*j)),c0);
			Cmul(c0,c0,cW2[i][j]);
			Im_cW[i][j]=.5/c0.i;
		}

	for(i=1;i<N_2;i++)
		for(j=0;j<i;j++)
		{   Im_cW[i][j]=Im_cW[j][i];
			cW2[i][j]=cW2[j][i];
		}

	return;
}

/****************************************************************************/
/********************************* FFT base *********************************/
/****************************************************************************/

void FFT(int i_max, real ** Im_cW, cplx ** cW2, cmatrix * cmDATA_ini,
		 cmatrix * cmDATA, cmatrix * cA, cmatrix * cB, cmatrix * cC)
{   int i,j,N_2=i_max/2;
	cmatrix c0,c1;

	Cadd(cmDATA_ini[0].e00,cmDATA_ini[N_2].e00,cB[0].e00);
	Cadd(cmDATA_ini[0].e01,cmDATA_ini[N_2].e01,cB[0].e01);
	Cadd(cmDATA_ini[0].e02,cmDATA_ini[N_2].e02,cB[0].e02);
	Cadd(cmDATA_ini[0].e11,cmDATA_ini[N_2].e11,cB[0].e11);
	Cadd(cmDATA_ini[0].e12,cmDATA_ini[N_2].e12,cB[0].e12);
	Cadd(cmDATA_ini[0].e22,cmDATA_ini[N_2].e22,cB[0].e22);

	cA[0].e00.r=0.; cA[0].e00.i=0.;
	cA[0].e01.r=0.; cA[0].e01.i=0.;
	cA[0].e02.r=0.; cA[0].e02.i=0.;
	cA[0].e11.r=0.; cA[0].e11.i=0.;
	cA[0].e12.r=0.; cA[0].e12.i=0.;
	cA[0].e22.r=0.; cA[0].e22.i=0.;

	for(i=1;i<N_2;i++)
	{   Cadd(cmDATA_ini[i].e00,cmDATA_ini[i+N_2].e00,cB[i].e00);
		Cadd(cmDATA_ini[i].e01,cmDATA_ini[i+N_2].e01,cB[i].e01);
		Cadd(cmDATA_ini[i].e02,cmDATA_ini[i+N_2].e02,cB[i].e02);
		Cadd(cmDATA_ini[i].e11,cmDATA_ini[i+N_2].e11,cB[i].e11);
		Cadd(cmDATA_ini[i].e12,cmDATA_ini[i+N_2].e12,cB[i].e12);
		Cadd(cmDATA_ini[i].e22,cmDATA_ini[i+N_2].e22,cB[i].e22);

		Csub(cmDATA_ini[i].e00,cmDATA_ini[i+N_2].e00,c0.e00);
		Csub(cmDATA_ini[i].e01,cmDATA_ini[i+N_2].e01,c0.e01);
		Csub(cmDATA_ini[i].e02,cmDATA_ini[i+N_2].e02,c0.e02);
		Csub(cmDATA_ini[i].e11,cmDATA_ini[i+N_2].e11,c0.e11);
		Csub(cmDATA_ini[i].e12,cmDATA_ini[i+N_2].e12,c0.e12);
		Csub(cmDATA_ini[i].e22,cmDATA_ini[i+N_2].e22,c0.e22);

		cA[i].e00.r=c0.e00.i*Im_cW[i][1];   cA[i].e00.i=-c0.e00.r*Im_cW[i][1];
		cA[i].e01.r=c0.e01.i*Im_cW[i][1];   cA[i].e01.i=-c0.e01.r*Im_cW[i][1];
		cA[i].e02.r=c0.e02.i*Im_cW[i][1];   cA[i].e02.i=-c0.e02.r*Im_cW[i][1];
		cA[i].e11.r=c0.e11.i*Im_cW[i][1];   cA[i].e11.i=-c0.e11.r*Im_cW[i][1];
		cA[i].e12.r=c0.e12.i*Im_cW[i][1];   cA[i].e12.i=-c0.e12.r*Im_cW[i][1];
		cA[i].e22.r=c0.e22.i*Im_cW[i][1];   cA[i].e22.i=-c0.e22.r*Im_cW[i][1];
	}

	for(i=0;i<N_2;i++)
	{   cmDATA[2*i].e00.r=0.;   cmDATA[2*i].e00.i=0.;
		cmDATA[2*i].e01.r=0.;   cmDATA[2*i].e01.i=0.;
		cmDATA[2*i].e02.r=0.;   cmDATA[2*i].e02.i=0.;
		cmDATA[2*i].e11.r=0.;   cmDATA[2*i].e11.i=0.;
		cmDATA[2*i].e12.r=0.;   cmDATA[2*i].e12.i=0.;
		cmDATA[2*i].e22.r=0.;   cmDATA[2*i].e22.i=0.;

		cC[i].e00.r=0.;         cC[i].e00.i=0.;
		cC[i].e01.r=0.;         cC[i].e01.i=0.;
		cC[i].e02.r=0.;         cC[i].e02.i=0.;
		cC[i].e11.r=0.;         cC[i].e11.i=0.;
		cC[i].e12.r=0.;         cC[i].e12.i=0.;
		cC[i].e22.r=0.;         cC[i].e22.i=0.;

		for(j=0;j<N_2;j++)
		{   Cmul(cB[j].e00,cW2[i][j],c0.e00);
			Cmul(cB[j].e01,cW2[i][j],c0.e01);
			Cmul(cB[j].e02,cW2[i][j],c0.e02);
			Cmul(cB[j].e11,cW2[i][j],c0.e11);
			Cmul(cB[j].e12,cW2[i][j],c0.e12);
			Cmul(cB[j].e22,cW2[i][j],c0.e22);

			Cadd(cmDATA[2*i].e00,c0.e00,cmDATA[2*i].e00);
			Cadd(cmDATA[2*i].e01,c0.e01,cmDATA[2*i].e01);
			Cadd(cmDATA[2*i].e02,c0.e02,cmDATA[2*i].e02);
			Cadd(cmDATA[2*i].e11,c0.e11,cmDATA[2*i].e11);
			Cadd(cmDATA[2*i].e12,c0.e12,cmDATA[2*i].e12);
			Cadd(cmDATA[2*i].e22,c0.e22,cmDATA[2*i].e22);

			Cmul(cA[j].e00,cW2[i][j],c0.e00);
			Cmul(cA[j].e01,cW2[i][j],c0.e01);
			Cmul(cA[j].e02,cW2[i][j],c0.e02);
			Cmul(cA[j].e11,cW2[i][j],c0.e11);
			Cmul(cA[j].e12,cW2[i][j],c0.e12);
			Cmul(cA[j].e22,cW2[i][j],c0.e22);

			Cadd(cC[i].e00,c0.e00,cC[i].e00);
			Cadd(cC[i].e01,c0.e01,cC[i].e01);
			Cadd(cC[i].e02,c0.e02,cC[i].e02);
			Cadd(cC[i].e11,c0.e11,cC[i].e11);
			Cadd(cC[i].e12,c0.e12,cC[i].e12);
			Cadd(cC[i].e22,c0.e22,cC[i].e22);
		}
	}

	cC[N_2]=cC[0];

	Csub(cmDATA_ini[0].e00,cmDATA_ini[N_2].e00,c0.e00);
	Csub(cmDATA_ini[0].e01,cmDATA_ini[N_2].e01,c0.e01);
	Csub(cmDATA_ini[0].e02,cmDATA_ini[N_2].e02,c0.e02);
	Csub(cmDATA_ini[0].e11,cmDATA_ini[N_2].e11,c0.e11);
	Csub(cmDATA_ini[0].e12,cmDATA_ini[N_2].e12,c0.e12);
	Csub(cmDATA_ini[0].e22,cmDATA_ini[N_2].e22,c0.e22);

	for(i=0;i<N_2;i++)
	{   Csub(cC[i+1].e00,cC[i].e00,c1.e00);
		Csub(cC[i+1].e01,cC[i].e01,c1.e01);
		Csub(cC[i+1].e02,cC[i].e02,c1.e02);
		Csub(cC[i+1].e11,cC[i].e11,c1.e11);
		Csub(cC[i+1].e12,cC[i].e12,c1.e12);
		Csub(cC[i+1].e22,cC[i].e22,c1.e22);

		Cadd(c1.e00,c0.e00,cmDATA[2*i+1].e00);
		Cadd(c1.e01,c0.e01,cmDATA[2*i+1].e01);
		Cadd(c1.e02,c0.e02,cmDATA[2*i+1].e02);
		Cadd(c1.e11,c0.e11,cmDATA[2*i+1].e11);
		Cadd(c1.e12,c0.e12,cmDATA[2*i+1].e12);
		Cadd(c1.e22,c0.e22,cmDATA[2*i+1].e22);
	}

	return;
}

void FFT(int i_max, real ** Im_cW, cplx ** cW2, cvector * cvDATA_ini,
		 cvector * cvDATA, cvector * cA, cvector * cB, cvector * cC)
{   int i,j,N_2=i_max/2;
	cvector c0,c1;

	Cadd(cvDATA_ini[0].x,cvDATA_ini[N_2].x,cB[0].x);
	Cadd(cvDATA_ini[0].y,cvDATA_ini[N_2].y,cB[0].y);
	Cadd(cvDATA_ini[0].z,cvDATA_ini[N_2].z,cB[0].z);

	cA[0].x.r=0.;   cA[0].x.i=0.;
	cA[0].y.r=0.;   cA[0].y.i=0.;
	cA[0].z.r=0.;   cA[0].z.i=0.;

	for(i=1;i<N_2;i++)
	{   Cadd(cvDATA_ini[i].x,cvDATA_ini[i+N_2].x,cB[i].x);
		Cadd(cvDATA_ini[i].y,cvDATA_ini[i+N_2].y,cB[i].y);
		Cadd(cvDATA_ini[i].z,cvDATA_ini[i+N_2].z,cB[i].z);

		Csub(cvDATA_ini[i].x,cvDATA_ini[i+N_2].x,c0.x);
		Csub(cvDATA_ini[i].y,cvDATA_ini[i+N_2].y,c0.y);
		Csub(cvDATA_ini[i].z,cvDATA_ini[i+N_2].z,c0.z);

		cA[i].x.r=c0.x.i*Im_cW[i][1];   cA[i].x.i=-c0.x.r*Im_cW[i][1];
		cA[i].y.r=c0.y.i*Im_cW[i][1];   cA[i].y.i=-c0.y.r*Im_cW[i][1];
		cA[i].z.r=c0.z.i*Im_cW[i][1];   cA[i].z.i=-c0.z.r*Im_cW[i][1];
	}

	for(i=0;i<N_2;i++)
	{   cvDATA[2*i].x.r=0.; cvDATA[2*i].x.i=0.;
		cvDATA[2*i].y.r=0.; cvDATA[2*i].y.i=0.;
		cvDATA[2*i].z.r=0.; cvDATA[2*i].z.i=0.;

		cC[i].x.r=0.;       cC[i].x.i=0.;
		cC[i].y.r=0.;       cC[i].y.i=0.;
		cC[i].z.r=0.;       cC[i].z.i=0.;

		for(j=0;j<N_2;j++)
		{   Cmul(cB[j].x,cW2[i][j],c0.x);
			Cmul(cB[j].y,cW2[i][j],c0.y);
			Cmul(cB[j].z,cW2[i][j],c0.z);

			Cadd(cvDATA[2*i].x,c0.x,cvDATA[2*i].x);
			Cadd(cvDATA[2*i].y,c0.y,cvDATA[2*i].y);
			Cadd(cvDATA[2*i].z,c0.z,cvDATA[2*i].z);

			Cmul(cA[j].x,cW2[i][j],c0.x);
			Cmul(cA[j].y,cW2[i][j],c0.y);
			Cmul(cA[j].z,cW2[i][j],c0.z);

			Cadd(cC[i].x,c0.x,cC[i].x);
			Cadd(cC[i].y,c0.y,cC[i].y);
			Cadd(cC[i].z,c0.z,cC[i].z);
		}
	}

	cC[N_2]=cC[0];

	Csub(cvDATA_ini[0].x,cvDATA_ini[N_2].x,c0.x);
	Csub(cvDATA_ini[0].y,cvDATA_ini[N_2].y,c0.y);
	Csub(cvDATA_ini[0].z,cvDATA_ini[N_2].z,c0.z);

	for(i=0;i<N_2;i++)
	{   Csub(cC[i+1].x,cC[i].x,c1.x);
		Csub(cC[i+1].y,cC[i].y,c1.y);
		Csub(cC[i+1].z,cC[i].z,c1.z);

		Cadd(c1.x,c0.x,cvDATA[2*i+1].x);
		Cadd(c1.y,c0.y,cvDATA[2*i+1].y);
		Cadd(c1.z,c0.z,cvDATA[2*i+1].z);
	}

	return;
}

/****************************************************************************/
/**************************** Inverse FFT base ******************************/
/****************************************************************************/

void iFFT(int i_max, real ** Im_cW, cplx ** cW2, cvector * cvDATA_ini,
		  cvector * cvDATA, cvector * cA, cvector * cB, cvector * cC)
{   int i;

	for(i=0;i<i_max;i++)
	{   cvDATA_ini[i].x=conj(cvDATA_ini[i].x);
		cvDATA_ini[i].y=conj(cvDATA_ini[i].y);
		cvDATA_ini[i].z=conj(cvDATA_ini[i].z);
	}

	FFT(i_max,Im_cW,cW2,cvDATA_ini,cvDATA,cA,cB,cC);

	for(i=0;i<i_max;i++)
	{   cvDATA[i].x=conj(cvDATA[i].x);
		cvDATA[i].y=conj(cvDATA[i].y);
		cvDATA[i].z=conj(cvDATA[i].z);
	}

	return;
}

/****************************************************************************/
/********************************** FFT 3D **********************************/
/****************************************************************************/

void FFT_3D(int i_max, int j_max, int l_max, real ** Im_cW_i, real ** Im_cW_j,
			real ** Im_cW_l, cplx ** cW2_i, cplx ** cW2_j, cplx ** cW2_l,
			cmatrix *** cmDATA)
{   int i,j,l,N;
	real rNorma=1./(i_max*j_max*l_max);
	cmatrix * cmWork,* cmResult,* cA,* cB,* cC;

	i=i_max;
	if(i<j_max) i=j_max;
	if(i<l_max) i=l_max;
	N=i/2;

	cA=(cmatrix *)malloc(N*sizeof(cmatrix));
	cB=(cmatrix *)malloc(N*sizeof(cmatrix));
	cC=(cmatrix *)malloc((N+1)*sizeof(cmatrix));
	cmWork=(cmatrix *)malloc(i*sizeof(cmatrix));
	cmResult=(cmatrix *)malloc(i*sizeof(cmatrix));
	if(cA==NULL||cB==NULL||cC==NULL||cmWork==NULL||cmResult==NULL)
	{   printf("\n No memory for array cA, cB, cC, cmWork[] or cmResult[]!");
		exit(1);
	}

/*********************************** z-FFT **********************************/

	for(i=0;i<i_max;i++)
		for(j=0;j<j_max;j++)
		{   for(l=0;l<l_max;l++)
				cmWork[l]=cmDATA[i][j][l];

			FFT(l_max,Im_cW_l,cW2_l,cmWork,cmResult,cA,cB,cC);

			for(l=0;l<l_max;l++)
				cmDATA[i][j][l]=cmResult[l];
		}

/*********************************** y-FFT **********************************/

	for(i=0;i<i_max;i++)
		for(l=0;l<l_max;l++)
		{   for(j=0;j<j_max;j++)
				cmWork[j]=cmDATA[i][j][l];

			FFT(j_max,Im_cW_j,cW2_j,cmWork,cmResult,cA,cB,cC);

			for(j=0;j<j_max;j++)
				cmDATA[i][j][l]=cmResult[j];
		}

/*********************************** x-FFT **********************************/

	for(j=0;j<j_max;j++)
		for(l=0;l<l_max;l++)
		{   for(i=0;i<i_max;i++)
				cmWork[i]=cmDATA[i][j][l];

			FFT(i_max,Im_cW_i,cW2_i,cmWork,cmResult,cA,cB,cC);

			for(i=0;i<i_max;i++)
				cmDATA[i][j][l]=cmResult[i];
		}

	for(i=0;i<i_max;i++)
		for(j=0;j<j_max;j++)
			for(l=0;l<l_max;l++)
				cmDATA[i][j][l]=rNorma*cmDATA[i][j][l];

	free(cmWork);   free(cmResult);   free(cA);   free(cB);   free(cC);
	return;
}

void FFT_3D(int i_max, int j_max, int l_max, real ** Im_cW_i, real ** Im_cW_j,
			real ** Im_cW_l, cplx ** cW2_i, cplx ** cW2_j, cplx ** cW2_l,
			cvector *** cvDATA)
{   int i,j,l,N,i_max_2=i_max/2,j_max_2=j_max/2;
	cvector * cvWork,* cvResult,* cA,* cB,* cC;

	i=i_max;
	if(i<j_max) i=j_max;
	if(i<l_max) i=l_max;
	N=i/2;

	cA=(cvector *)malloc(N*sizeof(cvector));
	cB=(cvector *)malloc(N*sizeof(cvector));
	cC=(cvector *)malloc((N+1)*sizeof(cvector));
	cvWork=(cvector *)malloc(i*sizeof(cvector));
	cvResult=(cvector *)malloc(i*sizeof(cvector));
	if(cA==NULL||cB==NULL||cC==NULL||cvWork==NULL||cvResult==NULL)
	{   printf("\n No memory for array cA, cB, cC, cvWork[] or cvResult[]!");
		exit(1);
	}

/*********************************** z-FFT **********************************/

	for(i=0;i<i_max_2;i++)
		for(j=0;j<j_max_2;j++)
		{   for(l=0;l<l_max;l++)
				cvWork[l]=cvDATA[i][j][l];

			FFT(l_max,Im_cW_l,cW2_l,cvWork,cvResult,cA,cB,cC);

			for(l=0;l<l_max;l++)
				cvDATA[i][j][l]=cvResult[l];
		}

/*********************************** y-FFT **********************************/

	for(i=0;i<i_max_2;i++)
		for(l=0;l<l_max;l++)
		{   for(j=0;j<j_max;j++)
				cvWork[j]=cvDATA[i][j][l];

			FFT(j_max,Im_cW_j,cW2_j,cvWork,cvResult,cA,cB,cC);

			for(j=0;j<j_max;j++)
				cvDATA[i][j][l]=cvResult[j];
		}

/*********************************** x-FFT **********************************/

	for(j=0;j<j_max;j++)
		for(l=0;l<l_max;l++)
		{   for(i=0;i<i_max;i++)
				cvWork[i]=cvDATA[i][j][l];

			FFT(i_max,Im_cW_i,cW2_i,cvWork,cvResult,cA,cB,cC);

			for(i=0;i<i_max;i++)
				cvDATA[i][j][l]=cvResult[i];
		}

	free(cvWork);   free(cvResult);   free(cA);   free(cB);   free(cC);
	return;
}

/****************************************************************************/
/****************************** Inverse FFT 3D ******************************/
/****************************************************************************/

void iFFT_3D(int i_max, int j_max, int l_max, real ** Im_cW_i, real ** Im_cW_j,
			 real ** Im_cW_l, cplx ** cW2_i, cplx ** cW2_j, cplx ** cW2_l,
			 cvector *** cvDATA)
{   int i,j,l,N,i_max_2=i_max/2,j_max_2=j_max/2;
	cvector * cvWork,* cvResult,* cA,* cB,* cC;

	i=i_max;
	if(i<j_max) i=j_max;
	if(i<l_max) i=l_max;
	N=i/2;

	cA=(cvector *)malloc(N*sizeof(cvector));
	cB=(cvector *)malloc(N*sizeof(cvector));
	cC=(cvector *)malloc((N+1)*sizeof(cvector));
	cvWork=(cvector *)malloc(i*sizeof(cvector));
	cvResult=(cvector *)malloc(i*sizeof(cvector));
	if(cA==NULL||cB==NULL||cC==NULL||cvWork==NULL||cvResult==NULL)
	{   printf("\n No memory for array cA, cB, cC, cvWork[] or cvResult[]!");
		exit(1);
	}

/*********************************** x-FFT **********************************/

	for(j=0;j<j_max;j++)
		for(l=0;l<l_max;l++)
		{   for(i=0;i<i_max;i++)
				cvWork[i]=cvDATA[i][j][l];

			iFFT(i_max,Im_cW_i,cW2_i,cvWork,cvResult,cA,cB,cC);

			for(i=0;i<i_max;i++)
				cvDATA[i][j][l]=cvResult[i];
		}

/*********************************** y-FFT **********************************/

	for(i=0;i<i_max_2;i++)
		for(l=0;l<l_max;l++)
		{   for(j=0;j<j_max;j++)
				cvWork[j]=cvDATA[i][j][l];

			iFFT(j_max,Im_cW_j,cW2_j,cvWork,cvResult,cA,cB,cC);

			for(j=0;j<j_max;j++)
				cvDATA[i][j][l]=cvResult[j];
		}

/*********************************** z-FFT **********************************/

	for(i=0;i<i_max_2;i++)
		for(j=0;j<j_max_2;j++)
		{   for(l=0;l<l_max;l++)
				cvWork[l]=cvDATA[i][j][l];

			iFFT(l_max,Im_cW_l,cW2_l,cvWork,cvResult,cA,cB,cC);

			for(l=0;l<l_max;l++)
				cvDATA[i][j][l]=cvResult[l];
		}

	free(cvWork);   free(cvResult);   free(cA);   free(cB);   free(cC);
	return;
}
