
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "complex.hpp"
#include "cvector.hpp"
#include "function.hpp"

int CrossTwoWorms(cylinder cl0, cylinder cl1);
cylinder worm(real r1, real theta1, real phi1, real r2, real theta2, real phi2,
			  real length, real radius, real R2);
char SeedToWorm(cylinder clWorm, vector vPoint);
int PointToWorm(cylinder clWorm, vector vPoint, real radius);

/****************************************************************************/
/************************* Generation of particles **************************/
/****************************************************************************/

/******************* Generation of anthrax-like particles *******************/

int iGenerateAnthrax(int i_max, int j_max, int l_max, real radius_worm,
					 real length_worm, int iWormsNumber, cylinder * worms,
					 char *** cAlpha)
{   int i,j,l,m,iTmp=0;
	real R=(real)i_max/2.,R2=(R-radius_worm)*(R-radius_worm);
	real r0,r1,r2,theta1,theta2,phi1,phi2;
	vector vCenter={(real)((real)i_max/2.-.5),(real)((real)j_max/2.-.5),(real)((real)l_max/2.-.5)},v0,v1;
	cylinder cl0;

	for(i=0;i<i_max;i++)
		for(j=0;j<j_max;j++)
			for(l=0;l<l_max;l++)
				cAlpha[i][j][l]=0;

	l=0;    m=0;
	while(m<iWormsNumber)
	{   cl0.r=-1.;
		while(cl0.r<0.)
		{   r1=R*(real)dda_rand()/(real)DDA_RAND_MAX;
			r2=r1;

			phi1=2.*PI*(real)dda_rand()/(real)DDA_RAND_MAX;
			phi2=2.*PI*(real)dda_rand()/(real)DDA_RAND_MAX;

			r0=(2.*(real)dda_rand()/(real)DDA_RAND_MAX-1.);
			theta1=acos(r0);
			r0=(2.*(real)dda_rand()/(real)DDA_RAND_MAX-1.);
			theta2=acos(r0);

			cl0=worm(r1,theta1,phi1,r2,theta2,phi2,length_worm,radius_worm,R2);
		}

		for(j=0;j<m;j++)            // begin of checking new worm for
		{   l=CrossTwoWorms(cl0,worms[j]);
			if(l==1)
				break;
		}

		if(l==0)
		{   worms[m]=cl0;
			m++;
		}
	}

	R2=R*R;
	for(i=0;i<i_max;i++)
	{   v0.x=(real)i;
		for(j=0;j<j_max;j++)
		{   v0.y=(real)j;
			for(l=0;l<l_max;l++)
			{   v0.z=(real)l;   v1=v0-vCenter;
				if((v1.x*v1.x)+(v1.y*v1.y)+(v1.z*v1.z)<=R2)
					for(m=0;m<iWormsNumber;m++)
					{   cAlpha[i][j][l]=SeedToWorm(worms[m],v1);
						if(cAlpha[i][j][l]!=0)
						{   iTmp++;
							break;
						}
					}
			}
		}
	}

	return iTmp;
}

/******************** Generation of hedgehog-like particles *****************/

int iGenerateHedgehog(int i_max, int j_max, int l_max, real radius_worm,
					  real length_worm, int iWormsNumber, cylinder * worms,
					  char *** cAlpha)
{   int i,j,l,m,iTmp=0;
	real R=(real)i_max/2.,R2=(R-radius_worm)*(R-radius_worm);
	real r1,r2,theta1,theta2,phi1,phi2;
	vector vCenter={(real)((real)i_max/2.-.5),(real)((real)j_max/2.-.5),(real)((real)l_max/2.-.5)},v0,v1;
	cylinder cl0;

	for(i=0;i<i_max;i++)
		for(j=0;j<j_max;j++)
			for(l=0;l<l_max;l++)
				cAlpha[i][j][l]=0;

	l=0;    m=0;
	while(m<iWormsNumber)
	{   cl0.r=-1.;
		while(cl0.r<0.)
		{   r1=R*(real)dda_rand()/(real)DDA_RAND_MAX;

			if((2.*(real)dda_rand()/(real)DDA_RAND_MAX-1.)<=0.)
				r2=r1-0.1;
			else
				r2=r1+0.1;

			phi1=2.*PI*(real)dda_rand()/(real)DDA_RAND_MAX;
			theta1=acos(2.*(real)dda_rand()/(real)DDA_RAND_MAX-1.);

			phi2=phi1;
			theta2=theta1;

			cl0=worm(r1,theta1,phi1,r2,theta2,phi2,length_worm,radius_worm,R2);
		}

		for(j=0;j<m;j++)            // begin of checking new worm for
		{   l=CrossTwoWorms(cl0,worms[j]);
			if(l==1)
				break;
		}

		if(l==0)
		{   worms[m]=cl0;
			m++;
		}
	}

	R2=R*R;
	for(i=0;i<i_max;i++)
	{   v0.x=(real)i;
		for(j=0;j<j_max;j++)
		{   v0.y=(real)j;
			for(l=0;l<l_max;l++)
			{   v0.z=(real)l;   v1=v0-vCenter;
				if((v1.x*v1.x)+(v1.y*v1.y)+(v1.z*v1.z)<=R2)
					for(m=0;m<iWormsNumber;m++)
					{   cAlpha[i][j][l]=SeedToWorm(worms[m],v1);
						if(cAlpha[i][j][l]!=0)
						{   iTmp++;
							break;
						}
					}
			}
		}
	}

	return iTmp;
}

/******************* Generation of dead anthrax-like particles **************/

int iGenerateDeadAnthrax(int i_max, int j_max, int l_max, real radius_worm,
						 real length_worm, int iWormsNumber, cylinder * worms,
						 char *** cAlpha)
{   int i,j,l,m,iTmp=0;
	real R=(real)i_max/2.,R2=(R-radius_worm)*(R-radius_worm);
	real r0,r1,r2,theta1,phi1;
	vector vCenter={(real)((real)i_max/2.-.5),(real)((real)j_max/2.-.5),(real)((real)l_max/2.-.5)},v0,v1;
	cylinder cl0;

	for(i=0;i<i_max;i++)
		for(j=0;j<j_max;j++)
			for(l=0;l<l_max;l++)
				cAlpha[i][j][l]=0;

	l=0;    m=0;
	while(m<iWormsNumber)
	{   cl0.r=-1.;
		while(cl0.r<0.)
		{   r1=R*(real)dda_rand()/(real)DDA_RAND_MAX;
			phi1=2.*PI*(real)dda_rand()/(real)DDA_RAND_MAX;
			r0=(2.*(real)dda_rand()/(real)DDA_RAND_MAX-1.);
			theta1=acos(r0);
			r2=sin(theta1);

			cl0.l=length_worm;
			cl0.n=vrand();

			cl0.b.x=r1*r2*cos(phi1);
			cl0.b.y=r1*r2*sin(phi1);
			cl0.b.z=r1*r0;

			cl0.e=cl0.b+cl0.l*cl0.n;

			// is end point of cylinder in the initial matrix?
			if(cl0.e*cl0.e<=R2)
				cl0.r=radius_worm;
		}

		for(j=0;j<m;j++)            // begin of checking new worm for
		{   l=CrossTwoWorms(cl0,worms[j]);
			if(l==1)
				break;
		}

		if(l==0)
		{   worms[m]=cl0;
			m++;
		}
	}

	R2=R*R;
	for(i=0;i<i_max;i++)
	{   v0.x=(real)i;
		for(j=0;j<j_max;j++)
		{   v0.y=(real)j;
			for(l=0;l<l_max;l++)
			{   v0.z=(real)l;   v1=v0-vCenter;
				if((v1.x*v1.x)+(v1.y*v1.y)+(v1.z*v1.z)<=R2)
					for(m=0;m<iWormsNumber;m++)
					{   cAlpha[i][j][l]=SeedToWorm(worms[m],v1);
						if(cAlpha[i][j][l]!=0)
						{   iTmp++;
							break;
						}
					}
			}
		}
	}

	return iTmp;
}

/***************************** Cross of two worms ***************************/

int CrossTwoWorms(cylinder cl0, cylinder cl1)
{   int i0,i1;
	real a11,a12,a21,a22,b1,b2,r0,r1,r2,d,R2;
	vector v0,v1,vNearPoint;

	v0=cl0.n^cl1.n;
	r0=abs(v0);

	if(r0<=1.e-20)      // cylinders are parallel
	{   return 1;
	}

	v1=cl1.b-cl0.b;
	d=fabs((v1*v0)/r0); // distance between two infinite cylinders

	r2=(cl0.r+cl1.r);   R2=r2*r2;
	if(d>r2)            // cylinders are not intersecting never
		return 0;
						// DEFINITION OF THE MOST NEAREST POINTS
						// coefficients of system linear algebraic equations
	b1=v1*cl0.n;        // a11*r0+a12*r1=b1;
	b2=v1*cl1.n;        // a21*r0+a22*r1=b2;

	a11=1.; a12=cl1.n*cl0.n; a21=-a12; a22=-a11;

	r2=a11*a22-a12*a21;
	r0=(b1*a22-b2*a12)/r2;  // parameter t in eq. r=r0+n*t
	r1=(b2*a11-b1*a21)/r2;  // for first and second worms

	if(r0>=0.&&r0<=cl0.l&&r1>=0.&&r1<=cl1.l)
		return 1;       // cylinders are intersecting

	if(r0>=0.&&r0<=cl0.l)
	{   vNearPoint=cl0.b+cl0.n*r0;
		i0=1;
	}
	else
		i0=0;

	if(r1>=0.&&r1<=cl1.l)
	{   vNearPoint=cl1.b+cl1.n*r1;
		i1=1;
	}
	else
		i1=0;

	if(PointToWorm(cl0,cl1.b,cl1.r)==1)
		return 1;   // cylinders are intersecting

	if(PointToWorm(cl0,cl1.e,cl1.r)==1)
		return 1;   // cylinders are intersecting

	if(PointToWorm(cl1,cl0.b,cl0.r)==1)
		return 1;   // cylinders are intersecting

	if(PointToWorm(cl1,cl0.e,cl0.r)==1)
		return 1;   // cylinders are intersecting

	if(i0==1&&PointToWorm(cl1,vNearPoint,cl0.r)==1)
		return 1;   // cylinders are intersecting

	if(i1==1&&PointToWorm(cl0,vNearPoint,cl1.r)==1)
		return 1;   // cylinders are intersecting

	return 0;   // cylinders are not intersecting
}

/**************************** Generation of worm ****************************/

cylinder worm(real r1, real theta1, real phi1, real r2, real theta2, real phi2,
			  real length, real radius, real R2)
{   real c_theta1=cos(theta1),s_theta1=sin(theta1),c_phi1=cos(phi1),s_phi1=sin(phi1);
	vector begin={r1,phi1,theta1},end={r2,phi2,theta2},direction=end-begin,v0;
	cylinder clTmp;

	// translation start point of cylinder
	// from spherical to Cartesian coordinates
	begin.x=r1*s_theta1*c_phi1;
	begin.y=r1*s_theta1*s_phi1;
	begin.z=r1*c_theta1;

	// is begin point of cylinder in the initial matrix?
	if(begin*begin>R2)
	{   clTmp.r=-1.;
		return clTmp;
	}

	// translation vector direction from spherical to Cartesian coordinates
	v0.x=s_theta1*c_phi1*direction.x-s_phi1*direction.y+c_theta1*c_phi1*direction.z;
	v0.y=s_theta1*s_phi1*direction.x+c_phi1*direction.y+c_theta1*s_phi1*direction.z;
	v0.z=c_theta1*direction.x-s_theta1*direction.z;
	direction=1./abs(v0)*v0;

	// determination end point of cylinder
	end=begin+length*direction;

	// is end point of cylinder in the initial matrix?
	if(end*end>R2)
	{   clTmp.r=-1.;
		return clTmp;
	}

	clTmp.r=radius;
	clTmp.l=length;
	clTmp.b=begin;
	clTmp.e=end;
	clTmp.n=direction;

	return clTmp;
}

/********************** Distance between seed and worm **********************/

char SeedToWorm(cylinder clWorm, vector vPoint)
{   real t=(vPoint-clWorm.b)*clWorm.n;
	vector v0;

	if(t>=0.&&t<=clWorm.l)
	{   v0=(clWorm.b-vPoint)+clWorm.n*t;
		if((v0*v0)<=(clWorm.r*clWorm.r))
			return 1;
		else
			return 0;
	}

	if(t>(-clWorm.r)&&t<0.)
	{   v0=clWorm.b-vPoint;
		if((v0*v0)<=(clWorm.r*clWorm.r))
			return 1;
		else
			return 0;
	}

	if(t>clWorm.l&&t<(clWorm.l+clWorm.r))
	{   v0=clWorm.e-vPoint;
		if((v0*v0)<=(clWorm.r*clWorm.r))
			return 1;
		else
			return 0;
	}

	return 0;
}

/********************* Distance between point and worm **********************/

int PointToWorm(cylinder clWorm, vector vPoint, real radius)
{   real t=(vPoint-clWorm.b)*clWorm.n,max_d=clWorm.r+radius;
	real d=abs((clWorm.b-vPoint)+clWorm.n*t);

	if(d>max_d)
		return 0;   // worms are not intersecting

	if(t>=0.&&t<=clWorm.l)
		return 1;   // worms are intersecting
	else
	{   if(t<0)
		{   d=abs(clWorm.b-vPoint);
			if(d>max_d)
				return 0;   // worms are not intersecting
			else
				return 1;   // worms are intersecting
		}
		else
		{   d=abs(clWorm.e-vPoint);
			if(d>max_d)
				return 0;   // worms are not intersecting
			else
				return 1;   // worms are intersecting
		}
	}
}

/************************* Generation of ellipsoid **************************/

int iGenerateEllipsoid(int i_max, int j_max, int l_max, real b_a, real c_a,
					   char *** cAlpha)
{   int i,j,l;
	vector vCenter={(real)((real)i_max/2.-.5),(real)((real)j_max/2.-.5),(real)((real)l_max/2.-.5)},v0,v1;
	real a,b,c;                         // axes of the ellipsoid
	int iTmp=0;

	a=(real)i_max;  b=a*b_a;  c=a*c_a;  // define axes of the ellipsoid

	if(b>(real)j_max||c>(real)l_max)
	{   printf("\n Ellipsoid does not be located in the space! ");
		printf("\n Please augment the size of space or press the ellipsoid! ");
		fflush(stdout);
		return -1;
	}

	a=a/2.; b=b/2.; c=c/2.; // semi-axes of the ellipsoid
	a=a*a;  b=b*b;  c=c*c;  // semi-axes^2 of the ellipsoid

	for(i=0;i<i_max;i++)
	{   v0.x=(real)i;
		for(j=0;j<j_max;j++)
		{   v0.y=(real)j;
			for(l=0;l<l_max;l++)
			{   v0.z=(real)l;   v1=v0-vCenter;
				if((v1.x*v1.x)/a+(v1.y*v1.y)/b+(v1.z*v1.z)/c<=1.)
				{   iTmp++;
					cAlpha[i][j][l]=1;
				}
				else
					cAlpha[i][j][l]=0;
			}
		}
	}

	return iTmp;
}

/**************************** Generate bisphere *****************************/

int iGenerateBisphere(int i_max, int j_max, int l_max, real R1_R2, char *** cAlpha)
{   int i,j,l;
	real R1,R2,R1_2,R2_2,d1,d2;
	vector vCenter1,vCenter2,v0,v1,v2;
	int iTmp=0;

	if(R1_R2>1.)
	{   R1=(real)j_max/2.;
		R2=R1/R1_R2;

		vCenter1.x=(real)i_max/4.-.5;
		vCenter2.x=(R1+R2)+vCenter1.x;
	}
	else
	{   R2=(real)j_max/2.;
		R1=R2*R1_R2;

		vCenter2.x=(real)i_max*3./4.-.5;
		vCenter1.x=vCenter2.x-(R1+R2);
	}

	vCenter1.y=vCenter2.y=(real)j_max/2.-.5;
	vCenter1.z=vCenter2.z=(real)l_max/2.-.5;
	R1_2=R1*R1; R2_2=R2*R2;

	for(i=0;i<i_max;i++)
	{   v0.x=(real)i;
		for(j=0;j<j_max;j++)
		{   v0.y=(real)j;
			for(l=0;l<l_max;l++)
			{   v0.z=(real)l;

				v1=v0-vCenter1; d1=(v1.x*v1.x)+(v1.y*v1.y)+(v1.z*v1.z);
				v2=v0-vCenter2; d2=(v2.x*v2.x)+(v2.y*v2.y)+(v2.z*v2.z);

				if(d1>R1_2&&d2>R2_2)
					cAlpha[i][j][l]=0;
				else
				{   iTmp++;
					if(d1<=R1_2)
						cAlpha[i][j][l]=1;
					if(d2<=R2_2)
						cAlpha[i][j][l]=2;
				}
			}
		}
	}

	return iTmp;
}

/************************* Generate parallelepiped **************************/

int iGenerateParallelepiped(int i_max, int j_max, int l_max, real b_a, real c_a,
							char *** cAlpha)
{   int i,j,l;
	real a,b,c;                         // axes of the parallelepiped
	int iTmp=0;

	a=(real)i_max;  b=a*b_a;  c=a*c_a;  // define axes of the parallelepiped

	if(b>(real)j_max||c>(real)l_max)
	{   printf("\n Parallelepiped does not be located in the space! ");
		printf("\n Please augment the size of space or press the parallelepiped! ");
		fflush(stdout);
		return -1;
	}

	for(i=0;i<i_max;i++)
		for(j=0;j<j_max;j++)
			for(l=0;l<l_max;l++)
				if((a-(real)i)>=0.001&&(b-(real)j)>=0.001&&(c-(real)l)>=0.001)
				{   iTmp++;
					cAlpha[i][j][l]=1;
				}
				else
					cAlpha[i][j][l]=0;

	return iTmp;
}

/************** Generate random particle in ellipsoidal matrix **************/

int iGenerateRandomEllipsoid(int i_max, int j_max, int l_max, real b_a, real c_a,
							 int iSorts, int * iSeeds, char *** cAlpha)
{   int i,j,l,m,iFlag,iiNoSeeds,iNoSeeds=0,iTmp=0;
	real a,b,c,dist,dist_min;
	vector vCenter={(real)((real)i_max/2.-.5),(real)((real)j_max/2.-.5),(real)((real)l_max/2.-.5)},v0,v1;
	vector * vSeeds;

	a=(real)i_max;  b=a*b_a;  c=a*c_a;  // define axes of the ellipsoid

	if(b>(real)j_max||c>(real)l_max)
	{   printf("\n Ellipsoid does not be located in the space! ");
		printf("\n Please augment the size of space or press the ellipsoid! ");
		fflush(stdout);
		return -1;
	}

	for(m=0;m<iSorts;m++)
		iNoSeeds+=iSeeds[m];

	vSeeds=(vector *)malloc(iNoSeeds*sizeof(vector));
	if(vSeeds==NULL)
	{   printf("\n No enough memory for array vSeeds! ");
		return 2;
	}

	a=a/2.; b=b/2.; c=c/2.; // semi-axes of the ellipsoid
	a=a*a;  b=b*b;  c=c*c;  // semi-axes^2 of the ellipsoid
							// Set the seeds to a random places in matrix
	for(m=0;m<iNoSeeds;m++)
	{   do
		{   iFlag=0;
			v0.x=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)i_max-1.));
			v0.y=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)j_max-1.));
			v0.z=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)l_max-1.));

			for(i=0;i<m;i++)
				if(abs(vSeeds[i]-v0)<ZERO)
				{   iFlag=1;
					break;
				}

			v1=v0-vCenter;
		}
		while(((v1.x*v1.x)/a+(v1.y*v1.y)/b+(v1.z*v1.z)/c>1.)&&iFlag==1);
		vSeeds[m]=v0;
	}

	for(i=0;i<i_max;i++)
	{   v0.x=(real)i;
		for(j=0;j<j_max;j++)
		{   v0.y=(real)j;
			for(l=0;l<l_max;l++)
			{   v0.z=(real)l;   v1=v0-vCenter;
				if((v1.x*v1.x)/a+(v1.y*v1.y)/b+(v1.z*v1.z)/c<=1.)
				{   iFlag=0;
					dist_min=sqrt((real)(i_max*i_max+j_max*j_max+l_max*l_max));
					for(m=0;m<iNoSeeds;m++)
					{   dist=abs(v0-vSeeds[m]);
						if(dist<dist_min)
						{   iFlag=m;    dist_min=dist;
						}
					}

					iiNoSeeds=iSeeds[0];
					if(iFlag<iiNoSeeds)
					{   cAlpha[i][j][l]=1;  iTmp++;
					}

					for(m=1;m<iSorts;m++)
					{   if(iFlag>=iiNoSeeds&&iFlag<(iiNoSeeds+iSeeds[m]))
						{   cAlpha[i][j][l]=(m+1);  break;
						}
						iiNoSeeds+=iSeeds[m];
					}
				}
				else
					cAlpha[i][j][l]=0;      // free space
			}
		}
	}

	free(vSeeds);
	return iTmp;
}

int iGenerateRandomParallelepiped(int i_max, int j_max, int l_max, real b_a,
								  real c_a, int iSorts, int * iSeeds,
								  char *** cAlpha)
{   int i,j,l,m,iFlag,iiNoSeeds,iNoSeeds=0,iTmp=0;
	real a,b,c,dist,dist_min;
	vector vCenter={(real)((real)i_max/2.-.5),(real)((real)j_max/2.-.5),(real)((real)l_max/2.-.5)},v0;
	vector * vSeeds;

	a=(real)i_max;  b=a*b_a;  c=a*c_a;  // define axes of the parallelepiped

	if(b>(real)j_max||c>(real)l_max)
	{   printf("\n Parallelepiped does not be located in the space! ");
		printf("\n Please augment the size of space or press the parallelepiped! ");
		fflush(stdout);
		return -1;
	}

	for(m=0;m<iSorts;m++)
		iNoSeeds+=iSeeds[m];

	vSeeds=(vector *)malloc(iNoSeeds*sizeof(vector));
	if(vSeeds==NULL)
	{   printf("\n No enough memory for array vSeeds! ");
		return 2;
	}

							// Set the seeds to a random places in matrix
	for(m=0;m<iNoSeeds;m++)
	{   do
		{   iFlag=0;
			v0.x=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)i_max-1.));
			v0.y=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)j_max-1.));
			v0.z=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)l_max-1.));

			for(i=0;i<m;i++)
				if(abs(vSeeds[i]-v0)<ZERO)
				{   iFlag=1;
					break;
				}
		}
		while(iFlag==1);
		vSeeds[m]=v0;
	}

	for(i=0;i<i_max;i++)
	{   v0.x=(real)i;
		for(j=0;j<j_max;j++)
		{   v0.y=(real)j;
			for(l=0;l<l_max;l++)
			{   v0.z=(real)l;

				iFlag=0;
				dist_min=sqrt((real)(i_max*i_max+j_max*j_max+l_max*l_max));
				for(m=0;m<iNoSeeds;m++)
				{   dist=abs(v0-vSeeds[m]);
					if(dist<dist_min)
					{   iFlag=m;    dist_min=dist;
					}
				}

				iiNoSeeds=iSeeds[0];
				if(iFlag<iiNoSeeds)
				{   cAlpha[i][j][l]=1;  iTmp++;
				}

				for(m=1;m<iSorts;m++)
				{   if(iFlag>=iiNoSeeds&&iFlag<(iiNoSeeds+iSeeds[m]))
					{   cAlpha[i][j][l]=(m+1);  break;
					}
					iiNoSeeds+=iSeeds[m];
				}

			}
		}
	}

	free(vSeeds);
	return iTmp;
}

int iGenerateNibbleEllipsoid(int i_max, int j_max, int l_max, real b_a, real c_a,
							 int iSorts, int * iSeeds, int iBites, real rDepth,
							 int * iListMV, char *** cAlpha)
{   int i,j,l,m,iFlag,iiNoSeeds,iNoSeeds=0,iTmp=0;
	real a,b,c,a1,b1,c1,dist,dist_min,v1x2,v1y2,v1z2;
	vector vCenter={(real)((real)i_max/2.-.5),(real)((real)j_max/2.-.5),(real)((real)l_max/2.-.5)},v0,v1;
	vector * vSeeds;

	a=(real)i_max;  b=a*b_a;  c=a*c_a;  // define axes of the ellipsoid
										// define axes of the small ellipsoid
	a1=(1.-rDepth)*a;   b1=(1.-rDepth)*b;   c1=(1.-rDepth)*c;

	if(b>(real)j_max||c>(real)l_max)
	{   printf("\n Ellipsoid does not be located in the space! ");
		printf("\n Please augment the size of space or press the ellipsoid! ");
		fflush(stdout);
		return -1;
	}

	for(m=0;m<iSorts;m++)
		iNoSeeds+=iSeeds[m];

	vSeeds=(vector *)malloc((iNoSeeds+iBites)*sizeof(vector));
	if(vSeeds==NULL)
	{   printf("\n No enough memory for array vSeeds! ");
		return 2;
	}

	a=a/2.; b=b/2.; c=c/2.;         // semi-axes of the ellipsoid
	a=a*a;  b=b*b;  c=c*c;          // semi-axes^2 of the ellipsoid

	a1=a1/2.; b1=b1/2.; c1=c1/2.;   // semi-axes of the small ellipsoid
	a1=a1*a1; b1=b1*b1; c1=c1*c1;   // semi-axes^2 of the small ellipsoid

							// Set the seeds to a random places in matrix
	for(m=0;m<iNoSeeds;m++)
	{   do
		{   iFlag=0;
			v0.x=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)i_max-1.));
			v0.y=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)j_max-1.));
			v0.z=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)l_max-1.));

			for(i=0;i<m;i++)
				if(abs(vSeeds[i]-v0)<ZERO)
				{   iFlag=1;
					break;
				}

			v1=v0-vCenter;
		}
		while(((v1.x*v1.x)/a1+(v1.y*v1.y)/b1+(v1.z*v1.z)/c1>1.)&&iFlag==1);
		vSeeds[m]=v0;
	}

	for(m=iNoSeeds;m<(iNoSeeds+iBites);m++)
	{   do
		{   iFlag=0;
			v0.x=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)i_max-1.));
			v0.y=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)j_max-1.));
			v0.z=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)l_max-1.));

			v1=v0-vCenter; v1x2=v1.x*v1.x; v1y2=v1.y*v1.y; v1z2=v1.z*v1.z;

			if((v1x2/a+v1y2/b+v1z2/c)<=1.&&(v1x2/a1+v1y2/b1+v1z2/c1)>1.)
			{   for(i=iNoSeeds;i<m;i++)
					if(abs(vSeeds[i]-v0)<ZERO)
					{   iFlag=1;
						break;
					}
			}
			else
				iFlag=1;
		}
		while(iFlag==1);

		vSeeds[m]=v0;
	}

	for(i=0;i<i_max;i++)
	{   v0.x=(real)i;
		for(j=0;j<j_max;j++)
		{   v0.y=(real)j;
			for(l=0;l<l_max;l++)
			{   v0.z=(real)l;   v1=v0-vCenter;
				if((v1.x*v1.x)/a+(v1.y*v1.y)/b+(v1.z*v1.z)/c<=1.)
				{   iFlag=0;
					dist_min=sqrt((real)(i_max*i_max+j_max*j_max+l_max*l_max));
					for(m=0;m<(iNoSeeds+iBites);m++)
					{   dist=abs(v0-vSeeds[m]);
						if(dist<dist_min)
						{   iFlag=m;    dist_min=dist;
						}
					}

					if(iFlag>=iNoSeeds&&iFlag<(iNoSeeds+iBites))
						cAlpha[i][j][l]=0;
					else
					{   iiNoSeeds=iSeeds[0];
						if(iFlag<iiNoSeeds)
						{   cAlpha[i][j][l]=1;
						}
						else
						{   for(m=1;m<iSorts;m++)
							{   if(iFlag>=iiNoSeeds&&iFlag<(iiNoSeeds+iSeeds[m]))
								{   cAlpha[i][j][l]=(m+1);  break;
								}
								iiNoSeeds+=iSeeds[m];
							}
						}
					}
				}
				else
					cAlpha[i][j][l]=0;      // free space

				m=(int)cAlpha[i][j][l];
				if(iListMV[m]!=0)
				{   iTmp++;
				}
			}
		}
	}

	free(vSeeds);
	return iTmp;
}

/******************* Generate random hierarchical particles *****************/

int iGenerateHierarchicalEllipsoid(int iHierarchicalLevels, int i_max, int j_max, int l_max,
								   real b_a, real c_a, int iSorts, int * iSeeds,
								   char *** cAlpha)
{   int i,j,l,m,iFlag,iiNoSeeds,iNoSeeds=0,iTmp=0,iLevelsCounter=1,i0;
	real a,b,c,dist,dist_min;
	vector vCenter={(real)((real)i_max/2.-.5),(real)((real)j_max/2.-.5),(real)((real)l_max/2.-.5)},v0,v1;
	vector * vSeeds;

	a=(real)i_max;  b=a*b_a;  c=a*c_a;  // define axes of the ellipsoid

	if(b>(real)j_max||c>(real)l_max)
	{   printf("\n Ellipsoid does not be located in the space! ");
		printf("\n Please augment the size of space or press the ellipsoid! ");
		fflush(stdout);
		return -1;
	}

	for(m=0;m<iSorts;m++)
		iNoSeeds+=iSeeds[m];

	vSeeds=(vector *)malloc(iNoSeeds*sizeof(vector));
	if(vSeeds==NULL)
	{   printf("\n No enough memory for array vSeeds! ");
		return 2;
	}

	a=a/2.; b=b/2.; c=c/2.; // semi-axes of the ellipsoid
	a=a*a;  b=b*b;  c=c*c;  // semi-axes^2 of the ellipsoid
							// Set the seeds to a random places in matrix
	for(m=0;m<iNoSeeds;m++)
	{   do
		{   iFlag=0;
			v0.x=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)i_max-1.));
			v0.y=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)j_max-1.));
			v0.z=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)l_max-1.));

			for(i=0;i<m;i++)
				if(abs(vSeeds[i]-v0)<ZERO)
				{   iFlag=1;
					break;
				}

			v1=v0-vCenter;
		}
		while(((v1.x*v1.x)/a+(v1.y*v1.y)/b+(v1.z*v1.z)/c>1.)&&iFlag==1);
		vSeeds[m]=v0;
	}

	for(i=0;i<i_max;i++)
	{   v0.x=(real)i;
		for(j=0;j<j_max;j++)
		{   v0.y=(real)j;
			for(l=0;l<l_max;l++)
			{   v0.z=(real)l;   v1=v0-vCenter;
				if((v1.x*v1.x)/a+(v1.y*v1.y)/b+(v1.z*v1.z)/c<=1.)
				{   iFlag=0;
					dist_min=sqrt((real)(i_max*i_max+j_max*j_max+l_max*l_max));
					for(m=0;m<iNoSeeds;m++)
					{   dist=abs(v0-vSeeds[m]);
						if(dist<dist_min)
						{   iFlag=m;    dist_min=dist;
						}
					}

					iiNoSeeds=iSeeds[0];
					if(iFlag<iiNoSeeds)
					{   cAlpha[i][j][l]=1;  iTmp++;
					}

					for(m=1;m<iSorts;m++)
					{   if(iFlag>=iiNoSeeds&&iFlag<(iiNoSeeds+iSeeds[m]))
						{   cAlpha[i][j][l]=(m+1);  break;
						}
						iiNoSeeds+=iSeeds[m];
					}
				}
				else
					cAlpha[i][j][l]=0;      // free space
			}
		}
	}

	if(iLevelsCounter<iHierarchicalLevels)
	{   iTmp=0;
		for(m=0;m<iNoSeeds;m++)
		{   do
			{   iFlag=0;
				do
				{   v0.x=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)i_max-1.));
					v0.y=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)j_max-1.));
					v0.z=rounding(((real)dda_rand()/(real)DDA_RAND_MAX)*((real)l_max-1.));
					i0=cAlpha[(int)v0.x][(int)v0.y][(int)v0.z];
				}
				while(i0!=1);

				for(i=0;i<m;i++)
					if(abs(vSeeds[i]-v0)<ZERO)
					{   iFlag=1;
						break;
					}
			}
			while(iFlag==1);
			vSeeds[m]=v0;
		}

		for(i=0;i<i_max;i++)
			for(j=0;j<j_max;j++)
				for(l=0;l<l_max;l++)
					if((int)cAlpha[i][j][l]==1)
					{   v0.x=(real)i;   v0.y=(real)j;   v0.z=(real)l;
						iFlag=0;
						dist_min=sqrt((real)(i_max*i_max+j_max*j_max+l_max*l_max));
						for(m=0;m<iNoSeeds;m++)
						{   dist=abs(v0-vSeeds[m]);
							if(dist<dist_min)
							{   iFlag=m;    dist_min=dist;
							}
						}

						iiNoSeeds=iSeeds[0];
						if(iFlag<iiNoSeeds)
						{   cAlpha[i][j][l]=1;  iTmp++;
						}

						for(m=1;m<iSorts;m++)
						{   if(iFlag>=iiNoSeeds&&iFlag<(iiNoSeeds+iSeeds[m]))
							{   cAlpha[i][j][l]=(m+1);  break;
							}
							iiNoSeeds+=iSeeds[m];
						}
					}

		iLevelsCounter++;
	}

	free(vSeeds);
	return iTmp;
}

/****** Generate sphere which do not fill all initial spherical volume ******/

int iGenerateSphereSpec(int i_max, int j_max, int l_max, real Radius, real rho,
					   char *** cAlpha)
{   int i,j,l;
	vector vCenter={(real)((real)i_max/2.-.5),(real)((real)j_max/2.-.5),(real)((real)l_max/2.-.5)},v0,v1;
	real radius,radius2;
	int iTmp=0;

	radius=pow((real)rho,(real)(1./3.))*Radius; radius2=radius*radius;

	for(i=0;i<i_max;i++)
	{   v0.x=(real)i;
		for(j=0;j<j_max;j++)
		{   v0.y=(real)j;
			for(l=0;l<l_max;l++)
			{   v0.z=(real)l;   v1=v0-vCenter;
				if((v1.x*v1.x+v1.y*v1.y+v1.z*v1.z)<=radius2)
				{   iTmp++;
					cAlpha[i][j][l]=1;
				}
				else
					cAlpha[i][j][l]=0;
			}
		}
	}

	return iTmp;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void Particle_Boundary(int iNoD, int i_max_ini, int j_max_ini, int l_max_ini,
					   char *** cAlpha, cplx * cListMaterials, int * X_min,
					   int * Y_min, int * Z_min, int * X_max, int * Y_max,
					   int * Z_max)
{   int i,j,l,m,n=0;
	int _X_min=1048576,_X_max=0,_Y_min=1048576,_Y_max=0,_Z_min=1048576,_Z_max=0;

	for(i=0;i<i_max_ini;i++)
		for(j=0;j<j_max_ini;j++)
			for(l=0;l<l_max_ini;l++)
			{   m=(int)cAlpha[i][j][l];

				if(cListMaterials[m].r>ZERO||cListMaterials[m].i>ZERO)
				{   n++;

					if(i<_X_min)
						_X_min=i;

					if(j<_Y_min)
						_Y_min=j;

					if(l<_Z_min)
						_Z_min=l;


					if(i>_X_max)
						_X_max=i;

					if(j>_Y_max)
						_Y_max=j;

					if(l>_Z_max)
						_Z_max=l;
				}
			}


	if(n!=iNoD)
	{   printf("\n The number of cells tested is not equal to total number of dipoles");
		fflush(stdout);
	}

	(* X_min)=_X_min;
	(* Y_min)=_Y_min;
	(* Z_min)=_Z_min;

	(* X_max)=_X_max;
	(* Y_max)=_Y_max;
	(* Z_max)=_Z_max;

	return;
}
