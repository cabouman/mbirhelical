/* version 6.0, P.Jin 07/13/2012 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sched.h>
#include <omp.h>
#include <sys/time.h>
#include <time.h>	/* for ctime */
#include "allocate.h"
#include "data.h"
#include "proj.h"
#include "prepro.h"
#include "solve.h"
#include "io.h"		/* for write image */
#include "quicksort.h" /* sjk: contains quickselect()  */
#include "icd.h"

void fillNeighbors(float *neighbors, int jx, int jy, int jz, ENTRY *X, struct ImgInfo *img_info)
{
	int plusx, minusx, plusy, minusy, plusz, minusz;
	int Nz,Nyz;

	/* use circular boundary condition */
	plusx = jx + 1;
	plusx = ((plusx < img_info->Nx) ? plusx : 0);
	minusx = jx - 1;
	minusx = ((minusx < 0) ? (img_info->Nx-1) : minusx);
	plusy = jy + 1;
	plusy = ((plusy < img_info->Ny) ? plusy : 0);
	minusy = jy - 1;
	minusy = ((minusy < 0) ? (img_info->Ny-1) : minusy);
	plusz = jz + 1;
	plusz = ((plusz < img_info->Nz) ? plusz : 0);
	minusz = jz - 1;
	minusz = ((minusz < 0) ? (img_info->Nz-1) : minusz);

	Nz = img_info->Nz;
	Nyz = img_info->Ny*img_info->Nz;

	/* jx-1 plane */
	neighbors[0] = X[minusx*Nyz+minusy*Nz+minusz];
	neighbors[1] = X[minusx*Nyz+minusy*Nz+jz];
	neighbors[2] = X[minusx*Nyz+minusy*Nz+plusz];

	neighbors[3] = X[minusx*Nyz+jy*Nz+minusz];
	neighbors[4] = X[minusx*Nyz+jy*Nz+jz];
	neighbors[5] = X[minusx*Nyz+jy*Nz+plusz];

	neighbors[6] = X[minusx*Nyz+plusy*Nz+minusz];
	neighbors[7] = X[minusx*Nyz+plusy*Nz+jz];
	neighbors[8] = X[minusx*Nyz+plusy*Nz+plusz];

	/* jx plane */
	neighbors[9] = X[jx*Nyz+minusy*Nz+minusz];
	neighbors[10]= X[jx*Nyz+minusy*Nz+jz];
	neighbors[11]= X[jx*Nyz+minusy*Nz+plusz];

	neighbors[12] = X[jx*Nyz+jy*Nz+minusz];
	neighbors[13] = X[jx*Nyz+jy*Nz+plusz];

	neighbors[14] = X[jx*Nyz+plusy*Nz+minusz];
	neighbors[15] = X[jx*Nyz+plusy*Nz+jz];
	neighbors[16] = X[jx*Nyz+plusy*Nz+plusz];

	/* jx+1 plane */
	neighbors[17] = X[plusx*Nyz+minusy*Nz+minusz];
	neighbors[18] = X[plusx*Nyz+minusy*Nz+jz];
	neighbors[19] = X[plusx*Nyz+minusy*Nz+plusz];

	neighbors[20] = X[plusx*Nyz+jy*Nz+minusz];
	neighbors[21] = X[plusx*Nyz+jy*Nz+jz];
	neighbors[22] = X[plusx*Nyz+jy*Nz+plusz];

	neighbors[23] = X[plusx*Nyz+plusy*Nz+minusz];
	neighbors[24] = X[plusx*Nyz+plusy*Nz+jz];
	neighbors[25] = X[plusx*Nyz+plusy*Nz+plusz];

  if(0)  /* old convention */
  {
	/* 6 closest neighbors */
	neighbors[0] = X[plusx*Nyz+jy*img_info->Nz+jz];
	neighbors[1] = X[minusx*Nyz+jy*img_info->Nz+jz];
	neighbors[2] = X[jx*Nyz+plusy*img_info->Nz+jz];
	neighbors[3] = X[jx*Nyz+minusy*img_info->Nz+jz];
	neighbors[4] = X[jx*Nyz+jy*img_info->Nz+plusz];
	neighbors[5] = X[jx*Nyz+jy*img_info->Nz+minusz];

	/* 12 closer neighbors */
	neighbors[6] = X[plusx*Nyz+plusy*img_info->Nz+jz];
	neighbors[7] = X[minusx*Nyz+plusy*img_info->Nz+jz];
	neighbors[8] = X[plusx*Nyz+minusy*img_info->Nz+jz];
	neighbors[9] = X[minusx*Nyz+minusy*img_info->Nz+jz];
	neighbors[10] = X[plusx*Nyz+jy*img_info->Nz+plusz];
	neighbors[11] = X[minusx*Nyz+jy*img_info->Nz+plusz];
	neighbors[12] = X[plusx*Nyz+jy*img_info->Nz+minusz];
	neighbors[13] = X[minusx*Nyz+jy*img_info->Nz+minusz];
	neighbors[14] = X[jx*Nyz+plusy*img_info->Nz+plusz];
	neighbors[15] = X[jx*Nyz+minusy*img_info->Nz+plusz];
	neighbors[16] = X[jx*Nyz+plusy*img_info->Nz+minusz];
	neighbors[17] = X[jx*Nyz+minusy*img_info->Nz+minusz];

	/* 8 close neighbors */
	neighbors[18] = X[plusx*Nyz+plusy*img_info->Nz+plusz];
	neighbors[19] = X[minusx*Nyz+plusy*img_info->Nz+plusz];
	neighbors[20] = X[plusx*Nyz+minusy*img_info->Nz+plusz];
	neighbors[21] = X[minusx*Nyz+minusy*img_info->Nz+plusz];
	neighbors[22] = X[plusx*Nyz+plusy*img_info->Nz+minusz];
	neighbors[23] = X[minusx*Nyz+plusy*img_info->Nz+minusz];
	neighbors[24] = X[plusx*Nyz+minusy*img_info->Nz+minusz];
	neighbors[25] = X[minusx*Nyz+minusy*img_info->Nz+minusz];
  }
}

void compNeighborWeight(float *nb_wt, struct ImgInfo *img_info)
{
	int j,jx,jy,jz;
	float sum,Dxy2,Dz2;

	Dxy2=img_info->Del_xy * img_info->Del_xy;
	Dz2= img_info->Del_z * img_info->Del_z;

	j=0;
	sum=0;
	for(jx=-1; jx<=1; jx++)   /* NOTE: order is critical here--has to match fillNeighbors convention */
	for(jy=-1; jy<=1; jy++)
	for(jz=-1; jz<=1; jz++)
	if((jx!=0)||(jy!=0)||(jz!=0)) 
	{
		nb_wt[j] = 1/sqrt((jx*jx*Dxy2)+(jy*jy*Dxy2)+(jz*jz*Dz2));
		sum += nb_wt[j];
		j++;
	}
	for(j=0; j<26; j++)
		nb_wt[j] /= sum;

/* sjk: boosting of z-regularization */
/*
nb_wt[12] *= 2;
nb_wt[13] *= 2;
*/

}


void compNeighborWeightNonhomogeneous(
	float *nb_wt, 
	struct ImgInfo *img_info, 
	struct GeomInfo *geom_info, 
	struct SourceLocInfo *source_loc_info,
	float x, 
	float y, 
	float z)
{
	int j,jx,jy,jz,iv;
	float sum,Dxy2,Dz2,xs,ys;
	float phi11,phi12,phi21,phi22,v1,v2,norm;
	float scale=1/0.5;  /* denominator is scaling in source direction */

	/* find source location at current z */
	/* source_loc_info->zs[iv] = iv*(geom_info->Del_zs); */
	iv = (int) z/geom_info->Del_zs + 0.5;
	if(iv < 0) iv=0;
	if(iv >= geom_info->Nv) iv=geom_info->Nv-1;
	xs=source_loc_info->xs[iv];
	ys=source_loc_info->ys[iv];

	/* want to shrink the weighting in the direction of pixel-to-source */
        /* compute vectors for minor/major axes (phi1/phi2) */
	phi11=xs-x;
	phi12=ys-y;
	norm=sqrt(phi11*phi11+phi12*phi12);
	phi11/=norm;
	phi12/=norm;
	phi21=phi12;  /* phi2 is just phi1 rotated 90 deg, [0 1; -1 0]*phi1 */
	phi22=-phi11;

	Dxy2=img_info->Del_xy * img_info->Del_xy;
	Dz2= img_info->Del_z * img_info->Del_z;

	j=0;
	sum=0;
	for(jx=-1; jx<=1; jx++)   /* NOTE: order is critical here--has to match fillNeighbors convention */
	for(jy=-1; jy<=1; jy++)
	for(jz=-1; jz<=1; jz++)
	if((jx!=0)||(jy!=0)||(jz!=0)) 
	{
		v1=scale*(jx*phi11 + jy*phi12);  /* project local deviation onto phi1 */
		v2=       jx*phi21 + jy*phi22;
		nb_wt[j] = 1/sqrt(Dxy2*(v1*v1+v2*v2)+(jz*jz*Dz2));
		sum += nb_wt[j];
		j++;
	}
	for(j=0; j<26; j++)
		nb_wt[j] /= sum;

}



void updateAX(ENTRY *AX, struct ACol *col_xyz, float diff)
{
	int r;

	for (r = 0; r < col_xyz->n_index; r++)
	{
		AX[col_xyz->index[r]] += (diff*col_xyz->val[r]);
	}
}

float negLogPosterior(
	ENTRY *X,
	ENTRY *Y,
	ENTRY *AX,
	ENTRY *D,  /* sjk */
	struct PriorInfo *prior_info,
	struct GeomInfo *geom_info,
	struct ImgInfo *img_info)
{
	int i;
	int j, jx, jy, jz;
	int plusx, minusx, plusy, minusy, plusz;
	int Nyz;
	int Nvcr;
	float d;
	float delta;
	float nloglike;
	/* sjk */
	float nlogprior;
	float neighbors[26];
	float nb_wt[26];
	int nb_list[13] = {2,5,8,11,13,15,16,18,19,21,22,24,25};  /* only these used for unique cliques */

	/* negative log likelihood term */
	nloglike = 0.0;
	Nvcr=(geom_info->Nv)*(geom_info->Nc)*(geom_info->Nr);  
	for (i = 0; i < Nvcr; i++)
	{
		if (geom_info->lambda0 == 0.0)
			nloglike += (Y[i]-AX[i])*(Y[i]-AX[i]);
		else
		{
			d = D[i];   /* sjk */
			nloglike += (Y[i]-AX[i])*d*(Y[i]-AX[i]);
		}
	}
	nloglike /= 2.0;

	/* for ML estimate, return here */
	if (prior_info->est == ML)
	{
		return nloglike;
	}

	/* unified q-GMMRF, including quad, p-GGMRF */
	/* negative log prior term */

        compNeighborWeight(nb_wt,img_info);  /* sjk */
	nlogprior = 0.0;

	Nyz = img_info->Ny*img_info->Nz;
	for (jx = 0; jx < img_info->Nx; jx++)
	{
		for (jy = 0; jy < img_info->Ny; jy++)
		{
			for (jz = 0; jz < img_info->Nz; jz++)
			{
				plusx = jx + 1;
				plusx = ((plusx < img_info->Nx) ? plusx : 0);
				minusx = jx - 1;
				minusx = ((minusx < 0) ? (img_info->Nx-1) : minusx);
				plusy = jy + 1;
				plusy = ((plusy < img_info->Ny) ? plusy : 0);
				minusy = jy - 1;
				minusy = ((minusy < 0) ? (img_info->Ny-1) : minusy);
				plusz = jz + 1;
				plusz = ((plusz < img_info->Nz) ? plusz : 0);

				j = jx*Nyz + jy*img_info->Nz + jz;

				fillNeighbors(neighbors,jx,jy,jz,X,img_info); 
				for(i=0; i<13; i++)
				{
					delta = X[j]-neighbors[nb_list[i]]; 
					nlogprior += nb_wt[nb_list[i]]*qpotential(delta,prior_info->q,prior_info->p,prior_info->c);
				}
			}
		}
	}
	nlogprior /= (2.0*pow(prior_info->sigma,prior_info->q));

	printf("nloglike = %.1f\tnlogprior = %.1f\n",nloglike,nlogprior);

	return (nloglike + nlogprior);
}

float ICDStep(
	struct ICDInfo *icd_info,
	struct PriorInfo *prior_info,
	int jx,
	int jy,
	int jz,
	ENTRY *X,			/* image, attenuation coefficients */
	ENTRY *Y,			/* line integral */
	ENTRY *AX,			/* current AX */
	unsigned short *AX_mask,	
	ENTRY *D,			/* sjk: noise matrix */
	float lambda0,
	struct ImgInfo *img_info,	/* image parameters */
	struct ACol *col_xyz)		/* A column for this (jx, jy, jz) pixel */
{
	int i, r;
	float d;
	/*float dsum=0.0;*/

	/* store old pixel value */
	icd_info->v = (float)X[jx*(img_info->Ny*img_info->Nz)+jy*img_info->Nz+jz];

	/* copmute theta1 and theta2 */
	icd_info->th1 = 0.0;
	icd_info->th2 = 0.0;
	if (lambda0 == 0.0)
	{
		for (r = 0; r < col_xyz->n_index; r++)
		{
			icd_info->th1 += (col_xyz->val[r]*(AX[col_xyz->index[r]]-Y[col_xyz->index[r]]));
			icd_info->th2 += (col_xyz->val[r]*col_xyz->val[r]);
		}
	}
	else
	{
		for (r = 0; r < col_xyz->n_index; r++)
		{
			d = D[col_xyz->index[r]];  
			/* square weighting if suspect metal is involved */
/* SJK
			if( AX_mask[col_xyz->index[r]] > 0)
				d *= d;  
			if( AX_mask[col_xyz->index[r]] > 10)
				d *= d;  
*/
			icd_info->th1 += (col_xyz->val[r]*d*(AX[col_xyz->index[r]]-Y[col_xyz->index[r]]));
			icd_info->th2 += (col_xyz->val[r]*d*col_xyz->val[r]);
/*
			dsum += d;
*/
		}
/* experimental-- renormalizing D entries locally */
/*
		icd_info->th1 *= ((float)col_xyz->n_index)/dsum;
		icd_info->th2 *= ((float)col_xyz->n_index)/dsum;
*/
	}

	/* fill the neighbors */
	/* sjk: I moved this operation before the function call to allow zero-skip check */
	/*fillNeighbors(icd_info->neighbors, jx, jy, jz, X, img_info); */

	/* for ML estimate, return here */
	if (prior_info->est == ML)
	{
		icd_info->ml = icd_info->v - icd_info->th1/icd_info->th2;
		return icd_info->ml;
	}

	/* for quadratic prior, return here */
	if (prior_info->est == QUAD)
	{
		return quadUpdate(icd_info);
	}

	/* for q-GGMRF, use substitute function */
	if (prior_info->est == QGGMRF)
	{
		return subFunctionUpdate(icd_info);
	}

	/* ELSE, p-GGMRF, use half interval search */

	/* sjk: moved this block down so it's only computed when needed */
	/* compute low and high */
	icd_info->ml = icd_info->v - icd_info->th1/icd_info->th2;
	icd_info->low = icd_info->ml;
	icd_info->high = icd_info->ml;
	for (i = 0; i < 26; i++)	/* num of neighbors = 26 */
	{
		if (icd_info->neighbors[i] < icd_info->low)
		{
			icd_info->low = icd_info->neighbors[i];
		}
		if (icd_info->neighbors[i] > icd_info->high)
		{
			icd_info->high = icd_info->neighbors[i];
		}
	}
	if (icd_info->low < 0.0)
	{
		icd_info->low = 0.0;
	}
	return halfIntervalUpdate(icd_info);

}

void ICDReconstruct(
	struct Image *image,		/* reconstructed image */
	struct Sinogram *sinogram, 	/* sinogram data */
	char **recon_mask,		/* XY reconstruction mask */
	struct PriorInfo *prior_info,	/* prior parameters */
	char *fname,			/* output file name */
	int Nit)			/* number of iterations */
{
	int it;				/* iteration index */
	float nlogpost;		/* negative log posterior */
	char name[100];
	char suffix[10];

	ENTRY *X;
	ENTRY *Y;
	ENTRY *AX;
	unsigned short *AX_mask;  /* indicator of projection hitting metal */
	ENTRY *D;  
	ENTRY **filter;
	int *order;

	int t;
	struct SourceLocInfo source_loc_info;


	FILE *fp;
	time_t now;

/* TEST mode */
/*
int j,jx,jy,jz;
struct ICDInfo icd_info;
createSourceLocInfo(&source_loc_info, &sinogram->geom_info);
compSourceLocInfo(&source_loc_info, &sinogram->geom_info);
printf("xs=%.2f, ys=%.2f\n\n",source_loc_info.xs[0],source_loc_info.ys[0]);
compNeighborWeight(icd_info.nb_wt,&image->img_info);
compNeighborWeightNonhomogeneous(icd_info.nb_wt,&image->img_info,&sinogram->geom_info,&source_loc_info,0,0,0);

j=0;
for(jx=-1; jx<=1; jx++)   
for(jy=-1; jy<=1; jy++)
for(jz=-1; jz<=1; jz++)
if((jx!=0)||(jy!=0)||(jz!=0)) 
{
  if((jz==-1)) 
    printf("%.3f ",icd_info.nb_wt[j]);
  j++;
}
printf("\n\n");
j=0;
for(jx=-1; jx<=1; jx++)   
for(jy=-1; jy<=1; jy++)
for(jz=-1; jz<=1; jz++)
if((jx!=0)||(jy!=0)||(jz!=0)) 
{
  if((jz==0)) 
    printf("%.3f ",icd_info.nb_wt[j]);
  j++;
}
printf("\n\n");
j=0;
for(jx=-1; jx<=1; jx++)   
for(jy=-1; jy<=1; jy++)
for(jz=-1; jz<=1; jz++)
if((jx!=0)||(jy!=0)||(jz!=0)) 
{
  if((jz==1)) 
    printf("%.3f ",icd_info.nb_wt[j]);
  j++;
}
printf("\n\n");
exit(0);

*/



	/* set X, Y, allocate memory for AX */
	X = image->img;
	Y = sinogram->sino;
	D = sinogram->D; 
	 
	AX = (ENTRY *)  get_spc((sinogram->geom_info.Nv)*(sinogram->geom_info.Nc)*(sinogram->geom_info.Nr), sizeof(ENTRY));
	
	AX_mask=(unsigned short *)get_spc((sinogram->geom_info.Nv)*(sinogram->geom_info.Nc)*(sinogram->geom_info.Nr), sizeof(unsigned short));
	


	/* generate random order index */
	
	order = (int *)get_spc(image->img_info.Nx*image->img_info.Ny, sizeof(int));
	shuffle(order, image->img_info.Nx*image->img_info.Ny);
	


	/* 5x5 Hamming filter */
	
	filter = (ENTRY **)get_img(5, 5, sizeof(ENTRY));
	filter[0][0] = 0.0013;
	filter[0][4] = 0.0013;
	filter[4][0] = 0.0013;
	filter[4][4] = 0.0013;
	filter[0][1] = 0.0086;
	filter[0][3] = 0.0086;
	filter[1][0] = 0.0086;
	filter[1][4] = 0.0086;
	filter[3][0] = 0.0086;
	filter[3][4] = 0.0086;
	filter[4][1] = 0.0086;
	filter[4][3] = 0.0086;
	filter[0][2] = 0.0159;
	filter[2][0] = 0.0159;
	filter[2][4] = 0.0159;
	filter[4][2] = 0.0159;
	filter[1][1] = 0.0581;
	filter[1][3] = 0.0581;
	filter[3][1] = 0.0581;
	filter[3][3] = 0.0581;
	filter[1][2] = 0.1076;
	filter[2][1] = 0.1076;
	filter[2][3] = 0.1076;
	filter[3][2] = 0.1076;
	filter[2][2] = 0.1993;
	


	/* clip X */
	
	clipImage(X, recon_mask, &(image->img_info)); 
	/* this is critical--projector skips pixels outside mask */
	


	/* forward project X to get AX */
	if(image->img_info.Nz <= omp_get_num_threads())
		omp_set_num_threads(image->img_info.Nz);
	if(image->img_info.Nz > omp_get_num_threads() && (image->img_info.Nz)%omp_get_num_threads() !=0)
	{
       		fprintf(stderr, "the number of slices must be a multiple of the thread numbers.\n");
                exit(-1);
	}	


	
	forwardProject(AX, X, AX_mask, recon_mask, &(sinogram->geom_info), &(image->img_info));
	printf("done projecting\n");
	



	/* set estimation type */
	
	if (prior_info->q == -1.0)
	{
		prior_info->est = ML;
		fprintf(stdout, "\nThis is ML estimation!\n");
	}
	else if (prior_info->q == 2.0 && prior_info->p == 2.0)
	{
		prior_info->est = QUAD;
		fprintf(stdout, "\nThis is MAP estimation with quadratic prior!\n");
	}
	else if (prior_info->q == 2.0 && prior_info->p < 2.0 && prior_info->p >= 1.0)
	{
		prior_info->est = QGGMRF;
		fprintf(stdout, "\nThis is MAP estimation using q-GGMRF!\n");
	}
	else if (prior_info->q >= 1.0 && prior_info->q < 2.0 && prior_info->q == prior_info->p)
	{
		prior_info->est = PGGMRF;
		fprintf(stdout, "\nThis is MAP estimation using p-GGMRF!\n");
	}
	else
	{
		fprintf(stdout, "\nSoftware does not support the prior parameters provided!\n");
		exit(1);
	}

	if (sinogram->geom_info.lambda0 == 0.0)
	{
		fprintf(stdout, "\nThis simulation uses D = I weighting matrix!\n");
	}
	else
	{
		if (sinogram->geom_info.evar == 0.0)
		{
			fprintf(stdout, "\nThis simulation uses D weighting matrix without electronic noise!\n");
		}
		else
		{
			fprintf(stdout, "\nThis simulation uses D weighting matrix with electronic noise!\n");
		}
	}

	createSourceLocInfo(&source_loc_info, &sinogram->geom_info);
	compSourceLocInfo(&source_loc_info, &sinogram->geom_info);


	/* ICD iteration starts here */
	fprintf(stdout, "\nstart ICD iterations...\n");
	fprintf(stdout, "it = 0\n");
	nlogpost = negLogPosterior(X,Y,AX,D,prior_info,&(sinogram->geom_info),&(image->img_info));
	fprintf(stdout, "negative log posterior = %.1f\n", nlogpost);

	for (it = 0; it < Nit; it++)
	{
		fprintf(stdout, "it = %d\n", it+1);
		#pragma omp parallel
		{
			paraICDReconstruct(&sinogram->geom_info,&image->img_info,&source_loc_info,prior_info,recon_mask,X,Y,AX,AX_mask,D,filter,order);
		}



		/* compute negative log posterior */
	
		nlogpost = negLogPosterior(X,Y,AX,D,prior_info,&(sinogram->geom_info),&(image->img_info));
		fprintf(stdout, "negative log posterior = %.1f\n", nlogpost);
		strcpy(name, "");
		strcpy(suffix, ".txt");
		strcat(name, fname);
		strcat(name, suffix);
		fp = fopen(name, "a");
		fprintf(fp, "it = %d, cost = %f\n", it+1, nlogpost);
		time(&now);
		fprintf(fp, "time = %s\n", ctime(&now));
		fclose(fp);

		/* write out the image after each iteration */	
		strcpy(name, "");
		sprintf(suffix, "_%d.vjk", it+1);
		strcat(name, fname);
		strcat(name, suffix);
		writeImage_short(name, image);
	}

	fprintf(stdout, "\nfinish ICD iterations!\n");

	strcpy(name, "");
	strcpy(suffix, "_AX.sino");
	strcat(name, fname);
	strcat(name, suffix);
	writeSinogram_float(name,AX,sinogram->geom_info.Nr,sinogram->geom_info.Nc,sinogram->geom_info.Nv);

	freeSourceLocInfo(&source_loc_info);
	free(AX);
	free(AX_mask);
        free(order);
	free_img((void **)filter);	
	
}

void *paraICDReconstruct(struct GeomInfo *geom_info,struct ImgInfo *img_info,struct SourceLocInfo *source_loc_info,struct PriorInfo *prior_info,char **recon_mask,ENTRY *X,ENTRY *Y,ENTRY *AX,unsigned short *AX_mask,ENTRY *D,ENTRY **filter,int *order)
{
	int tid = omp_get_thread_num();
	//struct GeomInfo *geom_info;
	//struct ImgInfo *img_info;
	//struct SourceLocInfo *source_loc_info;
	//struct PriorInfo *prior_info;
	//char **recon_mask;
	//ENTRY *X;
	//ENTRY *Y;
	//ENTRY *AX;
	//unsigned short *AX_mask;
	//ENTRY *D;  /* sjk */
	//ENTRY **filter;
	//int *order;
	//pthread_mutex_t	*AX_mutex;  /* sjk */
	//struct paraICDReconstructData *data;

	int offset;  /* sjk */
	int j, jx, jy, jz, jzmax, Nxy, Nyz;
	float x, y, z, pixel, diff;
	struct ACol col_xyz;
	struct ViewXYInfo view_xy_info;
	struct ViewXYZInfo view_xyz_info;
	struct ICDInfo icd_info;

	struct timeval start, end;

	int k, l, Nsubit, jjx, jjy;
	float gamma;
	char **check;
	ENTRY **UMM;
	ENTRY **VSC;
	/* sjk */
	int *voxel_list;
	float *key;
	char zero_skip_flag;  
	int n;  

	//data = (struct paraICDReconstructData *)input;
	//tid = data->tid;
	//geom_info = data->geom_info;
	//img_info = data->img_info;
	//source_loc_info = data->source_loc_info;
	//prior_info = data->prior_info;
	//recon_mask = data->recon_mask;
	//X = data->X;
	//Y = data->Y;
	//AX = data->AX;
	//AX_mask = data->AX_mask;
	//D = data->D;
	//filter = data->filter;
	//order = data->order;
	//AX_mutex = data->AX_mutex;  /* sjk */

	gettimeofday(&start, NULL);

	createACol(&col_xyz, COL_LEN);		/* TODO COL_LEN hard-coded */
	createViewXYInfo(&view_xy_info, geom_info);
	createViewXYZInfo(&view_xyz_info, geom_info);

	check = (char **)get_img(img_info->Ny, img_info->Nx, sizeof(char));
	UMM = (ENTRY **)get_img(img_info->Ny, img_info->Nx, sizeof(ENTRY));
	VSC = (ENTRY **)get_img(img_info->Ny, img_info->Nx, sizeof(ENTRY));
	/* sjk */
        voxel_list = (int *)get_spc(img_info->Nx*img_info->Ny,sizeof(int));
        key = (float *)get_spc(img_info->Nx*img_info->Ny,sizeof(float));

	icd_info.q = prior_info->q;
	icd_info.p = prior_info->p;
	icd_info.c = prior_info->c;
	icd_info.sigma = prior_info->sigma;
	icd_info.sigma2q = pow(prior_info->sigma, prior_info->q);
	icd_info.c2qmp = pow(prior_info->c, (prior_info->q - prior_info->p) );

	/* sjk: compute the distance-weights for the neighbors */
        compNeighborWeight(icd_info.nb_wt,img_info);

	jzmax = (tid+1)*img_info->Nz/omp_get_num_threads();
	if (tid == (omp_get_num_threads()-1))
	{
		jzmax = img_info->Nz;
	}
	
	Nxy = img_info->Nx*img_info->Ny;
	Nyz = img_info->Ny*img_info->Nz;


	/* homogeneous ICD */
	for (j = 0; j < Nxy; j++)
	{
		jx = order[j] / img_info->Ny;
		jy = order[j] % img_info->Ny;

		UMM[jx][jy] = 0.0;

		if (recon_mask[jx][jy])
		{
			x = img_info->x0 + jx*img_info->Del_xy;
			y = img_info->y0 + jy*img_info->Del_xy;
			compViewXYInfo(x, y, &view_xy_info, geom_info, img_info, source_loc_info);

			offset= jx*Nyz + jy*img_info->Nz;  /* sjk */
			for (jz = tid*img_info->Nz/omp_get_num_threads(); jz < jzmax; jz++)
			{
				/* sjk: moved this here to facilitate zero-skip check */
				fillNeighbors(icd_info.neighbors,jx,jy,jz,X,img_info); 

				/*** sjk: zero skip check ***/
				if( X[offset+jz] > 0)
				  zero_skip_flag=0;  /* don't skip */
				else
				{
				  zero_skip_flag=1;  /* maybe skip.. check neighbors */
				  for(n=0; n<26; n++)
				  {
				    if( icd_info.neighbors[n] > 0)
				    {
				      zero_skip_flag=0;
				      break;
				    }
				  }
				}

				if(zero_skip_flag==0)
				{

				z = img_info->z0 + jz*img_info->Del_z;

				compViewXYZInfo(z, &view_xyz_info, geom_info, img_info, source_loc_info, &view_xy_info);

				/* calculate A column on the fly*/
				compAColxyzOnFly(x, y, z, geom_info, source_loc_info, &view_xy_info, &view_xyz_info, &col_xyz);

				if (col_xyz.n_index > 0)
				{
					/* ICD update */
/*** for non-homogeneous prior, compute weights here ***/
/*
compNeighborWeightNonhomogeneous(icd_info.nb_wt,img_info,geom_info,source_loc_info,x,y,z);
*/

					pixel = ICDStep(&icd_info, prior_info, jx, jy, jz, X, Y, AX, AX_mask, D, geom_info->lambda0, img_info, &col_xyz);

					/* clip */
					X[offset+jz] = ((pixel < 0.0) ? 0.0 : pixel);  /* sjk */
/* sjk: try without positivity constraint */
if (recon_mask[jx][jy]==2)
	X[offset+jz] = pixel; 

					/* update AX */
					diff = X[offset+jz]-icd_info.v;  /* sjk */
					/* sjk: lock update of AX with a mutex to prevent write conflicts between threads */
					//pthread_mutex_lock(AX_mutex);
					#pragma omp critical
					{
					updateAX(AX, &(col_xyz), diff);
					}
					//pthread_mutex_unlock(AX_mutex);
					UMM[jx][jy] += abso(diff);
				}
				}  /* END zero_skip_flag */
			}

			freeViewXYInfoB(&view_xy_info, geom_info);
		}
	}

	/* NH-ICD */
	Nsubit = 10;
	for (k = 0; k < Nsubit; k++)
	{
		filterUMM(UMM, VSC, filter, img_info->Nx, img_info->Ny);

		gamma = 0.05;

		/* sjk: using quickselect() to find top VSC voxels */
                for (j = 0; j < Nxy; j++)
                {
                        jx = order[j] / img_info->Ny;
                        jy = order[j] % img_info->Ny;
                        key[j]=(float)(-VSC[jx][jy]);
                        voxel_list[j]=order[j];
                }
                /*quicksort(key,voxel_list,0,Nxy-1);*/
                quickselect(key,voxel_list,0,Nxy-1,(int)ceil(gamma*Nxy));

		for (l = 0; l < (int)ceil(gamma*Nxy); l++)
		{

			/* sjk */
                        jjx = voxel_list[l] / img_info->Ny;
                        jjy = voxel_list[l] % img_info->Ny;

			if (recon_mask[jjx][jjy])
			{
				UMM[jjx][jjy] = 0.0;
				x = img_info->x0 + jjx*img_info->Del_xy;
				y = img_info->y0 + jjy*img_info->Del_xy;
				compViewXYInfo(x, y, &view_xy_info, geom_info, img_info, source_loc_info);

				offset= jjx*Nyz + jjy*img_info->Nz;  /* sjk */
				for (jz = tid*img_info->Nz/omp_get_num_threads(); jz < jzmax; jz++)
				{
					/* sjk: moved this here to facilitate zero-skip check */
					fillNeighbors(icd_info.neighbors,jjx,jjy,jz,X,img_info); 

				        /*** sjk: zero skip check ***/
					if( X[offset+jz] > 0)
					  zero_skip_flag=0;  /* don't skip */
					else
					{
					  zero_skip_flag=1;  /* maybe skip.. check neighbors */
					  for(n=0; n<26; n++)
					  {
  					    if( icd_info.neighbors[n] > 0)
  					    {
    					      zero_skip_flag=0;
    					      break;
  					    }
					  }
					}

					if(zero_skip_flag==0)
					{

					z = img_info->z0 + jz*img_info->Del_z;

					compViewXYZInfo(z, &view_xyz_info, geom_info, img_info, source_loc_info, &view_xy_info);

					/* calculate A column on the fly*/
					compAColxyzOnFly(x, y, z, geom_info, source_loc_info, &view_xy_info, &view_xyz_info, &col_xyz);
					if (col_xyz.n_index > 0)
					{
						/* ICD update */
						pixel = ICDStep(&icd_info, prior_info, jjx, jjy, jz, X, Y, AX, AX_mask, D, geom_info->lambda0, img_info, &col_xyz);

						/* clip */
						X[offset+jz] = ((pixel < 0.0) ? 0.0 : pixel);  /* sjk */
						/* sjk: try without positivity constraint */
						if (recon_mask[jx][jy]==2)
							X[offset+jz] = pixel; 

						/* update AX */
						diff = X[offset+jz]-icd_info.v;  /* sjk */

                                        	#pragma omp critical
                                        	{
							updateAX(AX, &(col_xyz), diff);
						}
						UMM[jjx][jjy] += abso(diff);
					}
					}  /* END zero_skip_flag */
				}

				freeViewXYInfoB(&view_xy_info, geom_info);
			}
		}
	}

	freeViewXYZInfo(&view_xyz_info);
	freeViewXYInfo(&view_xy_info, geom_info);
	freeACol(&col_xyz);
	/* sjk */
	free((void *)voxel_list);
	free((void *)key);
	free_img((void **)check);
	free_img((void **)UMM);
	free_img((void **)VSC);

	gettimeofday(&end, NULL);
	/*fprintf(stdout, "thread %d: time = %ld secs\n", tid, (end.tv_sec - start.tv_sec));*/

	return 0;
}

float quadUpdate(struct ICDInfo *icd_info)
{
	int i;
	float num,denom,sum1,sum2;

	sum1=sum2=0;
	for (i = 0; i < 26; i++)
	{
		sum1 += icd_info->nb_wt[i]*icd_info->neighbors[i];
		sum2 += icd_info->nb_wt[i];
	}

	num = sum1/icd_info->sigma2q + icd_info->th2*icd_info->v - icd_info->th1;
	denom = icd_info->th2 + sum2/icd_info->sigma2q;

	return num/denom;

}

/* unified q-GMMRF potential function */
/* including GGMRF, quad */
float qpotential(float delta, float q, float p, float c)
{
	float ppot;
	float num;

	ppot = pow(fabs(delta), p);

	if (q == p)
	{
		/* p-GGMRF */
		return ppot;
	}

	num = pow(fabs(delta/c),q-p);

	/* q-GGMRF */
	return (ppot*num/(1.0+num));
}

/* not q-GMMRF influence function */
/* but rho'(Delta)/Delta */
float qinfluence(float delta, float q, float p, float c)
{
	float qp;
	float fact;
	float denom;
	float sum1, sum2;

	qp = q - p;
	fact = p*pow(c, qp);

	denom = 1.0 + pow(fabs(delta/c), qp);

	sum1 = 2.0/denom;
	sum2 = fabs(delta)*qp*pow(fabs(delta), 1.0-p)/(pow(c, qp)*denom*denom);

	return (sum1 - sum2)/fact;
}

/* substitute function update for QGGMRF.. assumes q=2.0 */
float subFunctionUpdate(struct ICDInfo *icd_info)
{
	int i;
	float num, denom, term, g;
	float delta;
	float v,c,p,c2qmp;

        v=icd_info->v;
        c=icd_info->c;
        p=icd_info->p;
        c2qmp=icd_info->c2qmp;

	num=0;
	denom=0;
	for (i = 0; i < 26; i++)
	{
		delta = (v - icd_info->neighbors[i]);
	        term = 1.0 + pow(fabs(delta/c),2.0-p);
		g = (2.0 - (2.0-p)*(term-1.0)/term) / term;

		num += icd_info->nb_wt[i] * g * icd_info->neighbors[i];
		denom += icd_info->nb_wt[i] * g;
	}
	num = num/(2.0* icd_info->c2qmp * icd_info->sigma2q) + icd_info->th2*icd_info->v - icd_info->th1;
	denom = denom/(2.0* icd_info->c2qmp * icd_info->sigma2q) + icd_info->th2;

	return(num/denom);
}

float pCost(float x, struct ICDInfo *icd_info)
{
	int i;
	float tmp,sum;

	sum=0;
	for (i = 0; i < 26; i++)
		sum += icd_info->nb_wt[i]*pow(fabs(x-icd_info->neighbors[i]), icd_info->p);

	tmp = x - icd_info->v;

	return (icd_info->th1*tmp + icd_info->th2*tmp*tmp/2.0 + sum/(icd_info->p*icd_info->sigma2q));

}

float pCostDerivative(float x, void *data)
{
	int i;
	struct ICDInfo *icd_info;
	float sum;

	icd_info = (struct ICDInfo *)data;

	sum=0;
	for (i = 0; i < 26; i++)
		sum += icd_info->nb_wt[i] * pow(fabs(x-icd_info->neighbors[i]), (icd_info->p-1))*sign(x,icd_info->neighbors[i]);

	return (icd_info->th1 + icd_info->th2*(x-icd_info->v) + sum/icd_info->sigma2q);

}

float halfIntervalUpdate(struct ICDInfo *icd_info)
{
	int status;
	float pixel;

	/* solve rooting equation */
	pixel = solve(pCostDerivative, icd_info->low, icd_info->high, TOL, &status, (void *)icd_info);

	if (status)
	{
		/* fprintf(stdout, "rooting equation has no solution!\n"); */
		if (pCost(icd_info->low, icd_info) < pCost(icd_info->high, icd_info))
		{
			pixel = icd_info->low;
		}
		else
		{
			pixel = icd_info->high;
		}
	}

	return pixel;
}
