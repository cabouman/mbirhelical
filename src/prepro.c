/* version 6.0, P.Jin 07/13/2012 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "data.h"
#include "allocate.h"
#include "proj.h"
#include "randlib.h"
#include "prepro.h"

void constImage(ENTRY *X, char **recon_mask, float value, struct ImgInfo *img_info)
{
	int j, jx, jy;
	int Nyz;
	int tmp;

	Nyz = img_info->Ny*img_info->Nz;

	for (j = 0; j < img_info->Nx*Nyz; j++)
	{
		jx = j/Nyz;
		tmp = j - jx*Nyz;
		jy = tmp/img_info->Nz;

		/* pixels out of the reconstruction scope are set to be 0 */
		if (recon_mask[jx][jy])
		{
			X[j] = value;
		}
		else
		{
			X[j] = 0.0;
		}
	}
}

void clipImage(ENTRY *X, char **recon_mask, struct ImgInfo *img_info)
{
	int j, jx, jy;
	int Nyz,Nxyz;
	int tmp;

	Nyz = img_info->Ny*img_info->Nz;
	Nxyz = Nyz*img_info->Nx;

	for (j = 0; j < Nxyz; j++)
	{
		jx = j/Nyz;
		tmp = j - jx*Nyz;
		jy = tmp/img_info->Nz;

		/* pixels out of the reconstruction scope are set to be 0 */
		/* pixels whose value is negative are set to be 0 */
//		if (recon_mask[jx][jy] == 0 || X[j] < 0.0)
//		{
//			X[j] = 0;
//		}
//
		if (recon_mask[jx][jy] == 0)
		{
			X[j] = 0;
		}

	}
}

void addWhiteNoise(ENTRY *data, int len, float noise_std)
{
	int i;

	fprintf(stdout, "\nadding white noise...\n");
	for (i = 0; i < len; i++)
	{
		data[i] += noise_std*normal();
	}
}

void addNoise(ENTRY *data, float lambda0, int len)
{
	int i;
	float cnt;

	fprintf(stdout, "\nadding noise...\n");
	for (i = 0; i < len; i++)
	{
		cnt = lambda0*exp(-data[i]);
		data[i] += sqrt(1.0/cnt)*normal();
	}
}

void createReconMask(char ***recon_mask, struct ImgInfo *img_info)
{
	*recon_mask = (char **)get_img(img_info->Ny, img_info->Nx, sizeof(char));
}

void compReconMask(char **recon_mask, struct ImgInfo *img_info)
{
	int jx, jy;
	float x, y;
	float r_sq;		/* x^2 + y^2 */
	float rI_sq;		/* rI^2 */

	rI_sq = img_info->rI*img_info->rI;

	for (jx = 0; jx < img_info->Nx; jx++)
	for (jy = 0; jy < img_info->Ny; jy++)
	{
		x = img_info->x0 + jx*img_info->Del_xy;
		y = img_info->y0 + jy*img_info->Del_xy;
		r_sq = x*x + y*y;

		recon_mask[jx][jy] = ((r_sq > rI_sq) ? 0 : 1);
	}
}

void freeReconMask(char **recon_mask)
{
	free_img((void *)recon_mask);
}

void initImage(
	struct GeomInfo *geom_info,
	struct Image *image,
	ENTRY value)
{
	int i, n_ext, len;
	float d;

	/* if only the number of good slices is specified, add buffer slices */
	if (image->img_info.Nz == 0)
	{
		d = geom_info->Nr*geom_info->Del_dr*(geom_info->r_si+geom_info->fov/2.0)/(2.0*geom_info->r_sd);
		n_ext = (int)ceil(d/image->img_info.Del_z);
		image->img_info.Nz = image->img_info.Nz_mid + 2*n_ext;
	}
	len = (image->img_info.Nx)*(image->img_info.Ny)*(image->img_info.Nz);
	image->img = (ENTRY *)get_spc(len, sizeof(ENTRY));

	for (i = 0; i < len; i++)
	{
		image->img[i] = value;
	}
}

/* randomly shuffle index array, generate randomized index order */
void shuffle(int *array, int len)
{
	int i, j, tmp;

	int sum1, sum2;

	sum1 = 0;
	for (i = 0; i < len; i++)
	{
		array[i] = i;
		sum1 += i;
	}

	srand(time(NULL));
	for (i = 0; i < len-1; i++)
	{
		j = i + rand()%(len-i);
		tmp = array[j];
		array[j] = array[i];
		array[i] = tmp;
	}

	sum2 = 0;
	for (i = 0; i < len; i++)
	{
		sum2 += array[i];
	}

	if(sum1 != sum2)
	{
		fprintf(stderr, "shuffle has bugs!!\n");
		exit(1);
	}
}
/*
void filterUMM(ENTRY **UMM, ENTRY **VSC, ENTRY **filter, int h, int w)
{
// use 5x5 Hamming window 

	int i, j, m, n, p, q;

	for (i = 0; i < h; i++)
	{
		for (j = 0; j < w; j++)
		{
			VSC[i][j] = 0.0;
			for (m = -2; m <= 2; m++)
			{
				for (n = -2; n <= 2; n++)
				{
					p = i - m;
					q = j - n;
					if (p >= 0 && p < h && q >= 0 && q < w)
					{
						VSC[i][j] += filter[m+2][n+2]*UMM[p][q];
					}
				}
			}
		}
	}
}

*/

/* TODO */
/*
void constDosage(ENTRY *LambdaT, float dos, struct GeomInfo *geom_info)
{
	int i;
	for (i = 0; i < (geom_info->Nv*geom_info->Nc*geom_info->Nr); i++)
	{
		LambdaT[i] = dos;
	}
}
*/
/*
void simulatePhotonCount(ENTRY *Lambda, ENTRY *LambdaT, ENTRY *AX, struct GeomInfo *geom_info)
{
	int i;
	float mean;

	FILE *fp;
	fp = fopen("debug.log", "w");

	for (i = 0; i < (geom_info->Nv*geom_info->Nc*geom_info->Nr); i++)
	{
		mean = LambdaT[i]*exp(-AX[i]);
		Lambda[i] = mean + normal()*sqrt(mean);
		if (Lambda[i] <= 0.0)
		{
			Lambda[i] = 10e-5;
		}
		if (AX[i] > 0.0)
		{
		fprintf(fp, "AX = %f\n", AX[i]);
		fprintf(fp, "mean = %f\n", mean);
		fprintf(fp, "count = %f\n", Lambda[i]);
		}
	}
	fclose(fp);
}
*/
/*
void compProjection(ENTRY *Y, ENTRY *Lambda, ENTRY *LambdaT, struct GeomInfo *geom_info)
{
	int i;

	for (i = 0; i < (geom_info->Nv*geom_info->Nc*geom_info->Nr); i++)
	{
		Y[i] = -log(Lambda[i]/LambdaT[i]);
	}
}
*/
