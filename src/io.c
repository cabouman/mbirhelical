/* version 6.0, P.Jin 07/13/2012 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>	/* readPrior */
#include "allocate.h"
#include "data.h"
#include "prepro.h"	/* readPrior */
#include "io.h"

/* byte-order possibilities */
/*
#define LittleEndian    0x4949
#define BigEndian   0x4d4d

#define is_bigendian() ( (*(char*)&check_i) == 0 )
#define HostByteOrder   ( is_bigendian() ? BigEndian : LittleEndian)
unsigned short int FileByteOrder = BigEndian;
const unsigned short int check_i=1;
*/

void readImgInfo(char *fname, struct ImgInfo *img_info)
{
	FILE *fp;
	char tag[100];

	if ((fp = fopen(fname, "r")) == NULL)
	{
		fprintf(stderr, "ERROR in readImgInfo: can't open file %s\n", fname);
		exit(1);
	}

	fgets(tag, 100, fp);
	fscanf(fp, "%d\n\n", &(img_info->Nx));
	fgets(tag, 100, fp);
	fscanf(fp, "%d\n\n", &(img_info->Ny));
	fgets(tag, 100, fp);
	fscanf(fp, "%d\n\n", &(img_info->Nz_mid));
	fgets(tag, 100, fp);
	fscanf(fp, "%d\n\n", &(img_info->Nz));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(img_info->xc));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(img_info->yc));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(img_info->zc));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(img_info->Del_xy));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(img_info->Del_z));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(img_info->rI));
	fgets(tag, 100, fp);
	fscanf(fp, "%s\n\n", img_info->imgFile);
	fgets(tag, 100, fp);
	fscanf(fp, "%s\n", img_info->maskFile);

	fclose(fp);
}

void writeImgInfo(char *fname, struct ImgInfo *img_info)
{
	FILE *fp;

	if ((fp = fopen(fname, "w")) == NULL)
	{
		fprintf(stderr, "ERROR in writeImgInfo: can't open file %s\n", fname);
		exit(1);
	}

	fprintf(fp, "number of voxels in x\n");
	fprintf(fp, "%d\n\n", img_info->Nx);
	fprintf(fp, "number of voxels in y\n");
	fprintf(fp, "%d\n\n", img_info->Ny);
	fprintf(fp, "number of voxels in z (good slices)\n");
	fprintf(fp, "%d\n\n", img_info->Nz_mid);
	fprintf(fp, "number of voxels in z (total)\n");
	fprintf(fp, "%d\n\n", img_info->Nz);
	fprintf(fp, "x coordinate of the center voxel (mm)\n");
	fprintf(fp, "%f\n\n", img_info->xc);
	fprintf(fp, "y coordinate of the center voxel (mm)\n");
	fprintf(fp, "%f\n\n", img_info->yc);
	fprintf(fp, "z coordinate of the center voxel (mm)\n");
	fprintf(fp, "%f\n\n", img_info->zc);
	fprintf(fp, "voxel spacing in xy (mm)\n");
	fprintf(fp, "%f\n\n", img_info->Del_xy);
	fprintf(fp, "voxel spacing in z (mm)\n");
	fprintf(fp, "%f\n\n", img_info->Del_z);
	fprintf(fp, "radius of the reconstruction mask (mm)\n");
	fprintf(fp, "%f\n\n", img_info->rI);
	fprintf(fp, "initial reconstruction image location\n");
	fprintf(fp, "%s\n", img_info->imgFile);

	fclose(fp);

}

void printImgInfo(struct ImgInfo *img_info)
{
	fprintf(stdout, "\nIMAGE PARAMETERS:\n");

	fprintf(stdout, "number of voxels in x: ");
	fprintf(stdout, "%d\n", img_info->Nx);
	fprintf(stdout, "number of voxels in y: ");
	fprintf(stdout, "%d\n", img_info->Ny);
	fprintf(stdout, "number of voxels in z (good slices): ");
	fprintf(stdout, "%d\n", img_info->Nz_mid);
	fprintf(stdout, "number of voxels in z (total): ");
	fprintf(stdout, "%d\n", img_info->Nz);
	fprintf(stdout, "x coordinate of the center voxel (mm): ");
	fprintf(stdout, "%f\n", img_info->xc);
	fprintf(stdout, "y coordinate of the center voxel (mm): ");
	fprintf(stdout, "%f\n", img_info->yc);
	fprintf(stdout, "z coordinate of the center voxel (mm): ");
	fprintf(stdout, "%f\n", img_info->zc);
	fprintf(stdout, "voxel spacing in xy (mm): ");
	fprintf(stdout, "%f\n", img_info->Del_xy);
	fprintf(stdout, "voxel spacing in z (mm): ");
	fprintf(stdout, "%f\n", img_info->Del_z);
	fprintf(stdout, "radius of the reconstruction mask (mm): ");
	fprintf(stdout, "%f\n", img_info->rI);
	fprintf(stdout, "initial reconstruction image location: ");
	fprintf(stdout, "%s\n", img_info->imgFile);
}

void readImage(char *fname, struct Image *image)
{
	FILE *fp;
	int i, Nx, Ny, Nz, len;
	ENTRY pixel;

	if ((fp = fopen(fname, "rb")) == NULL)
	{
		fprintf(stderr, "ERROR in readImage: can't open file %s\n", fname);
		exit(1);
	}

	fscanf(fp, "%d %d %d\n", &Nx, &Ny, &Nz);
	if (Nx != image->img_info.Nx || Ny != image->img_info.Ny || Nz != image->img_info.Nz)
	{
		fprintf(stderr, "ERROR in readImage: dimension does not match.\n");
		exit(1);
	}

	if (sizeof(ENTRY) == sizeof(float))
	{
		/* allocate memory to store image data */
		len = (image->img_info.Nx)*(image->img_info.Ny)*(image->img_info.Nz);
		image->img = (ENTRY *)get_spc(len, sizeof(ENTRY));

		/* TODO check endianness before read ? */
		for (i = 0; i < len; i++)
		{
			fscanf(fp, "%f ", &pixel);
			image->img[i] = hu2miu(pixel, MIU_AIR, MIU_WATER);	/* change to mm-1 */
		}
	}
	else
	{
		fprintf(stderr, "ERROR in readImage: unknown type error\n");
		exit(1);
	}

	fclose(fp);
}

void writeImage(char *fname, struct Image *image)
{
	FILE *fp;
	int i, len;
	ENTRY pixel;

	if ((fp = fopen(fname, "wb")) == NULL)
	{
		fprintf(stderr, "ERROR in writeImage: can't open file %s\n", fname);
		exit(1);
	}

	fprintf(fp, "%d %d %d\n", image->img_info.Nx, image->img_info.Ny, image->img_info.Nz);

	if (sizeof(ENTRY) == sizeof(float))
	{
		len = (image->img_info.Nx)*(image->img_info.Ny)*(image->img_info.Nz);

		/* TODO check endianness before write ? */
		for (i = 0; i < len; i++)
		{
			pixel = miu2hu(image->img[i], MIU_AIR, MIU_WATER);	/* change to HU */
			if(pixel==0)
				fprintf(fp, "0 ");
			else
				fprintf(fp, "%f ", pixel);
		}
	}
	else
	{
		fprintf(stderr, "ERROR in writeImage: unknown type error\n");
		exit(1);
	}

	fclose(fp);
}


void readImage_short(char *fname, struct Image *image)
{
	FILE *fp;
	int i, Nx, Ny, Nz, len,len2;
	float Dx,Dy,Dz,x0,y0,z0;
	unsigned short *input;

	if ((fp = fopen(fname, "r")) == NULL)
	{
		fprintf(stderr, "ERROR in readImage: can't open file %s\n", fname);
		exit(1);
	}

	fscanf(fp, "%d %d %d", &Nx, &Ny, &Nz);
/* this is new */
	fscanf(fp, "%f %f %f",&Dx,&Dy,&Dz);
	fscanf(fp, "%f %f %f",&x0,&y0,&z0);
	fseek(fp,1,SEEK_CUR);

	if (Nx != image->img_info.Nx || Ny != image->img_info.Ny || Nz != image->img_info.Nz)
	{
		fprintf(stderr, "ERROR in readImage: dimension does not match.\n");
		exit(1);
	}

	/* allocate memory to store image data */
	len = (image->img_info.Nx)*(image->img_info.Ny)*(image->img_info.Nz);
	image->img = (ENTRY *)get_spc(len, sizeof(ENTRY));
	input = (unsigned short *)get_spc(len, sizeof(unsigned short));

	len2=fread(input,sizeof(unsigned short),len,fp);
	if( len2 != len )
	{
		fprintf(stderr, "ERROR in readImage: file terminated early.\n");
		fprintf(stderr, "file size=%d, expected size=%d\n",len2,len);
		fclose(fp);
		exit(1);
	}

	for (i = 0; i < len; i++)
		image->img[i] = hu2miu((ENTRY)input[i], MIU_AIR, MIU_WATER);	/* change to mm-1 */

	fclose(fp);
	free(input);
}


void writeImage_short(char *fname, struct Image *image)
{
	FILE *fp;
	int i, len;
	unsigned short *output;
	ENTRY pixel;

	if ((fp = fopen(fname, "w")) == NULL)
	{
		fprintf(stderr, "ERROR in writeImage: can't open file %s\n", fname);
		exit(1);
	}

	fprintf(fp, "%03d %03d %03d\n", image->img_info.Nx, image->img_info.Ny, image->img_info.Nz);
/* this is new */
	fprintf(fp, "%f %f %f\n", image->img_info.Del_xy, image->img_info.Del_xy, image->img_info.Del_z);
	fprintf(fp, "%f %f %f\n", image->img_info.x0, image->img_info.y0, image->img_info.z0);

	len = (image->img_info.Nx)*(image->img_info.Ny)*(image->img_info.Nz);

	/* allocate memory for temporary storage */
	output = (unsigned short *)get_spc(len, sizeof(unsigned short));

	for (i = 0; i < len; i++)
	{
		pixel = max(0,miu2hu(image->img[i], MIU_AIR, MIU_WATER));	/* change to mm-1 */
		output[i] = (unsigned short)(pixel + 0.5);	/* round to integer and cast */
	}

	if( fwrite(output,sizeof(unsigned short),len,fp) != len )
	{
		fprintf(stderr, "ERROR in readImage: file not written successfully.\n");
		fclose(fp);
		exit(1);
	}

	fclose(fp);
	free(output);
}


void readReconMask(char **recon_mask, struct ImgInfo *img_info)
{
	FILE *fp;
	int jx,jy, Nx, Ny, Nz, temp;

	if ((fp = fopen(img_info->maskFile, "r")) == NULL)
	{
		fprintf(stderr, "ERROR in readReconMask: can't open file %s\n", img_info->maskFile);
		exit(1);
	}

	fscanf(fp, "%d %d %d\n", &Nx, &Ny, &Nz);
	if (Nx != img_info->Nx || Ny != img_info->Ny)
	{
		fprintf(stderr, "ERROR in readReconMask: dimension does not match.\n");
		exit(1);
	}

	/* allocate memory to store image data */
	for (jx = 0; jx < img_info->Nx; jx++)
	for (jy = 0; jy < img_info->Ny; jy++)
	{
		fscanf(fp, "%d ",&temp);
		recon_mask[jx][jy]=temp;
	}
	
	fclose(fp);
}


void readGeomInfo(char *fname, struct GeomInfo *geom_info)
{
	FILE *fp;
	char tag[100];

	if ((fp = fopen(fname, "r")) == NULL)
	{
		fprintf(stderr, "ERROR in readGeomInfo: can't open file %s\n", fname);
		exit(1);
	}

	fgets(tag, 100, fp);
	fscanf(fp, "%d\n\n", &(geom_info->Nr));
	fgets(tag, 100, fp);
	fscanf(fp, "%d\n\n", &(geom_info->Nc));
	fgets(tag, 100, fp);
	fscanf(fp, "%d\n\n", &(geom_info->Nv));
	fgets(tag, 100, fp);
	fscanf(fp, "%d\n\n", &(geom_info->Nvpr));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(geom_info->r_si));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(geom_info->r_sd));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(geom_info->u));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(geom_info->beta0));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(geom_info->Del_beta));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(geom_info->Del_alphac));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(geom_info->del_alphac));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(geom_info->Del_dr));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(geom_info->fov));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(geom_info->lambda0));
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &(geom_info->evar));
	fgets(tag, 100, fp);
	fscanf(fp, "%s\n\n", geom_info->sinoFile);
	fgets(tag, 100, fp);
	fscanf(fp, "%s\n\n", geom_info->doseFile);
	fgets(tag, 100, fp);
	fscanf(fp, "%s\n\n", geom_info->offsetFile);
	fgets(tag, 100, fp);
	fscanf(fp, "%s\n", geom_info->detectorsFile);

	fclose(fp);
}

void writeGeomInfo(char *fname, struct GeomInfo *geom_info)
{
	FILE *fp;

	if ((fp = fopen(fname, "r")) == NULL)
	{
		fprintf(stderr, "ERROR in readGeomInfo: can't open file %s\n", fname);
		exit(1);
	}

	fprintf(fp, "number of rows\n");
	fprintf(fp, "%d\n\n", geom_info->Nr);
	fprintf(fp, "number of channels\n");
	fprintf(fp, "%d\n\n", geom_info->Nc);
	fprintf(fp, "number of views\n");
	fprintf(fp, "%d\n\n", geom_info->Nv);
	fprintf(fp, "views per rotation\n");
	fprintf(fp, "%d\n\n", geom_info->Nvpr);
	fprintf(fp, "src to iso (mm)\n");
	fprintf(fp, "%f\n\n", geom_info->r_si);
	fprintf(fp, "src to det (mm)\n");
	fprintf(fp, "%f\n\n", geom_info->r_sd);
	fprintf(fp, "un-normalized pitch (det rows/rot)\n");
	fprintf(fp, "%f\n\n", geom_info->u);
	fprintf(fp, "initial view angle (rad)\n");
	fprintf(fp, "%f\n\n", geom_info->beta0);
	fprintf(fp, "view angle spacing (rad)\n");
	fprintf(fp, "%f\n\n", geom_info->Del_beta);
	fprintf(fp, "detector channel angle spacing (rad)\n");
	fprintf(fp, "%f\n\n", geom_info->Del_alphac);
	fprintf(fp, "detector channel offset (rad)\n");
	fprintf(fp, "%f\n\n", geom_info->del_alphac);
	fprintf(fp, "detector row width (mm)\n");
	fprintf(fp, "%f\n\n", geom_info->Del_dr);
	fprintf(fp, "diameter of field of view (mm)\n");
	fprintf(fp, "%f\n\n", geom_info->fov);
	fprintf(fp, "initial photon counts (0 for noiseless)\n");
	fprintf(fp, "%f\n\n", geom_info->lambda0);
	fprintf(fp, "electronic noise variance (0 for noiseless)\n");
	fprintf(fp, "%f\n\n", geom_info->evar);
	fprintf(fp, "sinogram location\n");
	fprintf(fp, "%s\n", geom_info->sinoFile);
	fprintf(fp, "dosage location\n");
	fprintf(fp, "%s\n", geom_info->doseFile);
	fprintf(fp, "offset data location\n");
	fprintf(fp, "%s\n", geom_info->offsetFile);

	fclose(fp);
}

void printGeomInfo(struct GeomInfo *geom_info)
{
	fprintf(stdout, "\nGEOMETRY PARAMETERS:\n");

	fprintf(stdout, "number of rows: ");
	fprintf(stdout, "%d\n", geom_info->Nr);
	fprintf(stdout, "number of channels: ");
	fprintf(stdout, "%d\n", geom_info->Nc);
	fprintf(stdout, "number of views: ");
	fprintf(stdout, "%d\n", geom_info->Nv);
	fprintf(stdout, "views per rotation: ");
	fprintf(stdout, "%d\n", geom_info->Nvpr);
	fprintf(stdout, "src to iso (mm): ");
	fprintf(stdout, "%f\n", geom_info->r_si);
	fprintf(stdout, "src to det (mm): ");
	fprintf(stdout, "%f\n", geom_info->r_sd);
	fprintf(stdout, "un-normalized pitch (det rows/rot): ");
	fprintf(stdout, "%f\n", geom_info->u);
	fprintf(stdout, "initial view angle (rad): ");
	fprintf(stdout, "%f\n", geom_info->beta0);
	fprintf(stdout, "view angle spacing (rad): ");
	fprintf(stdout, "%f\n", geom_info->Del_beta);
	fprintf(stdout, "detector channel angle spacing (rad): ");
	fprintf(stdout, "%f\n", geom_info->Del_alphac);
	fprintf(stdout, "detector channel offset (rad): ");
	fprintf(stdout, "%f\n", geom_info->del_alphac);
	fprintf(stdout, "detector row width (mm): ");
	fprintf(stdout, "%f\n", geom_info->Del_dr);
	fprintf(stdout, "diameter of field of view (mm): ");
	fprintf(stdout, "%f\n", geom_info->fov);
	fprintf(stdout, "initial photon counts (0 for noiseless): ");
	fprintf(stdout, "%f\n", geom_info->lambda0);
	fprintf(stdout, "electronic noise variance (0 for noiseless): ");
	fprintf(stdout, "%f\n", geom_info->evar);
	fprintf(stdout, "sinogram file location: ");
	fprintf(stdout, "%s\n", geom_info->sinoFile);
	fprintf(stdout, "dosage file location: ");
	fprintf(stdout, "%s\n", geom_info->doseFile);
	fprintf(stdout, "offset file location: ");
	fprintf(stdout, "%s\n", geom_info->offsetFile);
}


ENTRY* readSinogram_new(char *fname, int length)
{
	FILE *fp;
	int i, Nv, Nc, Nr;
	ENTRY *sino;

	if ((fp = fopen(fname, "rb")) == NULL)
	{
		fprintf(stderr, "ERROR in readSinogram: can't open file %s\n", fname);
		exit(1);
	}

	fscanf(fp, "%d %d %d\n", &Nr, &Nc, &Nv);
	if(length != Nr*Nc*Nv)
	{	
		fprintf(stderr, "ERROR in readSinogram: %s header doesn't match geometry.\n",fname);
		exit(-1);
	}

	sino = (ENTRY *)get_spc(length, sizeof(ENTRY));
	for (i = 0; i < length; i++)
		fscanf(fp, "%f ", sino+i);
	fclose(fp);

	return(sino);
}


ENTRY* readSinogram_float(char *fname, int length)
{
	FILE *fp;
	int Nv, Nc, Nr;
	char temp;
	ENTRY *sino;

	if ((fp = fopen(fname, "rb")) == NULL)
	{
		fprintf(stderr, "ERROR in readSinogram: can't open file %s\n", fname);
		exit(1);
	}

	fscanf(fp, "%d %d %d", &Nr, &Nc, &Nv);
	fread(&temp,1,1,fp);   /* skip past carriage return */
	if(length != Nr*Nc*Nv)
	{	
		fprintf(stderr, "ERROR in readSinogram: %s header doesn't match geometry.\n",fname);
		exit(-1);
	}

	sino = (ENTRY *)get_spc(length, sizeof(ENTRY));
	if( fread(sino,sizeof(ENTRY),length,fp) != length )
	{
		fprintf(stderr, "ERROR in readSinogram: file terminated early.\n");
		fclose(fp);
		exit(1);
	}
	return(sino);
}


void writeSinogram_float(char *fname, ENTRY *Y, int Nr, int Nc, int Nv)
{
        FILE *fp;
        int length;

        if ((fp = fopen(fname, "wb")) == NULL)
        {
                fprintf(stderr, "ERROR in writeSinogram_float: can't open file %s\n", fname);
                exit(1);
        }

        fprintf(fp, "%d %d %d\n", Nr, Nc, Nv);
        length=Nr*Nc*Nv;

        if( fwrite(Y,sizeof(ENTRY),length,fp) != length )
        {
                fprintf(stderr, "ERROR in writeSinogram_float: file not written successfully.\n");
                fclose(fp);
                exit(1);
        }
        fclose(fp);
}





/* reads/fills all available sinogram data (counts,dosage,etc.) */

void fillSinogramData(struct Sinogram *sinogram)
{
	int i,Nr,Nc,Nv,Nvpr,len,Nrc;
	float dose,counts,sum,offset;
	char offsetFlag=0;
	ENTRY *detector_mask;

	Nr=sinogram->geom_info.Nr;
	Nc=sinogram->geom_info.Nc;
	Nv=sinogram->geom_info.Nv;
	Nvpr=sinogram->geom_info.Nvpr;  /* views per rotation */
	Nrc = Nr*Nc;

	if (strcmp(sinogram->geom_info.detectorsFile, "NA") != 0)
		detector_mask= readSinogram_new(sinogram->geom_info.detectorsFile, Nrc);

	/* have dosage data */
	if (strcmp(sinogram->geom_info.doseFile, "NA") != 0)
	{
		sinogram->counts = readSinogram_float(sinogram->geom_info.sinoFile, Nr*Nc*Nv);
		sinogram->dose = readSinogram_float(sinogram->geom_info.doseFile, Nr*Nc*Nvpr);

		sinogram->geom_info.lambda0 = 1;   	/* basically used as flag in this case */

		if (strcmp(sinogram->geom_info.offsetFile, "NA") != 0)
		{
			sinogram->offset = readSinogram_float(sinogram->geom_info.offsetFile, Nr*Nc);
			offsetFlag=1;
		}

		/* fill in sinogram; compute D */
		len = Nr*Nc*Nv;
		sinogram->sino = (ENTRY *)get_spc(len, sizeof(ENTRY));
		sinogram->D = (ENTRY *)get_spc(len, sizeof(ENTRY));

		sum=0;
		for (i = 0; i < len; i++)
		{
			counts=sinogram->counts[i];
			dose = sinogram->dose[i%(Nr*Nc*Nvpr)];   /* dose[] contains 1 rotation's worth of data */
			if(offsetFlag)
			{
				offset=sinogram->offset[i%Nrc];
				counts -= offset;
				dose -= offset;
			}
			counts=max(1,counts);
			dose = max(1,dose);

			sinogram->sino[i]= log(dose/counts);
			if(dose >= counts)
				sinogram->D[i]= counts/dose;
			else
				sinogram->D[i]= 1.0;
			/*sinogram->D[i]= 1.0;*/
			/*sinogram->D[i]= 1.0/4.0 + 3.0/4.0*counts/dose;*/
			/*sinogram->D[i]= sqrt(counts/dose);*/
			/*sinogram->D[i]= pow(counts/dose,0.9);*/

			/* this is for weighting out "bad" detectors */
			if (strcmp(sinogram->geom_info.detectorsFile, "NA") != 0)
				if(detector_mask[i%Nrc]>0)
					sinogram->D[i]=0;
			sum += sinogram->D[i];
		}
		/* normalize so that trace(D)=trace(I) */
		for (i = 0; i < len; i++)
			sinogram->D[i] *= (len/sum);

		free(sinogram->counts);
		free(sinogram->dose);
		if(offsetFlag)
			free(sinogram->offset);
	}
	else  /* no dosage data */
	{
		sinogram->sino = readSinogram_float(sinogram->geom_info.sinoFile, Nr*Nc*Nv);

		/* fill in counts/dosage if assuming uniform dosage lambda0 */
		if(sinogram->geom_info.lambda0 > 0)
		{
			len = Nr*Nc*Nv;
			sinogram->D = (ENTRY *)get_spc(len, sizeof(ENTRY));
			/*sum=0; */
			for (i = 0; i < len; i++)
			{
				counts = sinogram->geom_info.lambda0 * exp(-sinogram->sino[i]);
				/*sinogram->D[i] = counts; */
				if( sinogram->geom_info.evar > 0 )
					sinogram->D[i] = counts/sinogram->geom_info.evar;
				else
					sinogram->D[i] = counts;
				/*sum += sinogram->D[i]; */
			}
			/* normalize so that trace(D)=trace(I) */
			
			/*for (i = 0; i < len; i++)
				sinogram->D[i] *= (len/sum);
			*/
		}
	}

}





void readSinogram(char *fname, struct Sinogram *sinogram)
{
	FILE *fp;
	int i, Nv, Nc, Nr, len;

	if ((fp = fopen(fname, "rb")) == NULL)
	{
		fprintf(stderr, "ERROR in readSinogram: can't open file %s\n", fname);
		exit(1);
	}

	fscanf(fp, "%d %d %d\n", &Nr, &Nc, &Nv);
	if (Nv != sinogram->geom_info.Nv || Nc != sinogram->geom_info.Nc || Nr != sinogram->geom_info.Nr)
	{
		fprintf(stderr, "ERROR in readSinogram: dimension does not match.\n");
		exit(1);
	}

	if (sizeof(ENTRY) == sizeof(float))
	{
		/* allocate memory to store image data */
		len = (sinogram->geom_info.Nr)*(sinogram->geom_info.Nc)*(sinogram->geom_info.Nv);
		sinogram->sino = (ENTRY *)get_spc(len, sizeof(ENTRY));
		/* sjk */
                if(sinogram->geom_info.lambda0 > 0)
                	sinogram->counts = (ENTRY *)get_spc(len, sizeof(ENTRY));

		/* TODO check endianness before read ? */
		/*fread(sinogram->sino, sizeof(ENTRY), len, fp);*/	/* maybe fread is faster */
		/**/
		for (i = 0; i < len; i++)
		{
			fscanf(fp, "%f ", sinogram->sino+i);
			/* sjk: precompute counts */
			if(sinogram->geom_info.lambda0 > 0)
				sinogram->counts[i]=(sinogram->geom_info.lambda0)*exp(-sinogram->sino[i]);
		}
		/**/
	}
	else
	{
		fprintf(stderr, "ERROR in readSinogram: unknown type error\n");
		exit(1);
	}

	fclose(fp);
}

void writeSinogram(char *fname, struct Sinogram *sinogram)
{
	FILE *fp;
	int i, len;

	if ((fp = fopen(fname, "wb")) == NULL)
	{
		fprintf(stderr, "ERROR in writeSinogram: can't open file %s\n", fname);
		exit(1);
	}

	fprintf(fp, "%d %d %d\n", sinogram->geom_info.Nr, sinogram->geom_info.Nc, sinogram->geom_info.Nv);

	if (sizeof(ENTRY) == sizeof(float))
	{
		/* allocate memory to store image data */
		len = (sinogram->geom_info.Nr)*(sinogram->geom_info.Nc)*(sinogram->geom_info.Nv);

		/* TODO check endianness before write ? */
		/*fwrite(sinogram->sino, sizeof(ENTRY), len, fp);*/	/* maybe fread is faster */
		/**/
		for (i = 0; i < len; i++)
		{
			if(sinogram->sino[i]==0)    /* sjk: can save a lot of bytes */
				fprintf(fp, "%d ", 0);
			else
				fprintf(fp, "%f ", sinogram->sino[i]);
		}
		/**/
	}
	else
	{
		fprintf(stderr, "ERROR in writeSinogram: unknown type error\n");
		exit(1);
	}

	fclose(fp);
}

void readPrior(char *fname, struct PriorInfo *prior_info)
{
	FILE *fp;
	float c;
	char tag[100];

	if ((fp = fopen(fname, "r")) == NULL)
	{
		fprintf(stderr, "ERROR in readPrior: can't open file %s\n", fname);
		exit(1);
	}

	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &prior_info->q);
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &prior_info->p);
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n\n", &c);
	fgets(tag, 100, fp);
	fscanf(fp, "%f\n", &prior_info->sigma);

	prior_info->c = hu2miu(c, MIU_AIR, MIU_WATER);

	fprintf(stdout, "\nPRIOR PARAMETERS:\n");
	fprintf(stdout, "q = %f\n", prior_info->q);
	fprintf(stdout, "p = %f\n", prior_info->p);
	fprintf(stdout, "c = %f mm-1 (%f HU)\n", prior_info->c, c);
	fprintf(stdout, "sigma = %f\n", prior_info->sigma);

	fclose(fp);
}
