/* version 6.0, P.Jin 07/13/2012 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "data.h"
#include "io.h"
#include "proj.h"
#include "randlib.h"
#include "prepro.h"
#include "allocate.h"

void error(char *name)
{
	fprintf(stdout, "usage: %s [image param] [geom param]\n", name);
}

// Update to use the correct signature for fowardProject, and call it with Y=0.  The output will be -sinogram.
// We can use this to project a phantom to give a sinogram.
int main(int argc, char *argv[])
{
	struct timeval start, end;

	struct Image image;
	struct Sinogram sinogram;
	char **recon_mask;
	unsigned short *AX_mask=NULL;
	int jx,jy;

	/* check arguments */
	if (argc != 3)
	{
		error(argv[0]);
		exit(1);
	}

	/* start timer */
	gettimeofday(&start, NULL);

	/* read image info */
	/* read image data, memory allocation for data array included */
	readImgInfo(argv[2], &(image.img_info));
	printImgInfo(&(image.img_info));
	readImage_short(image.img_info.imgFile, &image);

	/* read geomtry info */
	readGeomInfo(argv[1], &(sinogram.geom_info));
	printGeomInfo(&(sinogram.geom_info));

	/* fill in intermediate variables */
	fillImgInfo(&(image.img_info));
	fillGeomInfo(&(sinogram.geom_info));
        /*checkInfo(&(sinogram.geom_info),&(image.img_info));*/  /* sjk */

	/* allocate memory for sinogram */
	createSinogram(&sinogram);

	/* forward projection */
	createReconMask(&recon_mask, &(image.img_info));
	for(jx=0; jx<image.img_info.Nx; jx++)
	for(jy=0; jy<image.img_info.Ny; jy++)
		recon_mask[jx][jy]=1;
	/*forwardProject(sinogram.sino, image.img, recon_mask, &(sinogram.geom_info), &(image.img_info));*/
	forwardProject(sinogram.sino, image.img, AX_mask, recon_mask, &(sinogram.geom_info), &(image.img_info));

	/* add noise */
	if (sinogram.geom_info.lambda0 > 0.0)
	{
		srandom2(1);
		addNoise(sinogram.sino, sinogram.geom_info.lambda0, sinogram.geom_info.Nv*sinogram.geom_info.Nc*sinogram.geom_info.Nr);
	}

	/* write sinogram */
	/*writeSinogram(sinogram.geom_info.sinoFile, &sinogram);*/
        writeSinogram_float(sinogram.geom_info.sinoFile, sinogram.sino,
                     sinogram.geom_info.Nr, sinogram.geom_info.Nc, sinogram.geom_info.Nv);

	/* record running time */
	gettimeofday(&end, NULL);
	fprintf(stdout, "\nProcessing time = %ld secs\n", (end.tv_sec - start.tv_sec));

	/* free memory */
	freeImage(&image);
	freeSinogram(&sinogram);

	return 0;
}
