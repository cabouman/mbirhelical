#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "data.h"
#include "proj.h"
#include "icd.h"
#include "prepro.h"
#include "io.h"

void error(char *name)
{
	fprintf(stdout, "usage: %s [geom param] [image param] [prior param] [recon image] [iterations]\n", name);
}

int main(int argc, char *argv[])
{
	struct timeval start, end;
	int Nit;
	char **recon_mask;
	struct Image image;
	struct Sinogram sinogram;
	struct PriorInfo prior_info;


	/* check arguments */

	
	if (argc != 6)
	{
		error(argv[0]);
		exit(1);
	}
	

	/* start timer */
	
	gettimeofday(&start, NULL);
	
	/* read sinogram, memory allocation for data array included */
	
	readGeomInfo(argv[1], &(sinogram.geom_info));
	printGeomInfo(&(sinogram.geom_info));
	

	
	fprintf(stdout,"Reading sinogram...\n");

	fillSinogramData(&sinogram);
	   /* reads/fills all available sinogram data (counts,dosage,etc.) */


	/* read initialized image, memory allocation for data array included */
	
	readImgInfo(argv[2], &(image.img_info));
	if (strcmp(image.img_info.imgFile, "NA") == 0)
	{
		initImage(&(sinogram.geom_info), &(image), 500.0);
	}
	else
	{
		readImage_short(image.img_info.imgFile, &image);
	}
	printImgInfo(&(image.img_info));
	

	gettimeofday(&end, NULL);
	fprintf(stdout, "\nRead time = %ld secs\n", (end.tv_sec - start.tv_sec));
	


	/* read prior parameters */
	
	readPrior(argv[3], &prior_info);
	


	/* fill in intermediate variables */
	
	fillGeomInfo(&(sinogram.geom_info));
	fillImgInfo(&(image.img_info));
	


	/* create and compute reconstruction mask */
	
	createReconMask(&recon_mask, &(image.img_info));
	if (strcmp(image.img_info.maskFile, "NA") == 0)
		compReconMask(recon_mask, &(image.img_info));
	else
	{
		printf("Reading recon mask...\n");
		readReconMask(recon_mask, &(image.img_info));
	}
	


	/* read number of iterations */
	
	sscanf(argv[5], "%d", &Nit);
	fprintf(stdout, "\n# of iterations = %d\n", Nit);
	


	/* ICD algorithm and use q-GGMRF prior, will write out the reconstructed image */
	
	ICDReconstruct(&image, &sinogram, recon_mask, &prior_info, argv[4], Nit);
	


	/* record running time */
	
	gettimeofday(&end, NULL);
	fprintf(stdout, "\nProcessing time = %ld secs\n", (end.tv_sec - start.tv_sec));




	/* free memory */
	
	freeSinogram(&sinogram);
	freeImage(&image);
	freeReconMask(recon_mask);
	
	return 0;
}
