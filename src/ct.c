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
#include "mpi.h"

void error(char *name)
{
	fprintf(stdout, "usage: %s [executable file] [forward model geometry directory] [Number of Focal Spots] [Recon parameters] [prior parameters ] [consensus equilibrium parameters]  [reconstruction output] [Number of iterations] [Dual Energy Flag] [Debug Mode Flag] [Number of X-ray Sources]\n", name);
}

int main(int argc, char *argv[])
{
	int numprocs, myid;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	int Nit;
	int num_focal_spots;
	int DE_mode=0;
	int debug_mode=0;
	char **recon_mask;
	struct Image image;
	struct Image V;
	struct Image consensus_X;
	struct Image TildeV;
	struct Image VPrevious;
	struct Image Vmean;
	struct Sinogram sinogram;
	struct PriorInfo prior_info;
	struct CEInfo ce_info;

	/* check arguments */
	if (argc != 11)
	{
		error(argv[0]);
		exit(1);
	}



	sinogram.geom_info.num_focal_spots = atoi(argv[2]);
	sinogram.geom_info.num_sources = atoi(argv[10]);
	
	if((numprocs%(sinogram.geom_info.num_focal_spots * sinogram.geom_info.num_sources))!=0){
		fprintf(stdout,"numprocs must be a multiple of the number of focal spots times the number of sources \n");
		exit(1);	
	}


	/* read sinogram, memory allocation for data array included */
	readAll_GeomDirectory(argv[1],  myid,numprocs, &(sinogram.geom_info));
	printGeomInfo(&(sinogram.geom_info));
	

	fprintf(stdout,"Reading sinogram...\n");

	fillSinogramData(&sinogram,numprocs,myid);
	/* reads/fills all available sinogram data (counts,dosage,etc.) */

	/* read initialized image, memory allocation for data array included */

	readImgInfo(argv[3], &(image.img_info));
	readImgInfo(argv[3],&(V.img_info));
        readImgInfo(argv[3],&(consensus_X.img_info));
        readImgInfo(argv[3],&(TildeV.img_info));
        readImgInfo(argv[3],&(VPrevious.img_info));
	readImgInfo(argv[3],&(Vmean.img_info));


	char fname_consensus[300];
    	sprintf(fname_consensus,"%smyid0.vjk",image.img_info.imgFile);		



	fprintf(stdout,"fname_consensus is %s\n",fname_consensus);


	/* fill in intermediate variables */
	
	fillGeomInfo(&(sinogram.geom_info));
	fillImgInfo(&(image.img_info));
	fillImgInfo(&(V.img_info));
	fillImgInfo(&(consensus_X.img_info));
	fillImgInfo(&(TildeV.img_info));
	fillImgInfo(&(VPrevious.img_info));
	fillImgInfo(&(Vmean.img_info));


	/* create and compute reconstruction mask */
	
	createReconMask(&recon_mask, &(image.img_info));
	if (strcmp(image.img_info.maskFile, "NA") == 0)
		compReconMask(recon_mask, &(image.img_info));
	else
	{
		printf("Reading recon mask...\n");
		readReconMask(recon_mask, &(image.img_info));
	}
	


	if (strcmp(image.img_info.imgFile, "NA") == 0)
	{
		initImage(&(sinogram.geom_info), &(image), 0.0192);
                initImage(&(sinogram.geom_info), &V, 0.0192);
                initImage(&(sinogram.geom_info), &consensus_X, 0.0192);
                initImage(&(sinogram.geom_info), &TildeV, 0.0192);
                initImage(&(sinogram.geom_info), &VPrevious, 0.0192);
	        initImage(&(sinogram.geom_info), &Vmean, 0.0192);


		for (int i =0; i < image.img_info.Nx; i++){
			for(int j=0;j<image.img_info.Ny;j++){
				if(!recon_mask[i][j]){
					for(int k=0;k<image.img_info.Nz;k++){
						image.img[i*image.img_info.Ny *image.img_info.Nz + j*image.img_info.Nz+k]=0;
						TildeV.img[i*image.img_info.Ny *image.img_info.Nz + j*image.img_info.Nz+k]=0;
						V.img[i*image.img_info.Ny *image.img_info.Nz + j*image.img_info.Nz+k]=0;
						consensus_X.img[i*image.img_info.Ny *image.img_info.Nz + j*image.img_info.Nz+k]=0;
						VPrevious.img[i*image.img_info.Ny *image.img_info.Nz + j*image.img_info.Nz+k]=0;
						Vmean.img[i*image.img_info.Ny *image.img_info.Nz + j*image.img_info.Nz+k]=0;


					}
				}	
			}
		}
	
	}
	else
	{
		readImage(fname_consensus, &image);
                readImage(fname_consensus, &V);
                readImage(fname_consensus, &consensus_X);
                readImage(fname_consensus, &TildeV);
                readImage(fname_consensus, &VPrevious);
		readImage(fname_consensus, &Vmean);
	}
	printImgInfo(&(image.img_info));
		

	/* read prior parameters */
	
	readPrior(argv[4], &prior_info);
	
	readCE(argv[5],&ce_info);


	/*
	int i=0;	
	float penalizer =0.0;
	for(i=0;i<(image.img_info.Nx * image.img_info.Ny * image.img_info.Nz);i++)
    		penalizer += (image.img[i]-TildeV.img[i])*(image.img[i]-TildeV.img[i])/(2*ce_info.SigmaLambda*ce_info.SigmaLambda);
   
	fprintf(stdout,"initial: penalizer is %f \n",penalizer); 
	*/






	/* read number of iterations */
	
	sscanf(argv[7], "%d", &Nit);
	fprintf(stdout, "\n# of iterations = %d\n", Nit);
	
	sscanf(argv[8], "%d", &DE_mode);
	fprintf(stdout,"\n DE mode %d \n",DE_mode);

	//sscanf(argv[9], "%d", &debug_mode);
	//fprintf(stdout,"\n debug mode %d \n",debug_mode);


	int DE_bin= 0;
	int DE_id=0;
	MPI_Comm DE_comm;
	int DE_numprocs=0;
	
	if(DE_mode==4){
        // single focal spot
		DE_bin = myid/3; // 3 focal spots in a group
		MPI_Comm_split(MPI_COMM_WORLD, DE_bin, myid, &DE_comm); 
		MPI_Comm_size(DE_comm, &DE_numprocs);
		MPI_Comm_rank(DE_comm, &DE_id);
	
	}
	else if(DE_mode==3){
        // single focal spot
		DE_bin = myid; // independent recon
		MPI_Comm_split(MPI_COMM_WORLD, DE_bin, myid, &DE_comm); 
		MPI_Comm_size(DE_comm, &DE_numprocs);
		MPI_Comm_rank(DE_comm, &DE_id);
	
	}

	else if(DE_mode==2){
        // single focal spot
		DE_bin = myid % 2; // two sources with focal spot position 1 are in the same group
		MPI_Comm_split(MPI_COMM_WORLD, DE_bin, myid, &DE_comm); 
		MPI_Comm_size(DE_comm, &DE_numprocs);
		MPI_Comm_rank(DE_comm, &DE_id);
	
	}

	else if(DE_mode==1){
        // dual energy
		DE_bin = myid / 2; // focal spot 1 and 2 are in the same group
		MPI_Comm_split(MPI_COMM_WORLD, DE_bin, myid, &DE_comm); 
		MPI_Comm_size(DE_comm, &DE_numprocs);
		MPI_Comm_rank(DE_comm, &DE_id);
	
	}
	else{
	// normal recon DE_mode=0
	//
		DE_bin = myid / numprocs; 
		MPI_Comm_split(MPI_COMM_WORLD, DE_bin, myid, &DE_comm); 
		MPI_Comm_size(DE_comm, &DE_numprocs);
		MPI_Comm_rank(DE_comm, &DE_id);

	}


	/* ICD algorithm and use q-GGMRF prior, will write out the reconstructed image */
        MPI_Barrier(MPI_COMM_WORLD);
	
	ICDReconstruct(&image, &V, &consensus_X, &TildeV, &VPrevious, &Vmean, &sinogram, recon_mask, &prior_info, &ce_info, argv[6], Nit, myid,numprocs,DE_id,DE_numprocs,&DE_comm, DE_mode, debug_mode,argv[3]);

	MPI_Barrier(MPI_COMM_WORLD);	



	/* free memory */
	
	freeSinogram(&sinogram);

	fprintf(stdout,"after freeing sinogram \n");
	fflush(stdout);

	freeImage(&image);

	fprintf(stdout,"after freeing image \n");
	fflush(stdout);



	freeImage(&V);

	fprintf(stdout,"after freeing V \n");
	fflush(stdout);



        freeImage(&consensus_X);

	fprintf(stdout,"after freeing consensus_X \n");
	fflush(stdout);


        freeImage(&TildeV);

	fprintf(stdout,"after freeing TildeV\n");
	fflush(stdout);


        freeImage(&VPrevious);

	fprintf(stdout,"after freeing VPrevious \n");
	fflush(stdout);



	freeImage(&Vmean);

	fprintf(stdout,"after freeing Vmean \n");
	fflush(stdout);



	freeReconMask(recon_mask);


	fprintf(stdout,"after freeing recon mask \n");
	fflush(stdout);


	
	return 0;
}
