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
#include "mpi.h"
#include "tensorflow/c/c_api.h"

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
	float scalar,Dxy2,Dz2;

	Dxy2=img_info->Del_xy * img_info->Del_xy;
	Dz2= img_info->Del_z * img_info->Del_z;

	j=0;
	scalar=  1/sqrt(Dxy2);
	for(jx=-1; jx<=1; jx++)   /* NOTE: order is critical here--has to match fillNeighbors convention */
	for(jy=-1; jy<=1; jy++)
	for(jz=-1; jz<=1; jz++)
	if((jx!=0)||(jy!=0)||(jz!=0)) 
	{
		nb_wt[j] = 1/sqrt((jx*jx*Dxy2)+(jy*jy*Dxy2)+(jz*jz*Dz2));
		j++;
	}
	for(j=0; j<26; j++)
		nb_wt[j] = nb_wt[j]/scalar;   //square root the prior weight-Xiao Wang

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



void updateError(ENTRY *e, struct ACol *col_xyz, float diff)
{
	int r;

	for (r = 0; r < col_xyz->n_index; r++)
	{
		#pragma omp atomic
		e[col_xyz->index[r]] -= (diff*col_xyz->val[r]);
	}
}


void TwoGMinusIOperator_PnP(
                                    struct Image *TildeV,                   /* Output (2G-I)V, where G = merge, denoise & stack operator here */  
                                    struct Image *V,                        /* Input V */   
                                    struct Image *consensus_X,                         /* Output Z - merge V and denoise. If denoising is an optimization problem, Z acts as input (initial state) as well */
				    struct Image *Vmean,
        			    struct PriorInfo *prior_info, char **reconMask, float SigmaLambda, int *order,int myid,int NUM_PROCS,int it,int DE_id,int DE_numprocs,MPI_Comm *DE_comm,int debug_mode)
{

	/* In this function, if PnP=0, the consensus equilibrium method computes the mean of V, and then compute TildeV as 2 times mean - V */
	/* if PnP=1, the plug and play framework performs proximal map prior computation on a single node, and then broadcast the denoised result, consensus_X, to all nodes.
	 * Then compute TildeV as 2 times consensus_X -V */

    /* Merge V and Denoise. Store result in consensus_X */
    MergeAndDenoise(V, Vmean, consensus_X, prior_info, reconMask, SigmaLambda, order,myid,NUM_PROCS,it,DE_id,DE_numprocs,DE_comm,debug_mode);             

    /* (2G-I)V = 2*stacked(Z)-V, since G= merge, denoise and stack operator here */
    StackAndReflect(TildeV, V, consensus_X, Vmean);    
}


void StackAndReflect(
                            struct Image *TildeV,              /* Output = 2*Stacked(Z)-V , i.e. Reflect V about Z.  (Z is computed before hand) */   
                            struct Image *V,                       /* Input V */   
                            struct Image *consensus_X, struct Image *Vmean                        /* Input Z */)              
{

    int i;


    /* Stack Z*/
	if(PnP_mode){

    		for(i=0;i<V->img_info.Nx * V->img_info.Ny * V->img_info.Nz;i++)
    		{
        		TildeV->img[i] = 2*consensus_X->img[i]-V->img[i];
    		}

   	}
	else{
     		for(i=0;i<V->img_info.Nx * V->img_info.Ny * V->img_info.Nz;i++)
    		{
        		TildeV->img[i] = 2*Vmean->img[i]-V->img[i];
    		}
 	
   	}
/*
	
    ENTRY *avg_TildeV = (ENTRY *)get_spc(V->img_info.Nx * V->img_info.Ny * V->img_info.Nz, sizeof(ENTRY));
    int j=0;
    int Nx = V->img_info.Nx;
    int Ny = V->img_info.Ny;
    int Nz = V->img_info.Nz;
    for (i =0; i < Nx; i++){
	for(j=0;j< Ny;j++){
		MPI_Allreduce(&(TildeV->img[i*Ny*Nz+j*Nz]),&(avg_TildeV[i*Ny*Nz +j*Nz]), Nz, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	}
    }
    for (i =0;i<Nx*Ny*Nz;i++){
	avg_TildeV[i] = avg_TildeV[i]/4;
    }

    for (i=0;i<Nx*Ny*Nz;i++){
	if(avg_TildeV[i] != (2*consensus_X->img[i]-Vmean->img[i]))
		fprintf(stdout,"i %d avg_TildeV->img[i] %f consensus_X->img[i] %f Vmean %f \n",i,avg_TildeV[i],consensus_X->img[i],Vmean->img[i]);
    }
    free((void *)avg_TildeV);
*/

}


void MergeAndDenoise(   struct Image *V, struct Image *Vmean,                       /* Input V */   
                        struct Image *consensus_X,                       /* Output Z - merge V and denoise. If denoising is an optimization problem, Z acts as input (initial state) as well */
	        	struct PriorInfo *prior_info, char **reconMask, float SigmaLambda, int *order, int myid, int NUM_PROCS,int it,int DE_id,
int DE_numprocs,MPI_Comm *DE_comm, int debug_mode)              
{
    /* Allocate memory */

    /* No need for merging because only 1 node */
    
    /*Denoise Vmean */
    int Nx=V->img_info.Nx;
    int Ny=V->img_info.Ny;
    int Nz=V->img_info.Nz;

    
    int i=0;
    int j=0;
    for (i =0; i < Nx; i++){
	MPI_Allreduce(&(V->img[i*Ny*Nz]),&(Vmean->img[i*Ny*Nz]), Ny*Nz, MPI_FLOAT, MPI_SUM, *DE_comm);
    }
    for (i =0;i<Nx*Ny*Nz;i++){
	Vmean->img[i] = Vmean->img[i]/DE_numprocs;
    }



    if(PnP_mode){
    	if(!debug_mode){
    		if(DE_id==0){
    			SolveProximalMap_Prior(consensus_X, reconMask, Vmean, SigmaLambda, prior_info, order, it,myid);
    		
		}
	    	for (i=0;i<Nx;i++){
			MPI_Bcast(&consensus_X->img[i*Ny*Nz],Ny*Nz, MPI_FLOAT,0,*DE_comm);
    		}
	

    	}	
    	else{
		fprintf(stdout,"warning: debug mode is turned on ! \n");
		fflush(stdout);
		SolveProximalMap_Prior(consensus_X, reconMask, V, SigmaLambda, prior_info, order, it,myid);

    	}
    }


}


void NoOpDeallocator(void* data, size_t a, void* b) {}

void SolveProximalMap_Prior(struct Image *Image, 
                            char  ** reconMask, 
                            struct Image *ProximalMapInput, 
                            float  SigmaLambda, struct PriorInfo *prior_info, int *order,int it,int myid)
{
	#pragma omp parallel
	{  
	const int Nxy =Image->img_info.Nx * Image->img_info.Ny; /* image size */

  	TF_Graph* Graph = TF_NewGraph();
	TF_Status* Status = TF_NewStatus();

  	TF_SessionOptions* SessionOpts = TF_NewSessionOptions();
  	TF_Buffer* RunOpts = NULL;

  	// Get path to model directory from input
  	const char* saved_model_dir = prior_info->DL_File;
  	//fprintf(stdout,"Model: %s\n", saved_model_dir);
  	
	// model serve tag
  	const char* tags = "serve";
  	int ntags = 1;

  	TF_Session* Session = TF_LoadSessionFromSavedModel(SessionOpts, RunOpts, saved_model_dir, &tags, ntags, Graph, NULL, Status);

  	if(TF_GetCode(Status) != TF_OK)
    		fprintf(stderr,"%s\n",TF_Message(Status));

  	//****** Get input tensor
  	int NumInputs = 1;
  	TF_Output* Input = (TF_Output*)malloc(sizeof(TF_Output) * NumInputs);

  	TF_Output t0 = {TF_GraphOperationByName(Graph, "serving_default_input_1"), 0};
  	if(t0.oper == NULL)
    		fprintf(stderr,"ERROR: Failed TF_GraphOperationByName serving_default_input_1\n");

  	Input[0] = t0;

  	//********* Get Output tensor
  	int NumOutputs = 1;
  	TF_Output* Output = (TF_Output*)malloc(sizeof(TF_Output) * NumOutputs);

  	TF_Output t2 = {TF_GraphOperationByName(Graph, "StatefulPartitionedCall"), 0};
  	if(t2.oper == NULL)
   		fprintf(stderr,"ERROR: Failed TF_GraphOperationByName StatefulPartitionedCall\n");

  	Output[0] = t2;

  	//********* Allocate data for inputs & outputs
  	TF_Tensor** InputValues = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*NumInputs);
  	TF_Tensor** OutputValues = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*NumOutputs);

  	// set the dimensions of the images here
  	int ndims      = 4;  // TF requires an additional dim for batch size: 3+1 (nx,ny,nslice,nbatch)
  	int model_ipsize_x   = Image->img_info.Nx;
  	int model_ipsize_y   = Image->img_info.Ny;
  	int model_ipsize_z   = 5;
  	int batch_size = 1;  

  	// allocate TF arrays
  	int64_t dims[] = {batch_size,model_ipsize_x,model_ipsize_y,model_ipsize_z};

  	//fprintf(stdout,"batch_size %d,model_ipsize_x %d, model_ipsize_y %d, model_ipsize_z %d \n",batch_size,model_ipsize_x,model_ipsize_y,model_ipsize_z);

  	ENTRY   data[batch_size][model_ipsize_x][model_ipsize_y][model_ipsize_z];
  	int ndata = sizeof(ENTRY)*batch_size*model_ipsize_x*model_ipsize_y*model_ipsize_z; // number of bytes not number of elements

  // curate input data for given image slice
  // Image->img[nx][ny][slice]
  // NOTE: for batch_size = 1
  //
  //
  //
  
 


  
  	#pragma omp for
 	for (int k=0;k<Image->img_info.Nz;k++){     // this is the slice we are de-noising


  		for (int i=0; i<model_ipsize_x; i++) {
  			for(int j=0; j<model_ipsize_y; j++) {
  				for(int kk=0; kk<model_ipsize_z; kk++) {
          				if((k==0 || k==1) && (kk==0 || kk==1)){
            					data[0][i][j][kk] = ProximalMapInput->img[i*Image->img_info.Ny*Image->img_info.Nz + j*Image->img_info.Nz  + 0];
          				}
			
          				else if((k==(Image->img_info.Nz-1) || k==(Image->img_info.Nz-2)) && (kk==(model_ipsize_z-2) || kk==(model_ipsize_z-1))){
            					data[0][i][j][kk] = ProximalMapInput->img[i*Image->img_info.Ny*Image->img_info.Nz + j*Image->img_info.Nz  + Image->img_info.Nz-1];
          				}
         
          				else
            					data[0][i][j][kk] = ProximalMapInput->img[i*Image->img_info.Ny*Image->img_info.Nz + j*Image->img_info.Nz  + k+kk-2];
        			}
      			}
    		}


    		// no deallocator needed because data was put on a stack 
    		TF_Tensor* int_tensor = TF_NewTensor(TF_FLOAT, dims, ndims, data, ndata, &NoOpDeallocator, 0);

	    	if (int_tensor == NULL)
			fprintf(stderr,"ERROR: Failed TF_NewTensor\n");

	    	InputValues[0] = int_tensor;


	    	// ================================
	    	// Run forward pass of model
	    	// ================================
	    	TF_SessionRun(Session, NULL, Input, InputValues, NumInputs, Output, OutputValues, NumOutputs, NULL, 0,NULL , Status);

	    	if(TF_GetCode(Status) != TF_OK)
	    	{
	      		fprintf(stderr,"%s\n",TF_Message(Status));
	    	}


    		// write oiriginal Vmean image to the directory before subtracting noise from it
    		if(k==100){
	    		ENTRY* temp = (ENTRY *)  get_spc((Image->img_info.Ny)*(Image->img_info.Nx), sizeof(ENTRY));
      			for(int i=0; i<model_ipsize_x; i++) {
        			for (int j=0; j<model_ipsize_y; j++) {
          				temp[i*Image->img_info.Ny+j]=ProximalMapInput->img[i*Image->img_info.Ny*Image->img_info.Nz + j*Image->img_info.Nz + k]; 
        			}
      			}
      			if(myid==0){
	      			char errorFname[200];
	      			sprintf(errorFname,"/gpfs/alpine/gen006/proj-shared/xf9/recon/dcm134/input_%d",it);		    	  
        			writeSinogram_float(errorFname, (ENTRY *) temp, Image->img_info.Nx, Image->img_info.Ny, 1);
      			}
      			free(temp);
    		}



	    	// subtract noise (model output) from original image and then clip it for non-negativity
    		ENTRY* outvalues = (ENTRY*) TF_TensorData(OutputValues[0]);  // which is the extracted noise
	    	int counter = 0;
	    	for(int i=0; i<model_ipsize_x; i++) {
	      		for (int j=0; j<model_ipsize_y; j++) {
				if (reconMask[i][j]){

		  			ENTRY pixel = ProximalMapInput->img[i*Image->img_info.Ny*Image->img_info.Nz + j*Image->img_info.Nz + k]-outvalues[counter]; 
		  			//clip 
		  			if(positive_constraint ==1){
		    				Image->img[i*Image->img_info.Ny*Image->img_info.Nz + j*Image->img_info.Nz + k] = ((pixel < 0.0) ? 0.0 : pixel); 
		  			}
		  			else{
		    				Image->img[i*Image->img_info.Ny*Image->img_info.Nz + j*Image->img_info.Nz + k] = pixel;
		  			}
				}
				counter++;
	      		}
	    	}



    		// write ouput image
    		if(k==100){
	    		ENTRY* temp = (ENTRY *)  get_spc((Image->img_info.Ny)*(Image->img_info.Nx), sizeof(ENTRY));
      			for(int i=0; i<model_ipsize_x; i++) {
        			for (int j=0; j<model_ipsize_y; j++) {
          				temp[i*Image->img_info.Ny+j]=Image->img[i*Image->img_info.Ny*Image->img_info.Nz + j*Image->img_info.Nz + k]; 
        			}
      			}
      			if(myid==0){
	      			char errorFname[200];
	      			sprintf(errorFname,"/gpfs/alpine/gen006/proj-shared/xf9/recon/dcm134/output_%d",it);		    	  
        			writeSinogram_float(errorFname, (ENTRY *) temp, Image->img_info.Nx, Image->img_info.Ny, 1);
      			}
      			free(temp);
    		}


	}
		



	// Free memory
  	TF_DeleteGraph(Graph);
	TF_DeleteSession(Session, Status);
	TF_DeleteSessionOptions(SessionOpts);
	TF_DeleteStatus(Status);
	free(Input);
	free(Output);
	free(InputValues);
	free(OutputValues);
	
	}

}





void paraICD_Prior(struct ImgInfo *img_info,struct PriorInfo *prior_info,char **recon_mask, ENTRY *X,struct Image *ProximalMapInput, int *order, float SigmaLambda,int it,struct Image *Image)
{
	int tid = omp_get_thread_num();
	int offset;  /* sjk */
	int j, jx, jy, jz, jzmax, Nxy, Nyz;
	float x, y, z, pixel, diff;
	struct ACol col_xyz;
	struct ViewXYZInfo view_xyz_info;
	struct ICDInfo icd_info;

	int k, l, Nsubit, jjx, jjy;
	float gamma;
	ENTRY **VSC;
	/* sjk */
	int *voxel_list;
	float *key;
	char zero_skip_flag;  
	int n;  



	VSC = (ENTRY **)get_img(img_info->Ny, img_info->Nx, sizeof(ENTRY));
	/* sjk */
  voxel_list = (int *)get_spc(img_info->Nx*img_info->Ny,sizeof(int));
  key = (float *)get_spc(img_info->Nx*img_info->Ny,sizeof(float));

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

		VSC[jx][jy] = 0.0;

		if (recon_mask[jx][jy])
		{
			x = img_info->x0 + jx*img_info->Del_xy;
			y = img_info->y0 + jy*img_info->Del_xy;
			offset= jx*Nyz + jy*img_info->Nz;  /* sjk */
			for (jz = tid*img_info->Nz/omp_get_num_threads(); jz < jzmax; jz++)
			{

				/* sjk: moved this here to facilitate zero-skip check */
				fillNeighbors(icd_info.neighbors,jx,jy,jz,X,img_info); 

				/*** sjk: zero skip check ***/

				if( X[offset+jz] > 0  || positive_constraint ==0)
				  zero_skip_flag=0;  /* don't skip */
				else
				{
				  zero_skip_flag=1;  /* maybe skip.. check neighbors */
				 
				  if( (icd_info.neighbors[1] > 0) || (icd_info.neighbors[4] > 0) ||(icd_info.neighbors[7] > 0) ||(icd_info.neighbors[10] > 0) ||(icd_info.neighbors[15] > 0) ||(icd_info.neighbors[18] > 0) ||(icd_info.neighbors[21] > 0) ||(icd_info.neighbors[24] > 0) )
				  {
				    zero_skip_flag=0;
				  }
				}



				if(zero_skip_flag==0)
				{

				z = img_info->z0 + jz*img_info->Del_z;


				pixel = ICDStep_PriorOnly(&icd_info, prior_info, jx, jy, jz, X, ProximalMapInput, img_info, SigmaLambda);

					/* clip */
				if(positive_constraint ==1){
					X[offset+jz] = ((pixel < 0.0) ? 0.0 : pixel);  /* sjk */
				}
				else{
					X[offset+jz] = pixel; 
				}
				/* no clipping for now trying*/
				//X[offset+jz] = pixel;  /* sjk */


				/* sjk: try without positivity constraint */
				if (recon_mask[jx][jy]==2)
					X[offset+jz] = pixel; 

				diff = X[offset+jz]-icd_info.v;  /* sjk */
	
				VSC[jx][jy] += fabsf(diff);
				}  /* END zero_skip_flag */
			}
		}
		
	}


	/*
	int myid;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	char name[200];
	char suffix[100];
	char fname[200];
	strcpy(fname,"");
   	strcpy(fname,  "/scratch/rice/w/wang1698/lung_cancer/CE_method/recon");


        strcpy(name, "");
        sprintf(suffix, "myid%d_iteration_%d_Prior_After_Random.vjk",myid,it);
        strcat(name, fname);
        strcat(name, suffix);
	#pragma omp single
	{
        	writeImage(name, Image);

	}
	*/









	/* NH-ICD */
	Nsubit = 1;
	for (k = 0; k < Nsubit; k++)
	{


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

				VSC[jjx][jjy] = 0.0;
				x = img_info->x0 + jjx*img_info->Del_xy;
				y = img_info->y0 + jjy*img_info->Del_xy;

				offset= jjx*Nyz + jjy*img_info->Nz;  /* sjk */
				for (jz = tid*img_info->Nz/omp_get_num_threads(); jz < jzmax; jz++)
				{
					/* sjk: moved this here to facilitate zero-skip check */
					fillNeighbors(icd_info.neighbors,jjx,jjy,jz,X,img_info); 

				        /*** sjk: zero skip check ***/


					if( X[offset+jz] > 0 || positive_constraint ==0)
					  zero_skip_flag=0;  /* don't skip */
					else
					{
						zero_skip_flag=1;  /* maybe skip.. check neighbors */
				 
				  		if( (icd_info.neighbors[1] > 0) || (icd_info.neighbors[4] > 0) ||(icd_info.neighbors[7] > 0) ||(icd_info.neighbors[10] > 0) ||(icd_info.neighbors[15] > 0) ||(icd_info.neighbors[18] > 0) ||(icd_info.neighbors[21] > 0) ||(icd_info.neighbors[24] > 0) )
				  		{
				    			zero_skip_flag=0;
				  		}
					}

					if(zero_skip_flag==0)
					{

						z = img_info->z0 + jz*img_info->Del_z;

						/* calculate A column on the fly*/
						/* ICD update */

						pixel = ICDStep_PriorOnly(&icd_info, prior_info, jjx, jjy, jz, X, ProximalMapInput, img_info, SigmaLambda);

						/* clip */
						if (positive_constraint ==1){	
							X[offset+jz] = ((pixel < 0.0) ? 0.0 : pixel);  /* sjk */
						}
						else{
							X[offset+jz] = pixel; 

						}


						/* no clipping for now try */
						//X[offset+jz] = pixel; 



						/* sjk: try without positivity constraint */
						if (recon_mask[jjx][jjy]==2)
							X[offset+jz] = pixel; 
						
						diff = X[offset+jz]-icd_info.v;  /* sjk */
						VSC[jjx][jjy] += fabsf(diff);


					}  /* END zero_skip_flag */
				}

			}
		}
	}




	/*
 	strcpy(name, "");
        sprintf(suffix, "myid%d_iteration_%d_Prior_After_NHICD.vjk",myid,it);
        strcat(name, fname);
        strcat(name, suffix);
	#pragma omp single
	{
        	writeImage(name, Image);
	}
	*/

	/* sjk */
	free((void *)voxel_list);
	free((void *)key);
	free_img((void **)VSC);


}



float ICDStep_PriorOnly(
	struct ICDInfo *icd_info,
	struct PriorInfo *prior_info,
	int jx,
	int jy,
	int jz,
	ENTRY *X,			/* image, attenuation coefficients */
	struct Image *ProximalMapInput,
	struct ImgInfo *img_info,	/* image parameters */
	float SigmaLambda)
{
	int i, r;
	float d;
	/*float dsum=0.0;*/
	float proximalMapValue =0;
	/* store old pixel value */
	icd_info->v = (float)X[jx*(img_info->Ny*img_info->Nz)+jy*img_info->Nz+jz];
	proximalMapValue = (float)ProximalMapInput->img[jx*(img_info->Ny*img_info->Nz)+jy*img_info->Nz+jz];


	/* copmute theta1 and theta2 */
	icd_info->th1 = 0.0;
	icd_info->th2 = 0.0;
	
	if (prior_info->est == QGGMRF)
	{
    		float delta, SurrogateCoeff;
    		int j;
    		for (j = 0; j < 26; j++)
    		{
        		delta = icd_info->v - icd_info->neighbors[j];
        		SurrogateCoeff = QGGMRF_SurrogateCoeff(delta,prior_info);
        
        		icd_info->th1 += (icd_info->nb_wt[j]*SurrogateCoeff * delta);
        		icd_info->th2 += icd_info->nb_wt[j]*SurrogateCoeff;
    		}
		icd_info->th1 += (icd_info->v - proximalMapValue)/(SigmaLambda *SigmaLambda);
		icd_info->th2 += 1/(SigmaLambda *SigmaLambda);

 		return icd_info->v - icd_info->th1/icd_info->th2;
	}

	else
	{
        	printf("Error! Other prior modes are not implemented yet! \n");
		exit(0);
	}
}










float ProximalMapCostFunction3D_Prior(struct Image *Image, 
                                      struct PriorInfo *prior_info,int myid)
{
    int imglen, Nxy;
    int jx, jy, jz;
    int Nyz;
    float nlogprior;
    float nb_wt[26];
    int nb_list[13] = {2,5,8,11,13,15,16,18,19,21,22,24,25};  /* only these used for unique cliques */


    imglen = Image->img_info.Nx * Image->img_info.Ny * Image->img_info.Nz;
    Nxy = Image->img_info.Nx * Image->img_info.Ny;    


    /* unified q-GMMRF, including quad, p-GGMRF */
    /* negative log prior term */

    compNeighborWeight(nb_wt,&Image->img_info);  /* sjk */
    nlogprior = 0.0;

    Nyz = Image->img_info.Ny*Image->img_info.Nz;

    #pragma omp parallel for collapse(3) reduction(+:nlogprior)
    for (jx = 0; jx < Image->img_info.Nx; jx++)
    {
    	for (jy = 0; jy < Image->img_info.Ny; jy++)
      {
              for (jz = 0; jz < Image->img_info.Nz; jz++)
              {
  
                int j = jx*Nyz + jy*Image->img_info.Nz + jz;
                float neighbors[26];
                fillNeighbors(neighbors,jx,jy,jz,Image->img,&Image->img_info);
                int i=0;
                if(PnP_mode==0){
                  for(i=0; i<13; i++)
                  {
                      float delta = Image->img[j]-neighbors[nb_list[i]];
                      nlogprior += nb_wt[nb_list[i]]*QGGMRF_Potential(delta,prior_info);
                  }
                }
                else{
                  for(i=0; i<13; i++)
                  {
                      nlogprior += nb_wt[nb_list[i]]*(Image->img[j]-neighbors[nb_list[i]])*(Image->img[j]-neighbors[nb_list[i]]);
                  }
               
                }
              }
      }
    }
     
    return nlogprior;
}   

float Consensus_Cost(struct Image *consensus_X,struct Image *X, struct Image *Vmean,int NUM_PROCS,int DE_id,int DE_numprocs,MPI_Comm *DE_comm)
{
    int Nx = consensus_X->img_info.Nx;
    int Ny = consensus_X->img_info.Ny;
    int Nz = consensus_X->img_info.Nz;
 
    int i=0;
    int j=0;
    int k=0;
    float squared_error=0;

    if(PnP_mode){ /* X average compares with consensus_X*/
    	for (i =0; i < Nx; i++){
		MPI_Allreduce(&(X->img[i*Ny*Nz]),&(Vmean->img[i*Ny*Nz]), Ny*Nz, MPI_FLOAT, MPI_SUM, *DE_comm);
    	}
   	#pragma omp parallel for 
    	for (i =0;i<Nx*Ny*Nz;i++){
		Vmean->img[i] = Vmean->img[i]/DE_numprocs;
    	}

    	#pragma omp parallel for collapse(3)  reduction(+:squared_error)
    	for (i=0; i<Nx; i++){
		for(j=0;j< Ny;j++){
			for(k=0;k<Nz;k++){
				squared_error+=(Vmean->img[i*Ny*Nz+j*Nz+k]-consensus_X->img[i*Ny*Nz+j*Nz+k])*(Vmean->img[i*Ny*Nz+j*Nz+k]-consensus_X->img[i*Ny*Nz+j*Nz+k]);
			}
		}
    	}	
    }
    else{ /* Vmean compares with Xmean */
     	for (i =0; i < Nx; i++){
		MPI_Allreduce(&(X->img[i*Ny*Nz]),&(consensus_X->img[i*Ny*Nz]), Ny*Nz, MPI_FLOAT, MPI_SUM, *DE_comm);
    	}
   	#pragma omp parallel for 
    	for (i =0;i<Nx*Ny*Nz;i++){
		consensus_X->img[i] = consensus_X->img[i]/DE_numprocs;
    	}

    	#pragma omp parallel for collapse(3)  reduction(+:squared_error)
    	for (i=0; i<Nx; i++){
		for(j=0;j< Ny;j++){
			for(k=0;k<Nz;k++){
				squared_error+=(Vmean->img[i*Ny*Nz+j*Nz+k]-consensus_X->img[i*Ny*Nz+j*Nz+k])*(Vmean->img[i*Ny*Nz+j*Nz+k]-consensus_X->img[i*Ny*Nz+j*Nz+k]);
			}
		}
    	}
   
    }
    return squared_error;
}


float ProximalMapCostFunction3D_Likelihood(
	ENTRY *X,
	ENTRY *e,
	ENTRY *D,
	ENTRY *TildeV,
	struct PriorInfo *prior_info,
	struct GeomInfo *geom_info,
	struct ImgInfo *img_info,
	float SigmaLambda,int myid)
{
	int i;
	int j, jx, jy, jz;
	int Nvcr;
	int imglen = img_info->Nx * img_info->Ny * img_info->Nz;
	float d;
	float nloglike, penalizer;
	/* sjk */
	float neighbors[26];

	/* negative log likelihood term */
	nloglike = 0.0;
	penalizer= 0.0;
	Nvcr = (geom_info->Nv)*(geom_info->Nc)*(geom_info->Nr);  
	for (i = 0; i < Nvcr; i++)
	{
		if (geom_info->lambda0 == 0.0)
			nloglike += (e[i])*(e[i]);
		else
		{
			d = D[i];   /* sjk */
			nloglike += (e[i])*d*(e[i]);
		}
	}

	nloglike /= 2.0;

	//for(i=0;i<imglen;i++)
    	//	penalizer += (X[i]-TildeV[i])*(X[i]-TildeV[i])/(2*SigmaLambda*SigmaLambda);
   
	//fprintf(stdout,"myid %d: nloglike %f penalizer %f total %f \n",myid,nloglike,penalizer,nloglike+penalizer);
	//fflush(stdout);
	return nloglike;
	//return nloglike + penalizer; 

}

float ICDStep_Likelihood(
	struct ICDInfo *icd_info,
	struct PriorInfo *prior_info,
	int jx,
	int jy,
	int jz,
	ENTRY *X,			/* image, attenuation coefficients */
	ENTRY *e,			/* current AX */
	ENTRY *D,			/* sjk: noise matrix */
	ENTRY *TildeV,
	float lambda0,
	struct ImgInfo *img_info,	/* image parameters */
	struct ACol *col_xyz,
	float SigmaLambda,
	ENTRY *consensus_X, 
	int it,int DE_numprocs)		/* A column for this (jx, jy, jz) pixel */
{
	int i, r;
	float d;
	/*float dsum=0.0;*/

	/* store old pixel value */
	icd_info->v = (float)X[jx*(img_info->Ny*img_info->Nz)+jy*img_info->Nz+jz];
	float TildeV_value = (float)TildeV[jx*(img_info->Ny*img_info->Nz)+jy*img_info->Nz+jz];
 	float consensus_value = (float)consensus_X[jx*(img_info->Ny*img_info->Nz)+jy*img_info->Nz+jz];

	/* copmute theta1 and theta2 */
	icd_info->th1 = 0.0;
	icd_info->th2 = 0.0;



	if (lambda0 == 0.0)
	{
		for (r = 0; r < col_xyz->n_index; r++)
		{
			icd_info->th1 -= col_xyz->val[r]*e[col_xyz->index[r]];
			icd_info->th2 += col_xyz->val[r]*col_xyz->val[r];
		}
	}
	else
	{
		for (r = 0; r < col_xyz->n_index; r++)
		{
			d = D[col_xyz->index[r]];  


			icd_info->th1 -= (col_xyz->val[r]*d*e[col_xyz->index[r]]);
			icd_info->th2 += (col_xyz->val[r]*d*col_xyz->val[r]);
		}
	}

	if(PnP_mode){	

	/* if PnP_mode ==1, the algorithm uses plug and play deep learning prior model.
	 * adding the first derivative and second derivative for the penalizing term to th1 and th2 without computing the prior*/

		icd_info->th1 += (icd_info->v-TildeV_value)/(SigmaLambda*SigmaLambda);
		icd_info->th2 += 1/(SigmaLambda *SigmaLambda);
	}
	else{    
	/* if PnP_mode ==0, the algorithm is in the Consensus Equilibrium mode.
	 * Continue and modify th1 and th2 based on the QGGMRF statistsical prior model
	 */
		float oldTheta1=icd_info->th1;
		float oldTheta2=icd_info->th2;

		QGGMRF3D_UpdateICDParams(icd_info,prior_info);
		icd_info->th1 = oldTheta1 + 1.0/DE_numprocs * (icd_info->th1 - oldTheta1);
		icd_info->th2 = oldTheta2 + 1.0/DE_numprocs * (icd_info->th2 - oldTheta2);

		/* adding the first derivative and second derivative for the penalizing term to th1 and th2 */

		icd_info->th1 += (icd_info->v-TildeV_value)/(SigmaLambda*SigmaLambda);
		icd_info->th2 += 1/(SigmaLambda *SigmaLambda);

	}

        return icd_info->v - (icd_info->th1/icd_info->th2) ;

}






/* ICD update with the QGGMRF prior model */
/* Prior and neighborhood specific */
void QGGMRF3D_UpdateICDParams(struct ICDInfo *icd_info, struct PriorInfo *prior_info)
{
    int j; /* Neighbor relative position to Pixel being updated */
    float delta, SurrogateCoeff;
    
    for (j = 0; j < 26; j++)
    {
        delta = icd_info->v - icd_info->neighbors[j];
        SurrogateCoeff = QGGMRF_SurrogateCoeff(delta,prior_info);
        
        icd_info->th1 += (icd_info->nb_wt[j]*SurrogateCoeff * delta);
        icd_info->th2 += icd_info->nb_wt[j]*SurrogateCoeff;
    }
    
}



void paraICD_Likelihood(struct GeomInfo *geom_info,struct ImgInfo *img_info,struct SourceLocInfo *source_loc_info, struct PriorInfo *prior_info,char **recon_mask,ENTRY *X,ENTRY *e,ENTRY *D,ENTRY *TildeV,int *order, float SigmaLambda, int it, struct Image *Image, ENTRY *consensus_X,struct ViewXYInfo *stored_view_xy_info,int DE_numprocs)
{
	int tid = omp_get_thread_num();
	int offset;  /* sjk */
	int j, jx, jy, jz, jzmax, Nxy, Nyz;
	float x, y, z, pixel, diff;
	struct ACol col_xyz;
	struct ViewXYZInfo view_xyz_info;
	struct ICDInfo icd_info;
	//struct ViewXYInfo view_xy_info;




	int k, l, Nsubit, jjx, jjy;
	float gamma;
	ENTRY **VSC;
	int *voxel_list;
	float *key;
	char zero_skip_flag;  
	int n;  



	createACol(&col_xyz, COL_LEN);		/* TODO COL_LEN hard-coded */
	//createViewXYInfo(&view_xy_info, geom_info);
	//view_xy_info.ic_start=stored_view_xy_info->ic_start;
	//view_xy_info.ic_num = stored_view_xy_info->ic_num;
	//view_xy_info.Mag = stored_view_xy_info->Mag;
	//view_xy_info.B = stored_view_xy_info->B;


//	fprintf(stdout,"after reach createACol\n");
//	fflush(stdout);




	createViewXYZInfo(&view_xyz_info, geom_info);

	VSC = (ENTRY **)get_img(img_info->Ny, img_info->Nx, sizeof(ENTRY));
	/* sjk */
        voxel_list = (int *)get_spc(img_info->Nx*img_info->Ny,sizeof(int));
        key = (float *)get_spc(img_info->Nx*img_info->Ny,sizeof(float));


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

		VSC[jx][jy] = 0.0;


		if (recon_mask[jx][jy])
		{
			x = img_info->x0 + jx*img_info->Del_xy;
			y = img_info->y0 + jy*img_info->Del_xy;

			
			//:compViewXYInfo_OnTheFly(jx,jy,x, y, &view_xy_info, geom_info, img_info, source_loc_info);

			offset= jx*Nyz + jy*img_info->Nz; 
			for (jz = tid*img_info->Nz/omp_get_num_threads(); jz < jzmax; jz++)
			{
				fillNeighbors(icd_info.neighbors,jx,jy,jz,X,img_info); 

				if( X[offset+jz] > 0 || positive_constraint ==0)
				  zero_skip_flag=0; 
				else
				{
				  zero_skip_flag=1; 
				 
				  if( (icd_info.neighbors[1] > 0) || (icd_info.neighbors[4] > 0) ||(icd_info.neighbors[7] > 0) ||(icd_info.neighbors[10] > 0) ||(icd_info.neighbors[15] > 0) ||(icd_info.neighbors[18] > 0) ||(icd_info.neighbors[21] > 0) ||(icd_info.neighbors[24] > 0) )
				  {
				    zero_skip_flag=0;
				  }
				}


				if(zero_skip_flag==0)
				{

				z = img_info->z0 + jz*img_info->Del_z;



				compViewXYZInfo(jx,jy,jz, z, &view_xyz_info, geom_info, img_info, source_loc_info, stored_view_xy_info);
				
				compAColxyzOnFly(jx,jy,jz,x, y, z, geom_info, source_loc_info, stored_view_xy_info, &view_xyz_info, &col_xyz,img_info);




				if (col_xyz.n_index > 0)
				{

					pixel = ICDStep_Likelihood(&icd_info, prior_info, jx, jy, jz, X, e, D, TildeV,geom_info->lambda0,img_info, &col_xyz, SigmaLambda,consensus_X,it,DE_numprocs);
					pixel = X[offset+jz] + damping_constant*(pixel - X[offset+jz]);
					/* clip */
					if(positive_constraint ==1){
						X[offset+jz] = ((pixel < 0.0) ? 0.0 : pixel);
					}
					else{
						/* no clipping for now try */
						X[offset+jz] = pixel; 
					}

					diff = X[offset+jz]-icd_info.v; 
					updateError(e, &(col_xyz), diff);
					VSC[jx][jy] += fabsf(diff);

					/*
					int myid;
					MPI_Comm_rank(MPI_COMM_WORLD, &myid);

					if(tid==0 && myid ==0){
						fprintf(stdout,"jx %d jy %d jz %d VSC %f \n",jx,jy,jz,VSC[jx][jy]);
						fflush(stdout);
					}
					*/
				}
				} 
			}
	//		freeViewXYInfoB(&view_xy_info, geom_info);		
		}
		
	}



	/* NH-ICD */

	Nsubit = 2;
	for (k = 0; k < Nsubit; k++)
	{

		gamma = 0.05;

                for (j = 0; j < Nxy; j++)
                {
                        jx = order[j] / img_info->Ny;
                        jy = order[j] % img_info->Ny;
                        key[j]=(float)(-VSC[jx][jy]);
                        voxel_list[j]=order[j];


                }

                quickselect(key,voxel_list,0,Nxy-1,(int)ceil(gamma*Nxy));


		for (l = 0; l < (int)ceil(gamma*Nxy); l++)
		{

                        jjx = voxel_list[l] / img_info->Ny;
                        jjy = voxel_list[l] % img_info->Ny;



			




			if (recon_mask[jjx][jjy])
			{
		 
				float total_diff = 0.0;
				x = img_info->x0 + jjx*img_info->Del_xy;
				y = img_info->y0 + jjy*img_info->Del_xy;


				
				//compViewXYInfo_OnTheFly(jjx,jjy,x, y, &view_xy_info, geom_info, img_info, source_loc_info);

				
				
				
				offset= jjx*Nyz + jjy*img_info->Nz; 
				for (jz = tid*img_info->Nz/omp_get_num_threads(); jz < jzmax; jz++)
				{
					fillNeighbors(icd_info.neighbors,jjx,jjy,jz,X,img_info);

					if( X[offset+jz] > 0 || positive_constraint ==0)
					  zero_skip_flag=0; 
					else
					{
				  		zero_skip_flag=1; 
				 
				  		if( (icd_info.neighbors[1] > 0) || (icd_info.neighbors[4] > 0) ||(icd_info.neighbors[7] > 0) ||(icd_info.neighbors[10] > 0) ||(icd_info.neighbors[15] > 0) ||(icd_info.neighbors[18] > 0) ||(icd_info.neighbors[21] > 0) ||(icd_info.neighbors[24] > 0) )
				  		{
				    			zero_skip_flag=0;
				  		}
					}


					if(zero_skip_flag==0)
					{

					z = img_info->z0 + jz*img_info->Del_z;



					compViewXYZInfo(jjx,jjy,jz, z, &view_xyz_info, geom_info, img_info, source_loc_info, stored_view_xy_info);


					compAColxyzOnFly(jjx,jjy,jz,x, y, z, geom_info, source_loc_info, stored_view_xy_info, &view_xyz_info, &col_xyz,img_info);
					


						if (col_xyz.n_index > 0)
						{
						pixel = ICDStep_Likelihood(&icd_info, prior_info, jjx, jjy, jz, X, e, D, TildeV, geom_info->lambda0, img_info, &col_xyz, SigmaLambda,consensus_X,it,DE_numprocs);
						pixel = X[offset+jz] + damping_constant*(pixel - X[offset+jz]);


						/* clip */
						if(positive_constraint ==1){
							X[offset+jz] = ((pixel < 0.0) ? 0.0 : pixel); 
						}
						else{
						/* no clipping for now try */
							X[offset+jz] = pixel; 
						}

						diff = X[offset+jz]-icd_info.v;

						updateError(e, &(col_xyz), diff);

						total_diff += fabsf(diff);


						}
					} 
				}
				if(total_diff >0)
					VSC[jjx][jjy] = total_diff;

				
				
	//			freeViewXYInfoB(&view_xy_info, geom_info);
						
				
			}
			
		}
	}

	freeViewXYZInfo(&view_xyz_info);
	freeACol(&col_xyz);
	free((void *)voxel_list);
	free((void *)key);
	free_img((void **)VSC);


}







void ICDReconstruct(
	struct Image *image,		/* reconstructed image */
	struct Image *V,
	struct Image *consensus_X,
	struct Image *TildeV,
	struct Image *VPrevious,
	struct Image *Vmean,
	struct Sinogram *sinogram, 	/* sinogram data */
	char **recon_mask,		/* XY reconstruction mask */
	struct PriorInfo *prior_info,	/* prior parameters */
	struct CEInfo *ce_info,
	char *fname,			/* output file name */
	int Nit,
	int myid,
	int NUMPROCS,
	int DE_id,
	int DE_numprocs,
	MPI_Comm *DE_comm,
	int DE_mode,
	int debug_mode, char *info_recon_dir)			/* number of iterations */
{
	int   it;				/* iteration index */
	float nlogpost;		/* negative log posterior */
	char  name[200];
	char  suffix[200];

	ENTRY *X;
	ENTRY *Y;
	ENTRY *e;
	ENTRY *D;  
	int *order;

	int t,i,j;
	struct SourceLocInfo source_loc_info;


	FILE *fp;
	time_t now;

	/* set X, Y, allocate memory for AX */
	X = image->img;
	Y = sinogram->sino;
	D = sinogram->D; 
 
	e = (ENTRY *)  get_spc((sinogram->geom_info.Nv)*(sinogram->geom_info.Nc)*(sinogram->geom_info.Nr), sizeof(ENTRY));
	/* generate random order index */
	
	order = (int *)get_spc(image->img_info.Nx*image->img_info.Ny, sizeof(int));
	shuffle(order, image->img_info.Nx*image->img_info.Ny);
	


	



	//float penalizer = 0.0;
	//for(i=0;i<image->img_info.Nx * image->img_info.Ny * image->img_info.Nz;i++)
    		//penalizer += (X[i]-TildeV->img[i])*(X[i]-TildeV->img[i])/(2*ce_info->SigmaLambda*ce_info->SigmaLambda);
  	//fprintf(stdout,"initial: penalizer is %f \n",penalizer); 





	/* clip X */
	/*voxels outside of recon_mask are already set to 0 in ct.c So no need to clip it here*/
	//clipImage(X, recon_mask, &(image->img_info));  /// Don't clip V, or Vmean or consensus_X,or TildeV! They should be negative!


	/* this is critical--projector skips pixels outside mask */
	

/* set estimation type */
	
	if (PnP_mode==0){	
		if (prior_info->q == 2.0 && prior_info->p < 2.0 && prior_info->p >= 1.0)
		{
			prior_info->est = QGGMRF;
			fprintf(stdout, "\nThis is MAP estimation using q-GGMRF!\n");
		}
		else
		{
			fprintf(stdout, "\nSoftware does not support the prior parameters provided!\n");
			exit(1);
		}
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
	compSourceLocInfo(&source_loc_info, &sinogram->geom_info,myid,NUMPROCS);


	//char errorFname[200];
	//sprintf(errorFname,"%s_%d.test_InitialError",fname,myid);		    	  


	/* for writing intial error sinogram */
	
	//experiment_debug_forwardProject(e, X, Y, recon_mask, &(sinogram->geom_info), &(image->img_info),myid);
	//writeSinogram_float(errorFname, e, sinogram->geom_info.Nr, sinogram->geom_info.Nc, sinogram->geom_info.Nv);







	/* write out the image after each iteration */
	/*
	if(PnP_mode){	
		strcpy(name, "write_");
		sprintf(suffix, "myid%d.vjk",myid);
		strcat(name, fname);
		strcat(name, suffix);
		writeImage(name, consensus_X);
	}
	else{
		strcpy(name, "write_");
		sprintf(suffix, "myid%d.vjk",myid);
		strcat(name, fname);
		strcat(name, suffix);
		writeImage(name, Vmean);
	}
	
	strcpy(name, "");
	sprintf(suffix, "_write_myid%d_X.vjk", myid);
	strcat(name, fname);
	strcat(name, suffix);
	writeImage(name, image);
	*/











	struct ViewXYInfo stored_view_xy_info;
	forwardProject(e, X, Y, recon_mask, &(sinogram->geom_info), &(image->img_info),&stored_view_xy_info,myid,NUMPROCS);
	
	//writeSinogram_float(errorFname, e, sinogram->geom_info.Nr, sinogram->geom_info.Nc, sinogram->geom_info.Nv);


/*
 	 readSinogram_float(errorFname, e, sinogram->geom_info.Nr * sinogram->geom_info.Nc* sinogram->geom_info.Nv);
*/	

	free((void *) Y);
	Y=NULL;


	printf("done projecting\n");
	fflush( stdout );
	
	//nlogpost = ProximalMapCostFunction3D_Likelihood(X,e,D,TildeV->img,prior_info,&(sinogram->geom_info),&(image->img_info), ce_info->SigmaLambda,myid);

	//TwoGMinusIOperator_PnP(TildeV, V, consensus_X,  Vmean, prior_info, recon_mask, ce_info->SigmaLambda, filter, order, myid,NUMPROCS);

	int num_threads=0;

    	if(!debug_mode){

		for (i =0; i < image->img_info.Nx; i++){
                	MPI_Allreduce(&(V->img[i*image->img_info.Ny*image->img_info.Nz]),&(Vmean->img[i*image->img_info.Ny*image->img_info.Nz]), image->img_info.Ny * image->img_info.Nz, MPI_FLOAT, MPI_SUM, *DE_comm);
    		}


		#pragma omp parallel
		{
		num_threads = omp_get_num_threads();
		#pragma omp for
    		for (i =0;i<image->img_info.Nx*image->img_info.Ny*image->img_info.Nz;i++){
        		Vmean->img[i] = Vmean->img[i]/DE_numprocs;
    		}
		}
	}
	else{
		fprintf(stdout,"Debug Mode Turned on ! \n");
		fflush(stdout);

		#pragma omp parallel
		{
		num_threads = omp_get_num_threads();
		#pragma omp for
 		for (i =0;i<image->img_info.Nx*image->img_info.Ny*image->img_info.Nz;i++){
        		Vmean->img[i] = V->img[i];
    		}
		}
	}

		/* ICD iteration starts here */
	fprintf(stdout, "\nstart ICD iterations...\n");
	fflush(stdout);



	for (it = 0; it < Nit; it++)
	{

		int i=0;
	        shuffle(order, image->img_info.Nx*image->img_info.Ny);

		fprintf(stdout, "it = %d\n", it+1);
		fflush(stdout);

		for (i=0; i<V->img_info.Nz; i++){
    	 		memcpy(&VPrevious->img[i*image->img_info.Nx *image->img_info.Ny],&V->img[i*image->img_info.Nx*image->img_info.Ny],sizeof(float)*image->img_info.Nx * image->img_info.Ny);
    		}


		#pragma omp parallel
		{
			paraICD_Likelihood(&sinogram->geom_info,&image->img_info,&source_loc_info,prior_info,recon_mask,X,e,D,TildeV->img,order, ce_info->SigmaLambda,it, image,consensus_X->img,&stored_view_xy_info,DE_numprocs);
		}
		

		/* compute negative log posterior */
	
		//nlogpost = ProximalMapCostFunction3D_Likelihood(X,e,D,TildeV->img,prior_info,&(sinogram->geom_info),&(image->img_info),ce_info->SigmaLambda,myid);


                for (i=0; i<V->img_info.Nz * V->img_info.Ny * V->img_info.Nx; i++){

                	V->img[i] = 2*X[i] -TildeV->img[i];
		}



                for (i=0; i<V->img_info.Nz * V->img_info.Ny * V->img_info.Nx; i++){
                        V->img[i] = ce_info->consensus_rho*V->img[i] + (1-ce_info->consensus_rho)*VPrevious->img[i];
                }


		
		TwoGMinusIOperator_PnP(TildeV, V, consensus_X, Vmean,prior_info, recon_mask, ce_info->SigmaLambda, order, myid, NUMPROCS,it,DE_id,DE_numprocs,DE_comm,debug_mode);


		float cost = Consensus_Cost(consensus_X,image,Vmean,NUMPROCS,DE_id,DE_numprocs,DE_comm);
		
		//nlogpost = ProximalMapCostFunction3D_Likelihood(X,e,D,TildeV->img,prior_info,&(sinogram->geom_info),&(image->img_info),ce_info->SigmaLambda,myid);

		if(DE_id==0 && PnP_mode==1){
			fprintf(stdout,"myid %d: it %d cost %f prior %f\n",myid,it,cost,ProximalMapCostFunction3D_Prior(consensus_X,prior_info,myid));
			fflush(stdout);
		}
		else if(DE_id==0 && PnP_mode==0){
			fprintf(stdout,"myid %d: it %d cost %f prior %f\n",myid,it,cost,ProximalMapCostFunction3D_Prior(Vmean,prior_info,myid));
			fflush(stdout);
		}


	

		if (myid==0){
			/* write out the image after each iteration */
			if(PnP_mode){	
				strcpy(name, "");
				sprintf(suffix, "myid%d.vjk",myid);
				strcat(name, fname);
				strcat(name, suffix);
				writeImage(name, consensus_X);
			}
			else{
				strcpy(name, "");
				sprintf(suffix, "myid%d.vjk",myid);
				strcat(name, fname);
				strcat(name, suffix);
				writeImage(name, Vmean);
			}
		}

		MPI_Barrier(*DE_comm);

	}

	fprintf(stdout, "\nfinish ICD iterations!\n");
	fflush(stdout);

	freeSourceLocInfo(&source_loc_info);

	fprintf(stdout,"freeSource \n");
	fflush(stdout);

	free(e);

	fprintf(stdout,"free e\n");
	fflush(stdout);

        free(order);

	fprintf(stdout,"free order\n");
	fflush(stdout);
	

	freeViewXYInfo_stored(&stored_view_xy_info, &sinogram->geom_info,&image->img_info,recon_mask);

}


/* unified q-GMMRF potential function */
/* including GGMRF, quad */
/* the potential function of the QGGMRF prior model.  p << q <= 2 */
float QGGMRF_Potential(float delta, struct PriorInfo *Rparams)
{
    float p, q, T, SigmaX;
    float temp, GGMRF_Pot;
    
    p = Rparams->p;
    q = Rparams->q;
    T = Rparams->T;
    SigmaX = Rparams->SigmaX;
    
    GGMRF_Pot = pow(fabs(delta),p)/(p*pow(SigmaX,p));
    temp = pow(fabs(delta/(T*SigmaX)), q-p);
    
    return ( GGMRF_Pot * temp/(1.0+temp) );
}
/* Quadratic Surrogate Function for the log(prior model) */
/* For a given convex potential function rho(delta) ... */



/* The surrogate function defined about a point "delta_p", Q(delta ; delta_p), is given by ... */
/* Q(delta ; delta_p) = a(delta_p) * (delta^2/2), where the coefficient a(delta_p) is ... */
/* a(delta_p) = [ rho'(delta_p)/delta_p ]   ... */
/* for the case delta_current is Non-Zero and rho' is the 1st derivative of the potential function */
/* Return this coefficient a(delta_p) */
/* Prior-model specific, independent of neighborhood */

float QGGMRF_SurrogateCoeff(float delta, struct PriorInfo* reconparams)
{
    float p, q, T, SigmaX, qmp;
    float num, denom, temp;
    
    p = reconparams->p;
    q = reconparams->q;
    T = reconparams->T;
    SigmaX = reconparams->SigmaX;
    qmp = q - p;
    
    /* Refer to Chapter 7, MBIR Textbook by Prof Bouman, Page 151 */
    /* Table on Quadratic surrogate functions for different prior models */
    
    if (delta == 0.0)
    return 2.0/( p*pow(SigmaX,q)*pow(T,qmp) ) ; /* rho"(0) */
    
    temp = pow(fabs(delta/(T*SigmaX)), qmp);
    num = q/p + temp;
    denom = pow(SigmaX,p) * (1.0+temp) * (1.0+temp);
    num = num * pow(fabs(delta),p-2) * temp;
    
    return num/denom;
}


