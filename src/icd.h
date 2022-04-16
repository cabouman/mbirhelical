/* modified by P.J 03/10/2011 */

#ifndef _ICD_H_
#define _ICD_H_
#include "mpi.h"
//#include <pthread.h>

#define WA	0.06		/* weight for closest neighbors */
#define WB	0.04		/* weight for closer neighbors */
#define WC	0.02		/* weight for close neighbors */
#define TOL	10e-9		/* half-interval seaerch equation solver tolerance */
#define sign(a,b) ((a) > (b) ? (1.0) : (-1.0))		/* sign(a-b) */
#define abso(a) ((a) > 0 ? (a) : (-1.0*a))
#define PnP_mode 0
#define damping_constant 1.0

struct ICDInfo
{
	float v;		/* old center pixel value */
	float th1;
	float th2;
	float neighbors[26];	/* 3D neighbors, 0-5 closest, 6-17 closer, 18-25 close */
	float nb_wt[26];	/* weighting coefficiets, replacing WA,WB,WC because voxels not always isometric */
};

void fillNeighbors(float *neighbors, int jx, int jy, int jz, ENTRY *X, struct ImgInfo *img_info);

void compNeighborWeight(float *nb_wt, struct ImgInfo *img_info);

void compNeighborWeightNonhomogeneous(
        float *nb_wt,
        struct ImgInfo *img_info,
        struct GeomInfo *geom_info,
        struct SourceLocInfo *source_loc_info,
        float x,
        float y,
        float z);

void updateError(ENTRY *e, struct ACol *col_xyz, float diff);

void MergeAndDenoise(   struct Image *V,
			struct Image *VMean,                       /* Input V */
                        struct Image *consensus_X,                       /* Output Z - merge V and denoise. If denoising is an optimization problem, Z acts as input (initial state) as well */
                        struct PriorInfo *prior_info, char **reconMask, float SigmaLambda, int *order,int myid, int NUMPROCS, int it,int DE_id,int DE_numprocs,MPI_Comm *DE_comm,int debug_mode);


void StackAndReflect(
                            struct Image *TildeV,              /* Output = 2*Stacked(Z)-V , i.e. Reflect V about Z.  (Z is computed before hand) */
                            struct Image *V,                       /* Input V */
                            struct Image *consensus_X , struct Image *Vmean );


void SolveProximalMap_Prior(struct Image *Image,
                            char  ** reconMask,
                            struct Image *ProximalMapInput,
                            float  SigmaLambda, struct PriorInfo *prior_info, int *order,int it,int myid);



float ProximalMapCostFunction3D_Prior(struct Image *Image,
                                      struct PriorInfo *prior_info,int myid);


float QGGMRF_SurrogateCoeff(float delta, struct PriorInfo* reconparams);

float QGGMRF_Potential(float delta, struct PriorInfo *Rparams);



float ICDStep_PriorOnly(
        struct ICDInfo *icd_info,
        struct PriorInfo *prior_info,
        int jx,
        int jy,
        int jz,
        ENTRY *X,                       /* image, attenuation coefficients */
	struct Image *ProximalMapInput,
        struct ImgInfo *img_info,        /* image parameters */
	float SigmaLambda);


void paraICD_Prior(struct ImgInfo *img_info,struct PriorInfo *prior_info,char **recon_mask,float *X,struct Image *ProximalMapInput, int *order, float SigmaLambda,int it, struct Image * Image);


void QGGMRF3D_UpdateICDParams(struct ICDInfo *icd_info, struct PriorInfo *prior_info);




float negLogPosterior(
	ENTRY *X,
	ENTRY *Y,
	ENTRY *AX,
	ENTRY *D,
	struct PriorInfo *prior_info,
	struct GeomInfo *geom_info,
	struct ImgInfo *img_info);

float ICDStep_Likelihood(
	struct ICDInfo *icd_info,	/* ICD parameters */
	struct PriorInfo *prior_info,	/* prior parameters */
	int jx,
	int jy,
	int jz,
	ENTRY *X,			/* image, attenuation coefficients */
	ENTRY *e,			/* line integral */
	ENTRY *D,			/* noise matrix */
	ENTRY *TildeV,
	float lambda0,
	struct ImgInfo *img_info,	/* image parameters */
	struct ACol *col_xyz,
	float SigmaLambda,
	ENTRY *consensus_X,int it,int DE_numprocs);		/* A column for this (jx, jy, jz) pixel */



void ICDReconstruct(
	struct Image *image,		/* reconstructed image */
	struct Image *V,
	struct Image *consensus_X,
	struct Image *TildeV,
	struct Image *VPrevious,
	struct Image *Vmean,
	struct Sinogram *sinogram,	/* sinogram data */
	char **recon_mask,		/* XY reconstruction mask */
	struct PriorInfo *prior_info,	/* prior parameters */
	struct CEInfo *ce_info,
	char *fname,			/* output file name */
	int Nit,
	int myid,int NUMPROCS,
	int DE_id,int DE_numprocs,
	MPI_Comm *DE_comm, int DE_mode, int debug_mode,char *info_recon_dir);			/* number of iterations */



float quadUpdate(struct ICDInfo *icd_info);

float qpotential(float delta, float q, float p, float c);
float qinfluence(float delta, float q, float p, float c);
float subFunctionUpdate(struct ICDInfo *icd_info);

float pCost(float x, struct ICDInfo *icd_info);
float pCostDerivative(float x, void *data);
float halfIntervalUpdate(struct ICDInfo *icd_info);

void paraICD_Likelihood(struct GeomInfo *geom_info, struct ImgInfo *img_info,struct SourceLocInfo *source_loc_info, struct PriorInfo *prior_info,char **recon_mask,ENTRY *X,ENTRY *e,ENTRY *D, ENTRY *TildeV,int *order, float SigmaLambda,int it, struct Image* Image, ENTRY *consensus_X,struct ViewXYInfo *view_xy_info,int DE_numprocs);


void TwoGMinusIOperator_PnP(
                                    struct Image *TildeV,                   /* Output (2G-I)V, where G = merge, denoise & stack operator here */
                                    struct Image *V,                        /* Input V */
                                    struct Image *consensus_X,                         /* Output Z - merge V and denoise. If denoising is an optimization problem, Z acts as input (initial state) as well */
                                    struct Image *Vmean, struct PriorInfo *prior_info, char **reconMask, float SigmaLambda, int *order, int myid, int NUMPROCS,int it,int DE_id,int DE_numprocs,MPI_Comm *DE_comm,int debug_mode);


#endif
