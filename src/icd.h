/* modified by P.J 03/10/2011 */

#ifndef _ICD_H_
#define _ICD_H_

//#include <pthread.h>

#define WA	0.06		/* weight for closest neighbors */
#define WB	0.04		/* weight for closer neighbors */
#define WC	0.02		/* weight for close neighbors */
#define TOL	10e-9		/* half-interval seaerch equation solver tolerance */
#define sign(a,b) ((a) > (b) ? (1.0) : (-1.0))		/* sign(a-b) */
#define abso(a) ((a) > 0 ? (a) : (-1.0*a))

struct ICDInfo
{
	float q;		/* q-GGMRF power q */
	float p;		/* q-GGMRF power p */
	float c;		/* q-GGMRF c */
	float sigma;		/* q-GGMRF sigma */
	float sigma2q;		/* q-GGMRF sigma^p */
	float c2qmp;		/* q-GGMRF c^(q-p).. comp saver */
	float th1;		/* theta1 */
	float th2;		/* theta2 */
	float v;		/* old center pixel value */
	float ml;		/* ML estimate of the pixel value */
	float low;		/* lower bound for pixel update */
	float high;		/* upper bound for pixel update */
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

void updateAX(ENTRY *AX, struct ACol *col_xyz, float diff);

float negLogPosterior(
	ENTRY *X,
	ENTRY *Y,
	ENTRY *AX,
	ENTRY *D,
	struct PriorInfo *prior_info,
	struct GeomInfo *geom_info,
	struct ImgInfo *img_info);

float ICDStep(
	struct ICDInfo *icd_info,	/* ICD parameters */
	struct PriorInfo *prior_info,	/* prior parameters */
	int jx,
	int jy,
	int jz,
	ENTRY *X,			/* image, attenuation coefficients */
	ENTRY *Y,			/* line integral */
	ENTRY *AX,			/* current AX */
	unsigned short *AX_mask,
	ENTRY *D,			/* noise matrix */
	float lambda0,
	struct ImgInfo *img_info,	/* image parameters */
	struct ACol *col_xyz);		/* A column for this (jx, jy, jz) pixel */

void ICDReconstruct(
	struct Image *image,		/* reconstructed image */
	struct Sinogram *sinogram,	/* sinogram data */
	char **recon_mask,		/* XY reconstruction mask */
	struct PriorInfo *prior_info,	/* prior parameters */
	char *fname,			/* output file name */
	int Nit);			/* number of iterations */

float quadUpdate(struct ICDInfo *icd_info);

float qpotential(float delta, float q, float p, float c);
float qinfluence(float delta, float q, float p, float c);
float subFunctionUpdate(struct ICDInfo *icd_info);

float pCost(float x, struct ICDInfo *icd_info);
float pCostDerivative(float x, void *data);
float halfIntervalUpdate(struct ICDInfo *icd_info);

/* For parallel computing */
/*
struct paraICDReconstructData
{
	int tid;
	struct GeomInfo *geom_info;
	struct ImgInfo *img_info;
	struct SourceLocInfo *source_loc_info;
	struct PriorInfo *prior_info;
	char **recon_mask;
	ENTRY *X;
	ENTRY *Y;
	ENTRY *AX;
	unsigned short *AX_mask;
	ENTRY *D;  
	ENTRY **filter;
	int *order;
	pthread_mutex_t *AX_mutex;
};
*/

void *paraICDReconstruct(struct GeomInfo *geom_info,struct ImgInfo *img_info,struct SourceLocInfo *source_loc_info,struct PriorInfo *prior_info,char **recon_mask,ENTRY *X,ENTRY *Y,ENTRY *AX,unsigned short *AX_mask,ENTRY *D,ENTRY **filter,int *order);

#endif
