/* version 6.0, P.Jin 07/13/2012 */

#ifndef _DATA_H_
#define _DATA_H_

/* we need signed indicies to do comparison with 0 */
typedef short CHANNEL;		/* two bytes */
typedef char PROCHANNEL;	/* num of channels in projection of a voxel */
typedef char ROW;		/* one byte, TODO signed char enough? */
typedef char PROROW;		/* num of rows in projection of a voxel */
typedef int VIEW;		/* four bytes */
typedef short PROVIEW;
typedef float ENTRY;		/* four bytes, entries of arrays */

#define PI		3.1415926535897932384
#define EPSILON		10e-20				/* A entry tolerance */
#define COL_LEN		30000				/* A column num of non-zero entries TODO hard-coded */
#define min(a,b)	((a) < (b) ? (a) : (b))
#define max(a,b)	((a) > (b) ? (a) : (b))
#define clip(a,b,c)	(min(max((a), (b)), (c)))
#define adjust(a)	((a) > -PI ? ((a)-((int)(((a)+PI)/(2*PI)))*2*PI) : ((a)-((int)(((a)+PI)/(2*PI)))*2*PI)+2*PI)

struct GeomInfo
{
	/* -------------------------------------------- */
	/*		INPUT PARAMETERS		*/
	/* -------------------------------------------- */
	int Nr;			/* num of rows */
	int Nc;			/* num of channels */
	int Nv;			/* num of views */
	int Nvpr;		/* num of views per rotation */
	float r_si;		/* (mm) source to iso */
	float r_sd;		/* (mm) source to detector */
	float u;		/* pitch, rows per rotation */
	float beta0;		/* (rad) initial view angle */
	float Del_beta;	/* (rad) view angle spacing */
	float Del_alphac;	/* (rad) detector angle spacing in channel direction */
	float del_alphac;	/* (rad) sampling offset in channel direction */
	float Del_dr;		/* (mm) detector spacing in row direction */
	float fov;		/* diameter of circular field of view */
	float lambda0;		/* emission photon rate, used to calculate D, if lambda0 = 0, D = I */
	float evar;		/* electronic noise variance */
	char sinoFile[100];
	char doseFile[100];
	char offsetFile[100];
	char detectorsFile[100];

	/* -------------------------------------------- */
	/*		PRECOMPUTED PARAMETERS		*/
	/* -------------------------------------------- */
	float alphac0;		/* (rad) channel angle of the detector of channel 0, = -(Nc-1)*Del_alphac/2 + del_alphac */
	float Del_zs;		/* (mm) z shift per view, = u*Del_dr*Del_beta/(2*PI) */
	float detc;		/* (rad) radius of detector array in channel direction, = Nc*Del_alphac */
	float detr;		/* (mm) width of detector array in row direction, = Nr*Del_dr */ 
	float half_detr;	/* (mm) half width of detector array in row direction, = (Nr-1)*Del_dr/2 */
	float cone_zbuffer;	/* sjk: (mm) voxels w/ z outsize zs+/-cone_zbuffer will not fall within cone, =detr*(r_si+fov/2)/(2*r_sd) */
};

struct Sinogram
{
	struct GeomInfo geom_info;
	ENTRY *sino;		/* 1D array to store sinogram data */
	ENTRY *counts;		/* array to hold photon counts (if used) */
	ENTRY *dose;		/* array to hold dosage (if used) */
	ENTRY *offset;		/* sensor background current (if used) */
	ENTRY *D;		/* normalized noise matrix entries (if used) */
};

struct ImgInfo
{
	/* -------------------------------------------- */
	/*		INPUT PARAMETERS		*/
	/* -------------------------------------------- */
	int Nx;			/* num of voxels in x direction */
	int Ny;			/* num of voxels in y direction */
	int Nz_mid;
	int Nz;			/* num of voxels in z direction */
	float xc;		/* center voxel */
	float yc;
	float zc;
	float Del_xy;		/* (mm) voxel spacing in x, y direction, assume square in xy plane */
	float Del_z;		/* (mm) voxel spacing in z direction */
	float rI;		/* (mm) radius of reconstruction mask */
	char imgFile[100];
	char maskFile[100];	/* reconstruction mask */

	/* -------------------------------------------- */
	/*		PRECOMPUTED PARAMETERS		*/
	/* -------------------------------------------- */
	float x0;		/* = xc - Del_xy*(Nx-1)/2 */
	float y0;		/* = yc - Del_xy*(Ny-1)/2 */
	float z0;		/* = zc - Del_z*(Nz-1)/2 */
};

struct Image
{
	struct ImgInfo img_info;
	ENTRY *img;		/* (mm-1) 1D array to store image data */
};

struct SourceLocInfo		/* indepenedent of (x,y,z) */
{
	ENTRY *beta;		/* view angles, beta[iv] = beta0 + iv*Del_beta, 0 <= iv < Nv */
	ENTRY *xs;		/* x coord. of the source, xs[iv], 0 <= iv < Nv */
	ENTRY *ys;		/* y coord. of the source, ys[iv], 0 <= iv < Nv */
	ENTRY *zs;		/* z coord. of the source, zs[iv], 0 <= iv < Nv */
};

struct ViewXYInfo		/* depend on (x,y), independent of z */
{
	CHANNEL	*ic_start;	/* ic_start[iv], 0 <= iv < Nv */
	PROCHANNEL *ic_num;	/* ic_num[iv], 0 <= iv < Nv */
	ENTRY *Mag;		/* Mag[iv], 0 <= iv < Nv */
	ENTRY *Wr;		/* Wr[iv], 0 <= iv < Nv */
	ENTRY **B;		/* B[iv][p], 0 <= iv < Nv, 0 <= p < ic_num[iv] */
};

struct ViewXYZInfo		/* depend on (x,y,z) */
{
	VIEW iv_start;
	PROVIEW iv_num;
	ROW *ir_start;		/* ir_start[iv], 0 <= iv < Nv */
	PROROW *ir_num;		/* ir_num[iv], 0 <= iv < Nv */
};

struct ACol			/* column of A matrix */
{
	int n_index;		/* four bytes */
	int *index;
	ENTRY *val;
	int array_length;       /* sjk: keep note on full array size */
};

#define ETYPE	char		/* estimation type */
#define ML	-1		/* ML estimation */
#define QUAD	0		/* MAP estimation using quadratic prior, i.e. p = 2 */
#define QGGMRF	1		/* MAP estimation using q-GMMRF prior, i.e. 1 <= q < p = 2 */
#define PGGMRF	2		/* MAP estimation using p-GGMRF prior, i.e. 1 <= p < 2 */

struct PriorInfo
{
	ETYPE est;		/* estimation type */
	float q;		/* q-GMMRF q */
	float p;		/* q-GMMRF p */
	float c;		/* q-GMMRF c */
	float sigma;		/* q-GMMRF sigma */
};


#endif
