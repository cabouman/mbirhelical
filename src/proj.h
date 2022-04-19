/* version 6.0, P.Jin 07/13/2012 */

#ifndef _PROJ_H_
#define _PROJ_H_

//#include <pthread.h>

void createSinogram(struct Sinogram *sinogram);

void freeSinogram(struct Sinogram *sinogram);

void createImage(struct Image *image);

void freeImage(struct Image *image);

void fillGeomInfo(struct GeomInfo *geom_info);

void fillImgInfo(struct ImgInfo *img_info);

void checkInfo(struct GeomInfo *geom_info, struct ImgInfo *img_info);  /* sjk */

void createSourceLocInfo(
	struct SourceLocInfo *source_loc_info,
	struct GeomInfo *geom_info);

void compSourceLocInfo(
	struct SourceLocInfo *source_loc_info,
	struct GeomInfo *geom_info,int myid, int total_nodes);

void freeSourceLocInfo(struct SourceLocInfo *source_loc_info);

void createViewXYInfo(
	struct ViewXYInfo *view_xy_info,
	struct GeomInfo *gemo_info);

void createViewXYInfo_stored(
	struct ViewXYInfo *view_xy_info,
	struct GeomInfo *geom_info, struct ImgInfo * img_info,char **recon_mask);

void compViewXYInfo_OnTheFly(
	int jx,
	int jy,
	float x,
	float y,
	struct ViewXYInfo *view_xy_info,
	struct GeomInfo *geom_info,
	struct ImgInfo *img_info,
	struct SourceLocInfo *source_loc_info);


void compViewXYInfo(
	int jx,
	int jy,
	float x,
	float y,
	struct ViewXYInfo *view_xy_info,
	struct GeomInfo *geom_info,
	struct ImgInfo *img_info,
	struct SourceLocInfo *source_loc_info);

void freeViewXYInfoB(
	struct ViewXYInfo *view_xy_info,
	struct GeomInfo *geom_info);

void freeViewXYInfo(
	struct ViewXYInfo *view_xy_info,
	struct GeomInfo *geom_info);


void freeViewXYInfo_stored(
	struct ViewXYInfo *view_xy_info,struct GeomInfo *geom_info,
	struct ImgInfo *img_info,char ** recon_mask);



void createViewXYZInfo(
	struct ViewXYZInfo *view_xyz_info,
	struct GeomInfo *geom_info);


void compViewXYZInfo(
	int jx,
	int jy,
	int jz,
	float z,
	struct ViewXYZInfo *view_xyz_info,
	struct GeomInfo *geom_info,
	struct ImgInfo *img_info,
	struct SourceLocInfo *source_loc_info,
	struct ViewXYInfo *view_xy_info);


void freeViewXYZInfo(struct ViewXYZInfo *view_xyz_info);

void createACol(struct ACol *A_col, int length);

void increaseAColLength(struct ACol *A_col);  /* sjk */

void freeACol(struct ACol *A_col);


void compAColxyzOnFly(
	int jx,
	int jy,
	int jz,
	float x,
	float y,
	float z,
	struct GeomInfo *geom_info,
	struct SourceLocInfo *source_loc_info,
	struct ViewXYInfo *view_xy_info,
	struct ViewXYZInfo *view_xyz_info,
	struct ACol *A_col,struct ImgInfo *img_info);

/*void forwardProject(ENTRY *AX, ENTRY *X, char **recon_mask, struct GeomInfo *geom_info, struct ImgInfo *img_info);*/


void metal_forward(ENTRY *e, ENTRY *X, char **recon_mask, struct GeomInfo *geom_info, struct ImgInfo *img_info,struct ViewXYInfo *view_xy_info, int myid,int total_nodes);


void forwardProject(ENTRY *e, ENTRY *X, ENTRY *Y,  char **recon_mask, struct GeomInfo *geom_info, struct ImgInfo *img_info,struct ViewXYInfo *view_xy_info,int myid,int total_nodes);

/* For parallel computing */
/*struct paraForwardProjectData
{
	int tid;
	struct GeomInfo *geom_info;
	struct ImgInfo *img_info;
	struct SourceLocInfo *source_loc_info;
	ENTRY *X ;
	ENTRY *AX;
	unsigned short *AX_mask;
	char **recon_mask;
	pthread_mutex_t *AX_mutex;
};
*/

void paraForwardProject(struct GeomInfo *geom_info,struct ImgInfo *img_info,struct SourceLocInfo *source_loc_info,ENTRY *X,ENTRY *e,char **recon_mask,struct ViewXYInfo *view_xy_info);

void serialForwardProject(ENTRY *AX, ENTRY *X, struct GeomInfo *geom_info, struct ImgInfo *img_info);

void backProject(ENTRY *AX, ENTRY *X, struct GeomInfo *geom_info, struct ImgInfo *img_info,int myid);

void *paraBackProject(struct GeomInfo *geom_info,struct ImgInfo *img_info,struct SourceLocInfo *source_loc_info,ENTRY *X,ENTRY *AX);

#endif
