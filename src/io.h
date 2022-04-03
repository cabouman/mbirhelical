/* version 6.0, P.Jin 07/13/2012 */

#ifndef _IO_H_
#define _IO_H_

void readImgInfo(char *fname, struct ImgInfo *img_info);
void writeImgInfo(char *fname, struct ImgInfo *img_info);
void printImgInfo(struct ImgInfo *img_info);
void readImage(char *fname, struct Image *image);
void writeImage(char *fname, struct Image *image);
void readReconMask(char **recon_mask, struct ImgInfo *img_info);
void readGeomInfo(char *fname, int total_nodes,struct GeomInfo *geom_info);
void writeGeomInfo(char *fname, struct GeomInfo *geom_info);
void printGeomInfo(struct GeomInfo *geom_info);
void readSinogram(char *fname, struct Sinogram *sinogram);
ENTRY* readSinogram_new(char *fname, int length);
ENTRY* readSinogram_float(char *fname, ENTRY *e, int length);
void fillSinogramData(struct Sinogram *sinogram, int num_nodes, int myid);
void writeSinogram(char *fname, struct Sinogram *sinogram);
void writeSinogram_float(char *fname,ENTRY *Y,int Nr,int Nc,int Nv);
void readPrior(char *fname, struct PriorInfo *prior_info);
void readCE(char *fname, struct CEInfo *ce_info);
void readAll_GeomDirectory(char *fname, int my_node_index, int total_nodes, struct GeomInfo *geom_info);

#endif
