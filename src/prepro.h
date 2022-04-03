/* version 6.0, P.Jin 07/13/2012 */

#ifndef _PREPRO_H_
#define _PREPRO_H_

void constImage(ENTRY *X, char **recon_mask, float value, struct ImgInfo *img_info);
void clipImage(ENTRY *X, char **recon_mask, struct ImgInfo *img_info);
void addWhiteNoise(ENTRY *data, int len, float noise_std);
void addNoise(ENTRY *data, float lambda0, int len);

void createReconMask(char ***recon_mask, struct ImgInfo *img_info);
void compReconMask(char **recon_mask, struct ImgInfo *img_info);
void freeReconMask(char **recon_mask);

void initImage(struct GeomInfo *geom_info, struct Image *image, ENTRY value);

void shuffle(int *array, int len);

//void filterUMM(ENTRY **UMM, ENTRY **VSC, ENTRY **filter, int h, int w);

/* Bouman's HU definition, 0 <-> air, 1000 <-> water */

#endif
