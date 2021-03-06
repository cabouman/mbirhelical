#ifndef _ALLOCATE_H_
#define _ALLOCATE_H_

void *get_spc(int num, size_t size);
void *mget_spc(int num, size_t size);
void **get_img(int wd,int ht, size_t size);
void free_img(void **pt);
void *multialloc(size_t s, int d, ...);
void multifree(void *r,int d);

#endif
