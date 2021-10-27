#ifndef _RANDLIB_H_
#define _RANDLIB_H_

float random2();
int random3();
void srandom2(unsigned long num);
void readseed();
void writeseed();
float normal();
float dexprand();
int poisson(float lambda);

#endif
