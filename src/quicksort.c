/*
 *  program: quicksort.c
 *  written by:  Jordan Kisner, kisner@purdue.edu 
 *  updated: 5/20/12
 */ 

#include "quicksort.h"

/* external functions */
/* void quicksort(float *key, int *list, int p, int r); */
/* void quickselect(float *key, int *list, int p, int r, int i); */

/* internal function */
int partition(float *key, int *list, int p, int r);

/********** quicksort() ***********
* Sorts "list" array according according to the "key" array
* Note the partition routine is set up for ascending order.
* See "Intro to Algorithms" by Cormen, Leiserson, and Rivest.
* p and r are the starting and ending array indices. 
*/
void quicksort(key, list, p, r)
float *key;
int *list;
int p,r;
{
  int q;
  if(p<r) {
    q=partition(key, list, p, r);
    quicksort(key, list, p, q);
    quicksort(key, list, q+1, r);
  }
}


/********** quickselect() ***********
* Not a full sort but rearranges the arrays so that 
* the rank i element falls at key[i-1], and everything in
* key[i..r] are all >= key[i-1], and key[p..(i-1)] are
* all <= key[i-1].  Expected complexity is O(n).
*/
void quickselect(key,list,p,r,i)
float *key;
int *list;
int p,r,i;
{
  int q,k;

  if(p==r) 
    return;
  q=partition(key,list,p,r);
  k=q-p+1;
  if(i<=k)
    quickselect(key,list,p,q,i);
  else
    quickselect(key,list,q+1,r,i-k);
}


/********** partition() ***********
* Routine accompanying the quicksort function.
* It sorts two arrays at the same time.
*/
int partition(key, list, p, r)
float *key;
int *list;
int p,r;
{
  float x,tempf;
  int i,j,tempi;

  x=key[p];
  i=p-1;
  j=r+1;
  while(1) {
    do j--;
    while(key[j]>x);  /* ">" for ascending; "<" for decending */
    do i++;
    while(key[i]<x);  /* "<" for ascending; ">" for decending */
    if(i<j) {
      tempf=key[i]; 
      key[i]=key[j]; 
      key[j]=tempf; 
      tempi=list[i]; 
      list[i]=list[j]; 
      list[j]=tempi; 
    }
    else return j;
  }
}


