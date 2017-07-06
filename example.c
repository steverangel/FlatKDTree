#include <stdlib.h>
#include <stdio.h>
#include "flatKDTree.h"

void print_kd_arr_w_idx(REAL** A, size_t n, size_t *Idx) {
  size_t i,d;
  for(i=0;i<n;++i) {
    for(d=0;d<DIM;++d)
      printf("%f\t",A[d][Idx[i]]);
    printf("\n");
  }
}
void print_kd_arr(REAL** A, size_t n) {
  size_t i,d;
  for(i=0;i<n;++i) {
    for(d=0;d<DIM;++d)
      printf("%f\t",A[d][i]);
    printf("\n");
  }
}

int main(int argc, char** argv)
{
  printf("FlatKDTree Example\n");

  size_t n = 64;
  REAL **A = (REAL**)malloc(DIM*sizeof(REAL*));
  size_t *Idx = (size_t*)malloc(n*sizeof(size_t));

  unsigned short d;
  for(d=0;d<DIM;++d) {
    A[d] = (REAL*)malloc(n*sizeof(REAL)); 
    size_t i;
    for(i=0;i<n;++i)
      A[d][i] = (REAL)rand()/(REAL)RAND_MAX;
  }

  //print_kd_arr(A,n);

  size_t n_unique = n;
  n_unique = unique_rows(&A,DIM,n);
  size_t i;
  for(i=0;i<n_unique;++i) {
    Idx[i]=i;
  }
  tree_order(A,&Idx,0,0,n_unique-1);
  printf("%ld unique\n",n_unique);

  print_kd_arr_w_idx(A,n_unique,Idx);

  for(d=0;d<DIM;++d)
    free(A[d]); 
  return 0;
}
