#include <stdlib.h>
#include <stdio.h>
#include "spline.h"
#include<time.h>

int main(){
  double x[500];
  double y[2000][500];
  int n = 500;
  double c[500];
  double scratch[5000];
  
  FILE* fx = fopen("fx.dat", "r");
  FILE* x_grid = fopen("x_grid.dat", "r");
  
  for(int i = 0; i<2000; i++){
    for(int j = 0; j<500; j++){
      double entry;
      fscanf(fx,"%lf", &entry);
      y[i][j] = entry; 
    }
  }

  for(int i =0; i<500; i++){
    double entry;
    fscanf(x_grid, "%lf", &entry);
    x[i] = entry;
  }

  clock_t begin = clock();
  struct spline test;
  test.n = n;
  test.x = x;
  
  for(int i =0 ; i<2000; i++){
     test.y = y[i];
     test.c = c;
     spline_fit_natural(&test, scratch);
  }
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  

  printf("Time: %f", time_spent);
  
  return 0;

}
