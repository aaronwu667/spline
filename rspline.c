#include <stdio.h>
#include <stdlib.h>
#include "spline.h"
#include <assert.h>

/*
 * O(n) Tridiagonal system solver: Thomas algorithm.
 * Note: not guaranteed to be stable and destroys original input. 
 * Reference: http://www.industrial-maths.com/ms6021_thomas.pdf
 * x -- input vector, function returns solution
 * n -- number of equations
 * a -- subdiagonal
 * b -- main diagonal
 * c -- superdiagonal
 */
static void trilus(double *x, double *a, double *b, double *c, int n){
  // Forward sweep
  for(int i = 1; i<n; i++){
    double m = a[i]/b[i-1];
    b[i] = b[i] - (m*c[i-1]);
    x[i] = x[i] - (m*x[i-1]);
  }
  
  x[n-1] = x[n-1]/b[n-1];
  
  // Backwards sweep
  for(int i = n-2; i >-1; i--){
    x[i] = (x[i] - (c[i]*x[i+1]))/b[i];
  }
 
}

/* Returns last index in grid of length n that is less 
   than or equal to the given value. If the value is greater 
   than all points on grid, then second to last grid point 
   is returned. */
static int bsearch(double value, double* grid, int n){
  int idx = 0;
  int high = n-1;
  int low = 0;
  int mid = low + ((high-low)/2);
  while(low<=high){
    mid = low + ((high-low)/2);
    if(grid[mid] <= value){
      idx = mid;
      low = mid+1;
    }
    else if(grid[mid] > value){
      high = mid-1;
    }
  }
    
  if(idx == n-1){
    idx = n-2;
  }
  return idx;
}


void spline_fit_natural(struct spline *spl, double* scratch){
  int n;
  double *subdiag;
  double *diag;
  double *superdiag;
  assert(spl != NULL);
  assert(scratch != NULL);
  
  n = spl->n;
  subdiag = scratch;
  diag = scratch + n;
  superdiag = diag + n;

  // 1st and nth equations from boundary conditions
  diag[0] = 1;
  superdiag[0] = 0;
  diag[n-1] = 1;
  subdiag[n-2] = 0;

  // diagonal
  for(int i = 1; i < n-1; i++){
    diag[i] = 2 * ((spl->x[i] - spl->x[i-1]) + (spl->x[i+1] - spl->x[i]));
  }
  // subdiagonal
  for(int i=0; i < n-2; i++){
    subdiag[i] = spl->x[i+1]-spl->x[i];
  }
  // superdiagonal
  for(int i = 1; i< n-1;i++){
    superdiag[i]= spl->x[i+1]-spl->x[i];
  }
  
  spl->c[0] = 0;
  spl->c[n-1] = 0;
  for(int i =1; i< n-1; i++){
    spl->c[i] = 3*(((spl->y[i+1]-spl->y[i])/(spl->x[i+1]-spl->x[i]))-
		   ((spl->y[i]-spl->y[i-1])/(spl->x[i]-spl->x[i-1])));
  }
  
  trilus(spl->c, subdiag, diag, superdiag, spl->n);
}



void spline_fit_knot(struct spline *spl, double* scratch){
  int n;
  double *subdiag;
  double *diag;
  double *superdiag;
  assert(spl != NULL);
  assert(scratch != NULL);
  
  n = spl->n;
  subdiag = scratch;
  diag = scratch + n;
  superdiag = diag + n;


  // 1st equation from boundary conditions
  diag[0] = (3 * (spl->x[1] - spl->x[0])) + (2 * (spl->x[2]-spl->x[1])) +
    (((spl->x[1]-spl->x[0]) * (spl->x[1]-spl->x[0]))/(spl->x[2]-spl->x[1]));
  superdiag[0] = (spl->x[2]-spl->x[1]) -
    (((spl->x[1]-spl->x[0]) * (spl->x[1]-spl->x[0]))/(spl->x[2]-spl->x[1]));

  // n'th equation from boundary conditions
  subdiag[n-2] = (spl->x[n-2]-spl->x[n-3]) -
    (((spl->x[n-1]-spl->x[n-2]) * (spl->x[n-1]-spl->x[n-2]))/
     (spl->x[n-2]-spl->x[n-3]));
  diag[n-1] = (3 * (spl->x[n-1] - spl->x[n-2])) +
    (2 * (spl->x[n-2]-spl->x[n-3])) +
    (((spl->x[n-1]-spl->x[n-2]) * (spl->x[n-1]-spl->x[n-2]))/
     (spl->x[n-2]-spl->x[n-3]));

  for(int i=1; i<n-1; i++){
    diag[i] = 2 * ((spl->x[i] - spl->x[i-1]) + (spl->x[i+1] - spl->x[i]));
  }
  // subdiagonal
  for(int i=0; i<n-2; i++){
    subdiag[i] = spl->x[i+1]-spl->x[i];
  }
  // superdiagonal
  for(int i = 1; i< n-1;i++){
    superdiag[i]= spl->x[i+1]-spl->x[i];
  }
  
  spl->c[0] = 3*(((spl->y[2]-spl->y[1])/(spl->x[2]-spl->x[1]))-
	    ((spl->y[1]-spl->y[0])/(spl->x[1]-spl->x[0])));
  spl->c[n-1] = 3*(((spl->y[n-1]-spl->y[n-2])/(spl->x[n-1]-spl->x[n-2]))-
	    ((spl->y[n-2]-spl->y[n-3])/(spl->x[n-2]-spl->x[n-3])));
  for(int i =1; i< n-1; i++){
    spl->c[i] = 3*(((spl->y[i+1]-spl->y[i])/(spl->x[i+1]-spl->x[i]))-
	    ((spl->y[i]-spl->y[i-1])/(spl->x[i]-spl->x[i-1])));
  }

  trilus(spl->c, subdiag, diag, superdiag, spl->n);
  
}

  
void spline_eval(struct spline *spl, double* scratch, double* eval_pts,
		 int numEval, int isfit, double* values){
  assert(values != NULL);
  assert(eval_pts != NULL);
  assert(spl != NULL);
  assert(scratch != NULL);
  if(isfit){
    spline_fit_natural(spl, scratch);
  }
  
  for(int i = 0; i<numEval; i++){
    int idx = bsearch(eval_pts[i], spl->x, spl->n)    
    double dx = eval_pts[i] - spl->x[idx];
    double b_i = ((spl->y[idx+1] - spl->y[idx])/(spl->x[idx+1] - spl->x[idx])) -
      (((spl->x[idx+1] - spl->x[idx])*(spl->c[idx+1]+(2*spl->c[idx])))/3);
    double d_i = (spl->c[idx+1] - spl->c[idx])/(3*(spl->x[idx+1] - spl->x[idx]));
    values[i] = spl->y[idx] + dx*(b_i + dx*(spl->c[idx] + (dx*d_i)));
  }
}

void spline_eval_deriv(struct spline *spl, double* scratch, double* eval_pts,
		 int numEval, int isfit, double* values){
  if(isfit){
    spline_fit_natural(spl, scratch);
  }

  for(int i = 0; i<numEval; i++){
    int idx = bsearch(eval_pts[i], spl->x, spl->n)
    double dx = eval_pts[i] - spl->x[idx];
    double c_i = spl->c[idx];
    double b_i = ((spl->y[idx+1] - spl->y[idx])/(spl->x[idx+1] - spl->x[idx])) -
      (((spl->x[idx+1] - spl->x[idx])*(spl->c[idx+1]+(2*spl->c[idx])))/3);
    double d_i = (spl->c[idx+1] - spl->c[idx])/(3*(spl->x[idx+1] - spl->x[idx]));
    values[i] = b_i + dx*((2*c_i) + (dx*3*d_i));
  }
}


void spline_eval_deriv2(struct spline *spl, double* scratch, double* eval_pts,
		 int numEval, int isfit, double* values){
  if(isfit){
    spline_fit_natural(spl, scratch);
  }

  for(int i = 0; i<numEval; i++){
    int idx = bsearch(eval_pts[i], spl->x, spl->n)
    double dx = eval_pts[i] - spl->x[idx];
    double c_i = spl->c[idx];
    double b_i = ((spl->y[idx+1] - spl->y[idx])/(spl->x[idx+1] - spl->x[idx])) -
      (((spl->x[idx+1] - spl->x[idx])*(spl->c[idx+1]+(2*spl->c[idx])))/3);
    double d_i = (spl->c[idx+1] - spl->c[idx])/(3*(spl->x[idx+1] - spl->x[idx]));
    values[i] = (2*c_i) + (dx*6*d_i);
  }
}



  
  





