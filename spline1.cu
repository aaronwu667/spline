#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <time.h>

/* To-do: put this in a header file */
struct spline{
  int grid_size;
  int nqty;
  float *x;
  float *fx;
  float *coeff;
};


// Error checking macro due to Wenjie
#define cudaCheckErrors(msg)			\
  do {						\
    cudaError_t __err = cudaGetLastError();	\
    if (__err != cudaSuccess) {				       \
      fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n",       \
	      msg, cudaGetErrorString(__err),		       \
	      __FILE__, __LINE__);			       \
      fprintf(stderr, "*** FAILED - ABORTING\n");	       \
      exit(1);						       \
    }							       \
  } while(0)

__global__ void spline_knot_system(struct spline *spl, int grid_idx, int grid_size){
  int nqty_idx = blockIdx.x;
  int n = grid_size;
  if(grid_idx == 0){
    spl->coeff[0] = 3*(((spl->fx[(n * nqty_idx) + 2]-spl->fx[(n*nqty_idx) + 1])/
			(spl->x[2]-spl->x[1]))-
		       ((spl->fx[(n*nqty_idx)+1]-spl->fx[(n*nqty_idx)])/
			(spl->x[1]-spl->x[0])));
  }
  else if(grid_idx == n - 1){
    spl->coeff[n-1] = 3*(((spl->fx[(n*nqty_idx)+n-1]-spl->fx[(n*nqty_idx)+n-2])/
			  (spl->x[n-1]-spl->x[n-2]))-
			 ((spl->fx[(n*nqty_idx)+n-2]-spl->fx[(n*nqty_idx)+n-3])/
			  (spl->x[n-2]-spl->x[n-3])));
  }
  else{
    spl->coeff[grid_idx] = 3*(((spl->fx[(n*nqty_idx)+grid_idx+1]-spl->fx[(n*nqty_idx)+grid_idx])/(spl->x[grid_idx+1]-spl->x[grid_idx]))-
			      ((spl->fx[(n*nqty_idx)+grid_idx]-spl->fx[(n*nqty_idx)+grid_idx-1])/(spl->x[grid_idx]-spl->x[grid_idx-1])));
  }
}

__global__ void spline_knot_kernel(struct spline *spl, float* subdiag,
				   float* diag, float* superdiag){
  int grid_size = spl->grid_size;
  int grid_idx = blockIdx.x;
  if(grid_idx == 0){
    subdiag[0] = 0;
    
    diag[0] = (3 * (spl->x[1] - spl->x[0])) + (2 * (spl->x[2]-spl->x[1])) +
      (((spl->x[1]-spl->x[0]) * (spl->x[1]-spl->x[0]))/(spl->x[2]-spl->x[1]));
    
    superdiag[0] = (spl->x[2]-spl->x[1]) -
      (((spl->x[1]-spl->x[0]) * (spl->x[1]-spl->x[0]))/(spl->x[2]-spl->x[1]));
  }
  else if(grid_idx == grid_size-1){
    superdiag[grid_size-1] = 0;
    
    diag[grid_size - 1] = (3 * (spl->x[grid_size-1] - spl->x[grid_size-2])) +
      (2 * (spl->x[grid_size-2]-spl->x[grid_size-3])) +
    (((spl->x[grid_size-1]-spl->x[grid_size-2]) *
      (spl->x[grid_size-1]-spl->x[grid_size-2]))/
     (spl->x[grid_size-2]-spl->x[grid_size-3]));
    
    subdiag[grid_size - 1] = (spl->x[grid_size-2]-spl->x[grid_size-3]) -
      (((spl->x[grid_size-1]-spl->x[grid_size-2]) *
      (spl->x[grid_size-1]-spl->x[grid_size-2]))/
     (spl->x[grid_size-2]-spl->x[grid_size-3]));
  }
  else{
    diag[grid_idx] =  2 * ((spl->x[grid_idx] - spl->x[grid_idx-1]) +
			   (spl->x[grid_idx+1] - spl->x[grid_idx]));
    subdiag[grid_idx] = spl->x[grid_idx]-spl->x[grid_idx-1];
    superdiag[grid_idx]= spl->x[grid_idx+1]-spl->x[grid_idx];
  }
  
  spline_knot_system<<<spl->nqty,1>>>(spl, grid_idx, grid_size);
}


__global__ void eval_kernel_nqty(struct spline *spl, float* values, float* eval_pts,
				 int pt_idx, int floor_idx, int n, int numEval){

  int nqty_idx = blockIdx.x;
  float b_i = ((spl->fx[(n*nqty_idx)+floor_idx+1] - spl->fx[(n*nqty_idx)+floor_idx])/(spl->x[floor_idx+1] - spl->x[floor_idx])) -
    (((spl->x[floor_idx+1] - spl->x[floor_idx])*(spl->coeff[(n*nqty_idx)+floor_idx+1]+(2*spl->coeff[(n*nqty_idx)+floor_idx])))/3);
  
  float d_i = (spl->coeff[(n*nqty_idx)+floor_idx+1] - spl->coeff[(n*nqty_idx)+floor_idx])/(3*(spl->x[floor_idx+1] - spl->x[floor_idx]));
  
  values[(numEval*nqty_idx) + pt_idx] = spl->fx[(n*nqty_idx)+floor_idx] + (b_i*(eval_pts[pt_idx]-spl->x[floor_idx])) +
    (spl->coeff[(n*nqty_idx) + floor_idx]*((eval_pts[pt_idx]-spl->x[floor_idx])*(eval_pts[pt_idx]-spl->x[floor_idx]))) +
    (d_i*((eval_pts[pt_idx]-spl->x[floor_idx])*(eval_pts[pt_idx]-spl->x[floor_idx])*
	  (eval_pts[pt_idx]-spl->x[floor_idx]))); 
}



__global__ void spline_eval_kernel(struct spline *spl, float* values, float* eval_pts, int numEval){
  int pt_idx = blockIdx.x;
  int floor_idx;
  int n = spl->grid_size;
  /* linear search for index */
  if(eval_pts[pt_idx] <= spl->x[0]){
    floor_idx = 0;
  }
  else if(eval_pts[pt_idx] >= spl->x[n-1]){
    floor_idx = n-2;
  }
  for(int i = 0; i<n; i++){
    if(eval_pts[pt_idx]<spl->x[i]){
      floor_idx = i-1;
      break;
    }
  }
  eval_kernel_nqty<<<spl->nqty,1>>>(spl, values, eval_pts, pt_idx, floor_idx, n, numEval);
}


void spline_knot_fit(struct spline *spl, float* coeff, float* scratch, void* buffer, int grid_size, int nqty){


  float* subdiag = scratch;
  float* diag = subdiag + grid_size;
  float* superdiag = diag + grid_size;
  
  spline_knot_kernel<<<grid_size,1>>>(spl, subdiag, diag, superdiag);

  
  cusparseStatus_t status;
  cusparseHandle_t handle = 0;
  status = cusparseCreate(&handle);
  if(status != CUSPARSE_STATUS_SUCCESS){
    fprintf(stderr, "cuSparse failed to initialize");
  }
  
  
 
  cusparseSgtsv2_nopivot(handle, grid_size, nqty, subdiag,
			 diag, superdiag, coeff, grid_size, buffer);
 
  
}

void spline_eval(struct spline *spl, float* values, float* eval_pts,
		 int numEval, float* scratch, int fitFlag){
  /*
  if(fitFlag){
    spline_knot_fit(spl, scratch);
  }
  */
  spline_eval_kernel<<<numEval,1>>>(spl, values, eval_pts, numEval);
}




int main(){
  float x[500];
  float y[1000000];
  int n = 500;
  float c[1000000];
  float scratch[5000];

   
  FILE* fx = fopen("fx.dat", "r");
  FILE* x_grid = fopen("x_grid.dat", "r");
  for(int i = 0; i<1000000; i++){
    float entry;
    fscanf(fx,"%f", &entry);
    y[i] = entry; 
  }

  for(int i =0; i<500; i++){
    float entry;
    fscanf(x_grid, "%f", &entry);
    x[i] = entry;
  }
    

  
  struct spline test;
  test.grid_size = 500;
  test.nqty = 2000;
  
  struct spline *dtest;
  float* dx;
  float* dy;
  float* cd;
  float* dscratch;

  
  cudaMalloc((void **)&dtest, sizeof(struct spline));
  cudaMalloc((void**)&cd, sizeof(float)*1000000);
  cudaMalloc((void**)&dscratch, sizeof(float) * 5000);
  cudaMalloc((void **)&dx, sizeof(float) * 500);
  cudaMalloc((void **)&dy, sizeof(float) * 1000000);
  
  cudaMemcpy(dx, x, sizeof(float) * 500,cudaMemcpyHostToDevice);
  cudaMemcpy(dy, y, sizeof(float) * 1000000,cudaMemcpyHostToDevice);
 
 
  test.x = dx;
  test.fx = dy;
  test.coeff = cd;
  cudaMemcpy(dtest, &test, sizeof(struct spline), cudaMemcpyHostToDevice);
  void* buffer;
  cudaMalloc(&buffer, test.grid_size*(3+test.nqty)*sizeof(float));
  
  clock_t begin = clock();
  spline_knot_fit(dtest, cd, dscratch, buffer, test.grid_size,test.nqty);

  clock_t end = clock();
  float time_spent = (float)(end - begin) / CLOCKS_PER_SEC;

  
  printf("Time: %f", time_spent);
  
  return 0;
}
