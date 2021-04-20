#ifndef SPLINE
#define SPLINE

/* Real spline type. Client SHOULD NOT directly access
   any of these fields. */ 
struct spline{
  int n;

  double *x;
  
  double *y;
  
  double *c;
};
  
/* Given spline with values y on grid x, fits a natural spline
   and stores fitted second derivative coefficients in c. NOTE: 
   scratch should be at least 3*sizeof(double)*n, where n
   is the number of interpolation points. */
void spline_fit_natural(struct spline *spl, double* scratch);


/* Given spline with values y on grid x, fits a not-a-knot spline
   and stores fitted second derivative coefficients in c. NOTE: 
   scratch should be at least 3*sizeof(double)*n, where n
   is the number of interpolation points. */
void spline_fit_knot(struct spline *spl, double* scratch);

/* Evalutes spline at numEval given points. Stores result in values. 
   If isFit flag is set, function will fit the spline before evaluating. 
   As before, client is responsible for allocating enough scratch
   space: 3 * sizeof(double) * n if fitting is desired. If only evaluation
   is deisred, a null pointer can be passed. */
void spline_eval(struct spline *spl, double* scratch, double* eval_pts,
		 int numEval, int isfit, double* values);

/* Evalutes spline derivative at numEval given points. Stores result in values. 
   If isFit flag is set, function will fit the spline before evaluating. 
   As before, client is responsible for allocating enough scratch
   space: 3 * sizeof(double) * n if fitting is desired. If only evaluation
   is deisred, a null pointer can be passed. */
void spline_eval_deriv(struct spline *spl, double* scratch, double* eval_pts,
		       int numEval, int isfit, double* values);

/* Evalutes spline second derivatives at numEval given points. 
   Stores result in values. If isFit flag is set, function will fit 
   the spline before evaluating. As before, client is responsible for 
   allocating enough scratch space: 3 * sizeof(double) * n if fitting 
   is desired. If only evaluation is deisred, 
   a null pointer can be passed. */
void spline_eval_deriv2(struct spline *spl, double* scratch, double* eval_pts,
			int numEval, int isfit, double* values);


#endif
