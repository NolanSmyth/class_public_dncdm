#include "common.h"
#include "pseudo_spectral.h"

int main(){
  double *x, *a, *clog, *s, *d, *dp, *w;
  int n, p;
  ErrorMsg error_message;

  n = 8; // Number of collocation points
  p = 2; // Number of derivatives to compute

  class_alloc(x, n*sizeof(double), error_message);
  class_alloc(d, p*n*sizeof(double), error_message);
  class_alloc(dp,p*n*n*sizeof(double), error_message);
  
  class_alloc(a, n*sizeof(double), error_message);
  class_alloc(s, n*sizeof(double), error_message);
  class_alloc(clog, n*sizeof(double), error_message);

  for (int i = 0; i < n; i++)
  {
    x[i] = 10.*((double) i)/((double) (n-1));
  }

  class_call(init(n, p, x, a, clog, s, d, dp, alpha_welfert),
            error_message,
            error_message);
  class_call(diag(n, p, x, d),
            error_message,
            error_message);
  class_call(offdiag(n, p, x, d, dp),
            error_message,
            error_message);

  int istart = (p-1) * n * n;
  int iend = p * n * n;
  printf("derivative matrix: \n");
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      printf("%4.3f    ",dp[istart + i*n + j]); 
    }
    printf("\n");
  }

  printf("Signs:\n");
  for (int i = 0; i < n; i++)
  {
    printf("%4.3f    ",s[i]);
  }
  printf("\n");
  printf("Clogs:\n");
  for (int i = 0; i < n; i++)
  {
    printf("%4.3f    ",clog[i]);
  }
  printf("\n");

  free(x);
  free(d);
  //free(dp);

  free(s);
  free(a);
  free(clog);


  n = 6; // Number of collocation points
  p = 1; // Number of derivatives to compute

  class_alloc(x, n*sizeof(double), error_message);
  class_alloc(w, n*sizeof(double), error_message);
  class_alloc(d, p*n*sizeof(double), error_message);
  class_alloc(dp,p*n*n*sizeof(double), error_message);
  
  class_alloc(a, n*sizeof(double), error_message);
  class_alloc(s, n*sizeof(double), error_message);
  class_alloc(clog, n*sizeof(double), error_message);

  printf("Legendre collocation -------------------------------------------------------\n");
  // Note that xi and w have length n+1
  // The first elements are set to 0 by hand -- this is the end point that is part of the Laguerre collocation scheme, but not part of the Gauss-Laguerre 
  // quadrature
  x[0] = 0.0;
  w[0] = 0.0;

  // Compute the Laguerre quadrature abcissas (also the interior collocation points) and weights.
  // Fill the last n entries of the arrays
  class_call(gaulag(&x[0], &w[0], n-1, 0),
            error_message,
            error_message);

  class_call(init(n, p, x, a, clog, s, d, dp, alpha_laguerre),
            error_message,
            error_message);
  class_call(diag(n, p, x, d),
            error_message,
            error_message);
  class_call(offdiag(n, p, x, d, dp),
            error_message,
            error_message);

  istart = (p-1) * n * n;
  printf("derivative matrix: \n");
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      printf("%4.3f    ", dp[istart + i*n + j]); 
    }
    printf("\n");
  }

  
  return _SUCCESS_;
}


  
