#include "common.h"
#include "pseudo_spectral.h"

int main(){
  double x;
  double *xi, *c, *d1, *w;
  int n, m;
  ErrorMsg error_message;

  n = 5; // Number of collocation points - 1
  m = 1; // Number of derivatives to compute

  class_alloc(xi, (n+1)*sizeof(double), error_message);
  class_alloc(w, (n+1)*sizeof(double), error_message);
  class_alloc(c, (n+1)*(m+1)*sizeof(double), error_message);
  class_alloc(d1, (n+1)*(n+1)*sizeof(double), error_message);

  x = 0.5;

  printf("Chebyshev-Lobatto collocation (inludes endpoints) -------------------------------------------------------\n");
  // Chebyshev collocation pointis
  for (int i = 0; i <= n; i++){
    xi[i] = - cos(i*_PI_/n);
  }

  class_call(get_cardinal_and_diff_matrix(x, xi, c, n+1, m+1),
            error_message,
            error_message);

  //class_call(quadrature_in_rectangle(-2.0,2.0,1.0,4.0,&n,&x,&y,&w,error_message),
  //	     error_message,
  //	     error_message);

  printf("Cardinal functions evaluated at x = %e\n",x);
  for (int i = 0; i <= n; i++){
    printf("Collocation points, cardinal funcs xi = %e %e\n", xi[i], exp(-(x-xi[i])/2.)*c[i*(m+1) + 0]);
  }

  for (int i = 0; i <= n; i++)
  {
    class_call(get_cardinal_and_diff_matrix(xi[i], xi, c, n+1, m+1),
            error_message,
            error_message);
    for (int j = 0; j <= n; j++)
    {
      d1[i*(n+1) + j] = c[j*(m+1) + 1]; 
    }
  }
  printf("First derivative matrix: \n");
  for (int i = 0; i <= n; i++){
    for (int j = 0; j <= n; j++){
      printf("%4.3e    ",d1[i*(n+1) + j]); 
    }
    printf("\n");
  }

  printf("Legendre collocation -------------------------------------------------------\n");
  // Note that xi and w have length n+1
  // The first elements are set to 0 by hand -- this is the end point that is part of the Laguerre collocation scheme, but not part of the Gauss-Laguerre 
  // quadrature
  xi[0] = 0.0;
  w[0] = 0.0;

  // Compute the Laguerre quadrature abcissas (also the interior collocation points) and weights.
  // Fill the last n entries of the arrays
  class_call(gaulag(&xi[0], &w[0], n, 0),
            error_message,
            error_message);

  /*
  class_call(gaulag(xi, w, n, 1),
            error_message,
            error_message);
  */

  for (int i = 0; i < n+1; i++){
    printf("Collocation points, weights xi, w = %e %e\n", xi[i], w[i]);
  }

  class_call(get_cardinal_and_diff_matrix(x, xi, c, n+1, m+1),
            error_message,
            error_message);

  for (int i = 0; i <= n; i++){
    printf("Collocation points, cardinal funcs xi = %e %e\n", xi[i], c[i*(m+1) + 0]);
  }

  for (int i = 0; i <= n; i++)
  {
    class_call(get_cardinal_and_diff_matrix(xi[i], xi, c, n+1, m+1),
            error_message,
            error_message);
    for (int j = 0; j <= n; j++)
    {
      d1[i*(n+1) + j] = c[j*(m+1) + 1]; 
      if (i == j)
        d1[i*(n+1) + j] += -1./2.;
      d1[i*(n+1) + j] *= exp(-xi[i]/2.)/exp(-xi[j]/2.);
    }
  }

  printf("First derivative matrix: \n");
  for (int i = 0; i <= n; i++){
    for (int j = 0; j <= n; j++){
      printf("%4.3e    ",d1[i*(n+1) + j]); 
    }
    printf("\n");
  }


  free(xi);
  free(c);
  free(w);

  return _SUCCESS_;
}


  
