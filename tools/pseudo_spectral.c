#include "pseudo_spectral.h"


/* 
* Compute the differentiation matrix for pseudospectral collocation. Based on weights.h in C++ NR, which, in turn 
* is based on the pseudo-code in ``CALCULATION OF WEIGHTS IN FINITE DIFFERENCE FORMULAS'' by Fornberg (1998).
* x is the location at which the differentiation matrices are to be evaluated
* xi are the collocation points
* Np is the number of collocation points
* Nd is the number of derivatives to compute >= 0
* c is the output matrix
* This only gives the correct answer if the ``weight'' that goes into the polynomial approximation is 1.
* For Laguerre collocation, the differentiation matrix must be modified.
*/
int get_cardinal_and_diff_matrix(double x, double *xi, double *c, int Np, int Nd)
{
int n = Np - 1;
int m = Nd - 1;
double c1=1.0;
double c4=xi[0]-x;

for (int j=0;j<=n;j++)
 for (int k=0;k<=m;k++)
  c[j*Nd + k]=0.0;

c[0*Nd + 0]=1.0;

for (int i=1;i<=n;i++) {
 int mn = m < i ? m : i;
 double c2=1.0;
 double c5=c4;
 c4=xi[i]-x;

 for (int j=0;j<i;j++) {
  double c3=xi[i]-xi[j];
  c2=c2*c3;

  if (j == i-1) {
   for (int k=mn;k>0;k--)
    c[i*Nd + k]=c1*(k*c[(i-1)*Nd + k-1]-c5*c[(i-1)*Nd + k])/c2;

    c[i*Nd + 0]=-c1*c5*c[(i-1)*Nd + 0]/c2;
  }

  for (int k=mn;k>0;k--)
   c[j*Nd + k]=(c4*c[j*Nd + k]-k*c[j*Nd + k-1])/c3;

  c[j*Nd + 0]=c4*c[j*Nd + 0]/c3;
 }
 c1=c2;
}

return _SUCCESS_;
}

// Compute n Gauss-Laguerre  quadrature abcissas x and weights w
// alf = 0 for standard Laguerre polynomials, and alf > 0 for their generalized version (which are orthogona wrt exp(-x) x^alf)
// From NR; modified to return exp(x_i) * w_i, instead of x_i itself.
int gaulag(double *x, double *w, int n, double alf)
{
  double EPS = 3.0e-14;
  int MAXIT = 10;

  //double gammln(double xx);
  //void nrerror(char error_text[]);
  int i,its,j;
  double ai;
  double p1,p2,p3,pp,z,z1;

  for (i=1;i<=n;i++) {
    if (i == 1) {
      z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
    } else if (i == 2) {
      z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
    } else {
      ai=i-2;
      z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
        (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
    }
    for (its=1;its<=MAXIT;its++) {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {
        p3=p2;
        p2=p1;
        p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
      }
      pp=(n*p1-(n+alf)*p2)/z;
      z1=z;
      z=z1-p1/pp;
      if (fabs(z-z1) <= EPS) break;
    }
    //if (its > MAXIT) nrerror("too many iterations in gaulag");
    x[i]=z;
    w[i] = -exp(x[i] + lgamma(alf+n)-lgamma((double)n))/(pp*n*p2);
  }

  return _SUCCESS_; 
}

// The routines below are translations from 
// GENERATION OF PSEUDOSPECTRAL DIFFERENTIATION MATRICES by Bruno Welfert (1997)
// SIAM J. NUMER. ANAL. Vol. 34, No. 4, pp. 1640–1657 (1997)
// They return spectral differentiation matrices of arbitrary order.
int init(int n, int p, double *x, double *a, double *clog, double *s, double *d, double *dp, double (*alpha)(int,double))
{
  double ci, si, xi, xk, xik;
  
  // log of the alpha function
  for (int i = 0; i < n; i++){
    a[i] = alpha(0,x[i]);

    printf("x[i], a[i] = %e \t %e\n",x[i], a[i]);
  }
  // scaling matrix
  for (int i = 0; i < n; i++)
  {
    ci = 0.0;
    si = 1.0;
    xi = x[i];
    for (int k = 0; k < n; k++)
    {
      if (k == i)
        continue;
      xik = xi - x[k];
      if (xik < 0.)
        si = -si;
      ci = ci + log(fabs(xik));
    }
    clog[i] = ci + a[i];
    s[i] = si;
  }

  //initialize the derivative matrix
  for (int k = 0; k < n; k++)
  {
    xk = x[k];
    for (int m = 1; m <= p; m++)
    {
      d[k + (m - 1) * n] = alpha(m, xk);
    } 
  }
  return _SUCCESS_;
}

int diag(int n, int p, double *x, double *d)
{
  double xk, xkj;
  int j, k1;

  for (int k = 0; k < n; k++)
  {
    xk = x[k];
    j = k + 1;
    if (j >= n)
      j = 0;

    while (j != k){
      //printf("k, j = %d, %d\n",k, j);
      xkj = 1./(xk - x[j]);
      for (int m = p; m >= 2; m--)
      {
        k1 = k + (m - 1) * n;
        d[k1] = d[k1] + ((double) m) * xkj * d[k1-n];
      }
      d[k] = d[k] + xkj;
      j = j % n + 1;
      if (j >= n)
        j = 0;
    }
  }
  return _SUCCESS_;
}

int offdiag(int n, int p, double *x, double *d, double *dp)
{
  int i, j1, j2, n2, k1, k2, m1, m2;
  double xj, fm;
  for (int j = 0; j < n; j++)
  {
    j1 = j;//-1;
    xj = x[j];
    for (int k = 0; k < n; k++)
    {
      if (k == j)
        continue;
      i = k + j1*n;
      dp[i] = 1.0/(x[k] - xj);
    }
    dp[j + j1*n] = d[j];
  }
  if (p == 1)
    return _SUCCESS_;
  n2 = n*n;
  for (int m = 2; m <= p; m++)
  {
    m1 = (m-2)*n;
    m2 = m1 + n;
    fm = (double) m;
    for (int j = 0; j < n; j++)
    {
      //j1 = (m1 + j - 1) * n;
      j1 = (m1 + j) * n;
      j2 = j1 + n2;
      xj = x[j];
      for (int k = 0; k < n; k++)
      {
        if (k == j)
          continue;
        k1 = j1 + k;
        k2 = j2 + k;
        dp[k2] = fm * (d[k + m1] - dp[k1]) / (x[k] - xj);
      }
      dp[j + j2] = d[j + m2];
    } 
  }
  return _SUCCESS_;
}

double alpha_welfert(int m, double x)
{
  int n = 8;
  double ai;

  if (m == 0)
    return -((double) n) * log(1.0 + x); 
  else{
    ai = 1.;
    for (int i = n; i <= n + m -1; i++)
    {
      ai = -ai*((double) i) / (1.0 + x);
    }
    return ai;
  }
  
}

double alpha_laguerre(int m, double x)
{
  double ai;

  if (m == 0)
    return -x/2.; 
  else{
    return pow(-0.5,m)*exp(-x/2.);
  }
}
