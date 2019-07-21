#pragma once

#include "common.h"

int get_cardinal_and_diff_matrix(double x, double *xi, double *c, int n, int m);

int gaulag(double *x, double *w, int n, double alf);

int init(int n, int p, double *x, double *a, double *clog, double *s, double *d, double *dp, double (*alpha)(int,double));
int diag(int n, int p, double *x, double *d);
int offdiag(int n, int p, double *x, double *d, double *dp);
double alpha_welfert(int m, double x);
double alpha_laguerre(int m, double x);
