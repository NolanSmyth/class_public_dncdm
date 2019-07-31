#include <math.h>

double linear_interp(double xa[], double ya[], int n, double x)
{
	int i, ihi, ilo;
    double h,a,b;

    // Bracket x by table entries 
    ilo = 0;
    ihi = n-1;
    if (x < xa[0]) return ya[0];
    if (x > xa[n-1]) return ya[n-1];

    while (ihi-ilo > 1){
        i = (ihi + ilo) >> 1;
        if (xa[i] > x) ihi = i;
        else ilo = i;
    }
    h = xa[ihi] - xa[ilo];	
    if (h == 0.0) return ya[ihi];

    a = (ya[ihi] - ya[ilo])/h;
    b = (xa[ihi]*ya[ilo] - xa[ilo]*ya[ihi])/h;
    return a*x + b;
}

