#include <stdio.h>
#include <stdlib.h>

int main (void) {
	int i, n;
	printf("How many points?\nNumber of points=");
	scanf(" %d", &n);
	
	double* x=(double*)malloc(n*sizeof(double));
	double* y=(double*)malloc(n*sizeof(double));
	
	printf("Input their coordinates (indices range from 1 to %d):\n", n);
	for(i=0; i<n; i++) {
		printf("x[%d]=", i+1); scanf(" %lf", &x[i]); printf("\n");
		printf("y[%d]=", i+1); scanf(" %lf", &y[i]); printf("\n");
	}
	double Sx, Sy, Sxy, Sx2, Sy2;
	Sx=Sy=Sxy=Sx2=Sy2=0;
	for(i=0; i<n; i++) {
		Sx+=x[i];
		Sy+=y[i];
		Sxy+=x[i]*y[i];
		Sx2+=x[i]*x[i];
		Sy2+=y[i]*y[i];
	}
	double a=(n*Sxy-Sx*Sy)/(n*Sx2-Sx*Sx);
	double b=(Sy-a*Sx)/n;
	
	double e=0, eaux;
	for(i=0; i<n; i++) {
		eaux=y[i]-a*x[i]-b;
		e+=eaux*eaux;	
	}
	
	free(x);
	free(y);
	
	printf("best linear fit parameters are a=%lf and b=%lf\n", a, b);
	printf("total square error is: %lf\n", e);
	
	return 0;
}
