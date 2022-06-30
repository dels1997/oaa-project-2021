#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

//constant C corresponds to the cost of adding a new line
//extreme cases correspond to C=0.0 in which all the errors are zero since we can choose lines through each two points with no error
//and also to C greater than e.g. the sum of all errors so the optimal solution is one line through all the points
#define C 0.1

//we use the structure point consisting of two doubles since this will allow us to apply the built-in qsort directly
typedef struct {
	double x;
	double y;
} point;

//compare function for quicksort to sort in ascending order of x coordinates
int cmpfun (const void* a, const void* b) {
	if (((point*)a)->x>((point*)b)->x) return 1;
	return -1;	
}

//function to find the e_{i,j} matrix of least square errors for all the combinations of points
void fill2D (double** e, double** a, double ** b, int n, const point* p) {
	int i, j, k;
	double Sx, Sy, Sxy, Sx2, Sy2;
	double atemp, btemp;

	for(i=0; i<n; i++) {
		Sx=Sy=Sxy=Sx2=Sy2=0;
		k=0;
		for(j=i; j<n; j++) {
			
			Sx+=(p+j)->x;
			Sy+=(p+j)->y;
			Sxy+=(p+j)->x*(p+j)->y;
			Sx2+=(p+j)->x*(p+j)->x;
			Sy2+=(p+j)->y*(p+j)->y;
						
			atemp=((j-i+1)*Sxy-Sx*Sy)/((j-i+1)*Sx2-Sx*Sx);
			btemp=(Sy-atemp*Sx)/(j-i+1);
			
			if(k==0) {
				a[i][k]=0.0;
				b[i][k]=0.0;
			}
			else {
				a[i][k]=atemp;
				b[i][k]=btemp;
			}
						
			if(j==i || j==i+1) {
				e[i][k++]=0.0;
				continue;
			}
			
			e[i][k++]=Sy2+atemp*atemp*Sx2+btemp*btemp*(j-i+1)-2*atemp*Sxy-2*btemp*Sy+2*atemp*btemp*Sx;
		}
	}
}

//function to print the matrix of least square errors for all the combinations of points
void print2D (double** e, double** a, double** b, int n) {
	int i, j;
	
	//matrix of least square errors
	for(i=0; i<n; i++) {
		printf("|");
		for(j=0; j<n-i; j++) {
			printf(" %lf |", e[i][j]);
		}
		printf("\n");
	}
	
	//matrix of parameters of best line fits in format (a,b) for y=ax+b
	for(i=0; i<n; i++) {
		printf("|");
		for(j=0; j<n-i; j++) {
			printf(" (%9.6lf,%9.6lf) |", a[i][j], b[i][j]);
		}
		printf("\n");
	}
}

//function to find exact partition of points corresponding to the optimal solution
void findSeg (int i, double** e, double** a, double** b, double M[]) {
	if(i==0)
		return;
		
	int j, flag=1;
	double temp, min;
	
	for(j=1; j<=i; j++) {
		temp=e[j-1][i-j]+C+M[j-1];
		if(j==1)
			min=temp;
		if(temp<min) {
			min=temp;
			flag=j;
		}
	}
	printf("(%d,%d)\n", flag, i);
	printf("(%lf,%lf)\n", a[flag-1][i-flag], b[flag-1][i-flag]);
	findSeg(flag-1, e, a, b, M);
}

int main (void) {
	double min;
	int i, j, n;
	
	printf("How many points do you wish to apply the algorithm to?\nNumber of points=");
	scanf(" %d", &n);
	
	point* p=(point*)malloc(n*sizeof(point));
	
	printf("Input their coordinates (indices range from 1 to %d):\n", n);
	for(i=0; i<n; i++) {
		printf("x[%d]=", i+1);
		scanf(" %lf", &(p+i)->x); printf("\n");
		//we can keep the following line to force the user to input the sorted set of points,
		//however for large n this may be intractable
		//if(i) if((p+i)->x<(p+i-1)->x) {printf("Ascending order of x coordinates!\n"); exit(0);} 
		printf("y[%d]=", i+1);
		scanf(" %lf", &(p+i)->y); printf("\n");
	}	
	
	//start the stopwatch for evaluation time
	//struct timeval  tv1, tv2;
	//gettimeofday(&tv1, NULL);
	
	//we sort the points w.r.t. their x coordinates to prepare them for algorithm evaluation
	qsort(p, n, sizeof(point), cmpfun);	
	

	//dynamic memory allocation for 2D array in which the minimal errors of lines through different contiguous subsets of points will be stored
	double** e=(double**)malloc((n)*sizeof(double*));
	double** a=(double**)malloc((n)*sizeof(double*));
	double** b=(double**)malloc((n)*sizeof(double*));
	
	for(i=0; i<n; i++) {
		e[i]=(double*)malloc((n-i)*sizeof(double));
		a[i]=(double*)malloc((n-i)*sizeof(double));
		b[i]=(double*)malloc((n-i)*sizeof(double));
	}
	
	//the following line invokes the function which calculates the minimal quadratic errors for all combinations of starting and ending points
	fill2D(e, a, b, n, p);
	
	//with the following line, the table of relevant e values can be printed, note that the jth input in the ith row corresponds
	//to error of best linear fit for points from ith one to the one that is j places distant (both included)
	//print2D(e, a, b, n);
	
	//array M for storing optimal solutions of subproblems has n+1 elements since we use the initial condition M[0]=0.0 when there are no points
	double *M;
	M=(double*)malloc((n+1)*sizeof(double));
	
	M[0]=0.0;
	min=0.0;
	
	double temp;
	
	for(i=1; i<=n; i++) {
		M[i]=e[0][i-1]+C+M[0];
		min=e[0][i-1]+C+M[0];
						
		for(j=2; j<=i; j++) {
			temp=e[j-1][i-j]+C+M[j-1];
			if(temp<min)
				M[i]=min=temp;
		}		
	}
	
	//printing out the optimal results for all subproblems
	printf("The optimal solutions M[i] for the subproblems consisting of i points are:\n");
	for(i=0; i<=n; i++)
		printf("\nM[%d]=%lf\n", i, M[i]);
	
	//we now print the segments in the optimal partitioning
	printf("Segments in the optimal partitioning are (in format: (firstpoint, lastpoint)newline(a,b) for y=ax+b linear fit; indices begin at 1):\n");
	findSeg(n, e, a, b, M);
	
	//freeing the allocated memory
	for(i=0; i<n; i++) {
		free(e[i]);
		free(a[i]);
		free(b[i]);
	}
		
	free(e);
	free(a);
	free(b);
	free(p);

	//gettimeofday(&tv2, NULL);
	//the previous and the following line are used to calculate the time needed for execution
	//printf("\nTotal time = %f seconds\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec));
		
	return 0;	
}
