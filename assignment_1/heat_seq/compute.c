#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>	/* for uint64 definition */
//#include "img.h"
#include "compute.h"

#define BILLION 1000000000L
#define FPOPS_PER_POINT_PER_ITERATION (                 \
        1     /* current point 1 mul */ +               \
        3 + 1 /* direct neighbors 3 adds + 1 mul */ +   \
        3 + 1 /* diagonal neighbors 3 adds + 1 mul */ + \
        2     /* final add */ +                         \
        1     /* difference old/new */                  \
        )
/// Replace/erase the following line:
#include "ref1.c"

void do_compute(const struct parameters* p, struct results *r)
{
	//make txt file for results
	FILE *file;
	char size[20]; 
	sprintf(size, "%ldx%ld.txt", p->N, p->M);
	file = fopen(size, "w");
	char data[200];
	sprintf(data,"Output:\n\n"
               "%13s %13s %13s %13s %13s %13s %13s\n",
               //PACKAGE_NAME, PACKAGE_VERSION, PACKAGE_BUGREPORT,
               "Iterations",
               "T(min)", "T(max)", "T(diff)", "T(avg)", "Time", "FLOP/s");
	fputs(data, file);

	//time of computing
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);	/* mark start time */

	//initiate values
	int N = p->N;
	int M = p->M;
	printf("%d %d\n", N,M );
	int maxiter = p->maxiter;
	float T[(N+2) *M];
	//store tinit in T
	for (int j = 0; j < M; j ++) {
		T[j] = p -> tinit[j];
		T[(N + 1) * M + j] = p -> tinit[(N-1) * M + j];
	}

	for(int i = 0; i<N;i++) {
		for(int j = 0; j<M; j++){
			//T[i][j] = *(*(p->tinit + i)+j);
			T[(i + 1)*M + j] = p->tinit[i*M +j];
		}		
	}
	
	float weight;
	float joint_weight;
	float direct_w = (sqrt(2.0000)/(sqrt(2.0000)+1.0000));
	float sum_direct;
	float indirect_w = (1.0000/(sqrt(2.0000)+1.0000));
	float sum_indirect;
	float min_weight;
	float old_temp;
	float new_temp;
	int k = 0;
	float Tmin;
	float Tmax;
	float maxdiff;
	float sum_temp;
	float new_T[N*M];
	//threshold condition
	bool thresh_cond = false;
	//printreports
	int report = p->period;

	//while (k < maxiter or diff between old and new (mean) temp < p->threshold)
	while (k < maxiter && !thresh_cond) {
		//begin_picture (k, M , N, p->io_tmin , p->io_tmax );
		Tmin = p->io_tmax;
		Tmax = p->io_tmin;
		maxdiff = 0;
		sum_temp = 0;
		thresh_cond = true;
		for(int i = 1; i<(N+1);i++) {
			for(int j = 0; j<M; j++){
				// calculate new temperature
				weight = p->conductivity[(i-1)*M + j];
				min_weight = 1 - weight;
				old_temp = T[i * M + j];
				sum_direct = (T[(i-1)*M +j] + T[(i+1)*M + j] + T[i*M + (((j-1)+M)%M)] + T[i*M + ((j+1)%M)])/4.0000;
				sum_direct = sum_direct * direct_w;
				sum_indirect = (T[(i-1)*M +(((j-1)+M)%M)] + T[(i-1)*M + ((j+1)%M)] + T[(i+1)*M + (((j-1)+M)%M)] + T[(i+1)*M + ((j+1)%M)])/4.0000;
				sum_indirect = sum_indirect * indirect_w;
				new_temp = (old_temp * weight) + (min_weight * (sum_direct + sum_indirect));
				new_T[(i-1)*M + j] = new_temp; 
				//printf("%f\t", old_temp);
				if (fabsf(old_temp - new_temp) > p->threshold ) {
					thresh_cond = false;
				}
				//update results
				//Tmin
				if (new_temp < Tmin){
					Tmin = new_temp;

					//printf("min %f %d %d \n", Tmin, i , j);
				}
				//Tmax
				if (new_temp > Tmax){
					Tmax = new_temp;
					//printf("max %f %d %d \n", Tmax, i , j);
				}
				//max difference
				if (fabsf(old_temp - new_temp) > maxdiff) {
					maxdiff = fabsf(old_temp - new_temp);
				}
				//avg temp
				sum_temp += new_temp;
			}	
		}
		r->niter = k;
		r->tmin = Tmin;
		r->tmax = Tmax;
		r->tavg = sum_temp/(N*M);
		r->maxdiff = maxdiff;
		for(int i = 0; i<N;i++) {
			for(int j = 0; j<M; j++){
				T[(i+1)*M +j] = new_T[i*M +j];
				//draw_point ( j , i, T[(i + 1)*M + j]);
			}
		}
		
		
		//end_picture ();

		//update k
		
		clock_gettime(CLOCK_MONOTONIC, &end);	/* mark the end time */
		r->time = BILLION * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
		if(k%report == 0 && p->printreports == 1){
			sprintf(data, "%-13zu % .6e % .6e % .6e % .6e % .6e % .6e\n",
           r->niter,
           r->tmin,
           r->tmax,
           r->maxdiff,
           r->tavg,
           r->time,
           (double)p->N * (double)p->M * 
           (double)(r->niter * FPOPS_PER_POINT_PER_ITERATION +
                    (double)r->niter / p->period) / r->time);
    		fputs(data,file);
			//report results
			report_results(p, r);
			
		}
		k++;
		
	}
	fclose(file);
}
