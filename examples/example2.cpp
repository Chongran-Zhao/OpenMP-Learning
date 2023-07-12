#include<stdio.h>
#include"omp.h"
#define NUM_THREADS 1
static long num_steps = 1000000;
double dx;

int main()
{
  int i, nthreads;
  long double pi, sum[NUM_THREADS];
  long double start_time, end_time;
  dx = 1.0 / (long double) num_steps;
  omp_set_num_threads(NUM_THREADS);
  start_time = omp_get_wtime();
#pragma omp parallel
  {
    int i, id;
    long double x;
    id = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    for (i = id, sum[id] = 0.0; i < num_steps; i = i+nthreads)
    {
      x = (i + 0.5) * dx;
      sum[id] += 4.0 / (1.0 + x*x);
    }
  }
  for (i = 0,pi = 0.0; i < NUM_THREADS; i++) pi += sum[i] * dx;
  end_time = omp_get_wtime();
  printf("The numerical integration of pi is %.10Lf.\n", pi);
  printf("Execution time: %Lf seconds.", end_time - start_time);
  return 0;
} 
