#include<stdio.h>
#include"omp.h"
#define NUM_THREADS 1
static long num_steps = 1000000;
double dx;

int main()
{
  int i, num_threads;
  long double pi = 0.0;
  long double sum;
  long double start_time, end_time;
  dx = 1.0 / (long double) num_steps;
  omp_set_num_threads(NUM_THREADS);
  start_time = omp_get_wtime();
#pragma omp parallel
  {
    int id;
    long double x;
    id = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    for (i = id, sum = 0.0; i < num_steps; i = i+num_threads)  
    {
      x = (i + 0.5) * dx;
      sum += 4.0 / (1.0 + x*x);
    }
    #pragma omp critical
    pi += sum * dx;
  }
  end_time = omp_get_wtime();
  printf("The numerical integration of pi is %.10Lf.\n", pi);
  printf("Execution time: %Lf seconds.", end_time - start_time);
  return 0;
} 
