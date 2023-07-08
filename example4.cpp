#include<stdio.h>
#include"omp.h"
#define NUM_THREADS 1
static long num_steps = 1000000;
double dx;

int main()
{
  int i;
  long double pi = 0.0;
  long double sum = 0.0;
  long double start_time, end_time;
  dx = 1.0 / (long double) num_steps;
  omp_set_num_threads(NUM_THREADS);
  start_time = omp_get_wtime();
  long double x;
#pragma omp parallel for reduction(+:sum)
  for (i = 0; i < num_steps; i++)
  {
    x = (i + 0.5) * dx;
    sum += 4.0 / (1.0 + x*x);
  }
  pi = sum * dx;
  end_time = omp_get_wtime();
  printf("The numerical integration of pi is %.10Lf.\n", pi);
  printf("Execution time: %Lf seconds.", end_time - start_time);
  return 0;
} 
