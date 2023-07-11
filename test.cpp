#include <omp.h>
#include <stdio.h>

int main() {
  if (omp_in_parallel()) {
    printf("Thread is executing in a parallel region.\n");
  } else {
    printf("Thread is executing in a serial region.\n");
  }
  #pragma omp parallel
  {
    if (omp_in_parallel()) {
      printf("Thread is executing in a parallel region.\n");
    } else {
      printf("Thread is executing in a serial region.\n");
    }
  }
  return 0;
}
