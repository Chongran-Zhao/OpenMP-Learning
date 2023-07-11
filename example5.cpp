#include <omp.h>
#include <stdio.h>

int main() {
    // Check if in parallel region
    if (omp_in_parallel()) {
        printf("Currently inside a parallel region.\n");
    } else {
        printf("Currently not inside a parallel region.\n");
    }

    // Get dynamic adjustment status
    int dynamic = omp_get_dynamic();
    printf("Dynamic adjustment of threads is currently %s.\n", dynamic ? "enabled" : "disabled");

    // Set dynamic adjustment to disabled
    omp_set_dynamic(1);
    printf("Dynamic adjustment of threads has been enabled.\n");

    // Get the number of available processors
    int num_procs = omp_get_num_procs();
    printf("Number of available processors: %d\n", num_procs);

    // Check if dynamic adjustment is still enabled
    dynamic = omp_get_dynamic();
    printf("Dynamic adjustment of threads is currently %s.\n", dynamic ? "enabled" : "disabled");

    return 0;
}

