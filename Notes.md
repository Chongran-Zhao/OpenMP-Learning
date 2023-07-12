# Shared Address Based Architecture: SMP & NUMA

## Symmetric Multiprocessing
In SMP systems, multiple processors are connected to a single shared memory system. All processors have equal access to the shared memory, which means that any processor can access any memory location without any significant difference in memory access time. This symmetric access pattern is the key characteristic of SMP. Some other key points are as following:
- SMP systems typically have a small number of processors, ranging from a few to a couple of dozen.
- The processors in SMP are usually identical in terms of architecture and capabilities.
- SMP systems operate under a single operating system, and the workload is evenly distributed among the processors.
- The shared memory allows processes or threads running on different processors to communicate with each other efficiently.
- SMP systems are designed to provide high throughput and parallel processing capabilities, which are suitable for a wide range of applications.
## Non-Uniform Memory Access
NUMA is an architecture where multiple processors are connected to a system that consists of several memory modules. In a NUMA system, each processor has its own local memory (referred to as local memory nodes), and it can also access the remote memory nodes attached to other processors over a shared interconnect. Some key points are as following:
- NUMA systems can have a larger number of processors compared to SMP systems, ranging from several dozen to hundreds or even thousands.
- Processors in NUMA systems can have different capabilities and may not be identical.
- NUMA systems are designed to optimize memory access by reducing latency. Local memory access is faster than accessing remote memory.
- Each processor in a NUMA system has faster access to its own local memory and slower access to remote memory.
- NUMA systems are commonly used in large-scale servers or high-performance computing environments where memory access latency can significantly impact performance.
## Differences between SMP and NUMA
1.  Memory Access: In SMP, all processors have equal and symmetric access to the shared memory. In NUMA, each processor has its own local memory, and accessing local memory is faster than accessing remote memory.
    
2.  Scalability: SMP systems typically have a smaller number of processors, while NUMA systems can scale to a larger number of processors.
    
3.  Processor Similarity: In SMP, the processors are generally identical in terms of architecture and capabilities. In NUMA, processors can have different capabilities and may not be identical.
    
4.  Memory Latency: NUMA systems are designed to reduce memory access latency by providing faster access to local memory. SMP systems do not differentiate memory access based on locality.
    
5.  Workload Distribution: In SMP, the workload is evenly distributed among all processors. In NUMA, workload distribution can be optimized to minimize remote memory accesses and reduce latency.

## Others
There are other architectural models where processors do not share a common address space. Here are a few examples:

1.  Distributed Memory Architecture: In this model, each processor has its own private memory, and there is no shared memory accessible by all processors. To communicate or share data, explicit message passing is required between processors. Examples of distributed memory architectures include clusters and massively parallel processing (MPP) systems.
    
2.  Hybrid Architectures: Some systems combine shared memory and distributed memory models. These hybrid architectures have a combination of shared-memory nodes and distributed-memory nodes. The shared-memory nodes typically follow SMP or NUMA architectures, while the distributed-memory nodes communicate through message passing. Examples include systems like heterogeneous clusters or hybrid supercomputers.
    
3.  SIMD (Single Instruction, Multiple Data): SIMD architectures execute the same instruction on multiple data elements simultaneously. They often feature a single control unit (processor) that broadcasts instructions to multiple processing units (vector processors) that operate on different data elements. In SIMD architectures, the processors do not share a common address space.
    
4.  MIMD (Multiple Instruction, Multiple Data): MIMD architectures encompass systems where multiple processors independently execute different instructions on different data sets. Examples of MIMD architectures include clusters, grid computing, and multi-computer systems. In these architectures, processors typically have their own memory and do not share a common address space.
# OpenMP Overview
- OpenMP is a multi-threading, shared address model.
- Threads communicate by sharing variables.
- Unintended sharing of data causes ***race conditions***.
- Race Condition: when the program’s outcome changes as the threads are scheduled differently.
- To control race conditions, use ***synchronization*** to protect data conflicts, which is expensive.
## Example 1
```c#
#include<stdio.h>
#include"omp.h"
int main()
{
    #pragma omp parallel
    {
      int ID = omp_get_thread_num();
      printf("hello(%d)", ID);
      printf("world(%d)\n", ID);
    }
} 
```
**Notes**: I didn't use ```omp_set_num_threads(num_threads)``` to set the number of threads *OpenMP* will use. By default, OpenMP will use the number of threads specified by the environment variable ```OMP_NUM_THREADS``` or, if not set, it will use the default number of threads provided by the system.

# Fork-join Parallelism

<img src="./figures/IMG_1.png" alt="IMG_1" style="zoom:45%;" />
Fork-join parallelism is a programming model and execution pattern that allows for the efficient execution of parallel tasks. It consists of two main phases: the "fork" phase and the "join" phase.

1.  Fork Phase: In the fork phase, a task or a computation is divided into smaller subtasks, creating a parallel execution hierarchy. Each subtask is assigned to a separate thread or worker, and these threads execute their respective subtasks concurrently. This phase is called "fork" because the parent task spawns multiple child tasks, creating a parallel execution flow.
2.  Join Phase: In the join phase, the threads or workers wait for their subtasks to complete. Once all the subtasks have finished executing, the threads synchronize and join together, consolidating the results of their computations. This synchronization point ensures that all subtasks have completed before proceeding further. This phase is called "join" because the parallel execution flow is consolidated back into a single flow.

The fork-join parallelism model is often implemented using parallel programming frameworks or APIs such as OpenMP, Java's Fork/Join framework, or the Cilk programming language.

## Example2

$$
\int_0^1 \frac{4}{1+x^2}dx=4\arctan x\big|_0^1=\pi\approx 3.14159=\sum_{i=0}^n\frac{4}{1+x^2}\Delta x,\quad x=(i+0.5)\Delta x,\quad \Delta x=\frac{1}{n}
$$

```c#
#include<stdio.h>
#include"omp.h"
#define NUM_THREADS 4
static long num_steps = 10000000;
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
```
**Remark**: You need to manually assign each thread the number that needs to be executed in the loop, based on the index of threads. The ```long``` type  variable is used to output a more accurate execution time.

**Summary**: 
- Use ```include![ IMG_2](/Users/chongran/OpenMP-Learning/figures/ IMG_2.png)<omp.h> ``` to insert the library of *OpenMP*.
- Use ```omp_set_num_threads(NUM_THREADS);``` to set the number of threads.
- Use ```omp_get_threads_num();``` to get the *id* of threads.
- Use ```omp_get_num_threads();``` to get the *number* of threads.

# False Sharing

False sharing is a phenomenon that occurs in parallel programming when multiple threads or processors inadvertently share the same cache line, resulting in performance degradation. It is a performance issue rather than an actual sharing of data.

In modern computer architectures, the memory system is typically organized into cache lines. A cache line is a fixed-size block of memory (e.g., 64 bytes) that is loaded from main memory into the cache. When a processor accesses a memory location, it brings the entire cache line containing that location into its local cache.

False sharing arises when multiple threads or processors access different variables that happen to reside on the same cache line. Even though these variables are logically distinct and unrelated, the sharing of the cache line causes unnecessary cache invalidations and coherence traffic between the processors, degrading performance.

False sharing can significantly impact performance, as it introduces cache contention and increases memory access latency. It is particularly problematic in parallel applications where multiple threads or processors frequently access shared data.

Mitigating false sharing typically involves techniques such as:

- Padding: Inserting additional padding or dummy variables to separate variables that are prone to false sharing, ensuring they reside on different cache lines.
- Thread/Processor Affinity: Assigning threads or processors to specific cores or sockets, reducing the chances of false sharing due to cache line conflicts.
- Data Replication: Making copies of data to ensure each thread or processor works on its private copy, eliminating sharing and false sharing altogether.
- Compiler and Language Optimizations: Compiler optimizations and programming techniques, such as thread-local storage or data alignment directives, can help mitigate false sharing.

Detecting false sharing requires careful performance profiling, monitoring cache behavior, and examining cache coherence traffic. Specialized profiling tools and performance counters provided by the system or development environments can assist in identifying false sharing issues.

In example2, I also use `omp_get_time()` to output the time of running codes in order to reveal the fact of false sharing, here is the results of using different number of threads:


| Number of Threads  | Consuming Time (s) (Macbook Air)|Consuming Time (s) (Macbook Pro)|
| :---: | :---: | :---: |
|  1   |   0.025067   | 0.019296|
|   2   |   0.013962   | 0.010428|
|   4  |   0.007185  | 0.005739|
|8     |    0.005292    |0.003425|
|16 | 0.004983|0.003287|

The time consuming can’t decrease linearly as the number of threads gets larger.

## Example3 (way to solve *False Sharing*)
```c#
#include<stdio.h>
#include"omp.h"
#define NUM_THREADS 16
#define PAD 8 // assume 64 byte L1 cache line size
static long num_steps = 10000000;
double dx;

int main()
{
  int i, nthreads;
  long double pi, sum[NUM_THREADS][PAD];
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
    for (i = id, sum[id][0] = 0.0; i < num_steps; i = i+nthreads)
    {
      x = (i + 0.5) * dx;
      sum[id][0] += 4.0 / (1.0 + x*x);
    }
  }
  for (i = 0,pi = 0.0; i < NUM_THREADS; i++) pi += sum[i][0] * dx;
  end_time = omp_get_wtime();
  printf("The numerical integration of pi is %.10Lf.\n", pi);
  printf("Execution time: %Lf seconds.", end_time - start_time);
  return 0;
} 
```
**Notes**: Use ```#define PAD 8``` to avoid the conflicts caused by different threads reading the same block of memory.
**Remark**: The method is ugly.

MacBook Pro CPU information:
	*machdep.cpu.cores_per_package: 12
	machdep.cpu.core_count: 12
	machdep.cpu.logical_per_package: 12
	machdep.cpu.thread_count: 12
	machdep.cpu.brand_string: Apple M2 Pro*

| Number of Threads  | Consuming Time (s) (Example3)|Consuming Time (s) (Example2)|
| :---: | :---: | :---: |
|  1   |   0.019267   | 0.019296|
|   2   |   0.010408   | 0.010428|
|   4  |   0.005772  | 0.005739|
|8     |    0.003344    |0.003425|
|16 | 0.003456|0.003287|

**Comment**: The example is too simple to demonstrate the effectiveness of CPU's advanced performance in resolving false sharing issues. Another reason is that , the device nowadays has large memory, and *Fasle Sharing* hardly happens when processing some simple data.

# Synchronization

## Barrier Synchronization
Each thread wait at the *barrier* until all threads arrive (figure below).

```c#
  #pragma omp parallel
  {
  	int id = omp_get_thread_num();
    A[id] = big_cal1(id); // big_cal1 is a function
    #pragma omp barrier
    B[id] = big_cal2(id,A);
  }
```
​                                                          <img src="./figures/ IMG_2.png" alt=" IMG_2" style="zoom:60%;" />
In OpenMP, both **implicit** and **explicit barriers** are synchronization mechanisms used to ensure that all threads reach a particular point in the code before proceeding further. The difference between them lies in how they are triggered and their visibility to the programmer.

- Implicit Barrier:
  - An implicit barrier is automatically inserted at the end of certain OpenMP constructs, such as `paralle`l, `for`, and `sections`.
  - The implicit barrier ensures that all threads participating in the parallel region or loop have completed their assigned work before proceeding to the next sequential code outside the construct.
  - The implicit barrier is added by the OpenMP runtime system and is not explicitly specified by the programmer.
  - The implicit barrier is a fundamental feature of OpenMP that simplifies parallel programming by ensuring that all threads synchronize implicitly at specific points in the program.
- Explicit Barrier:
    - An explicit barrier is a synchronization point that can be explicitly specified by the programmer using the ``#pragma omp barrier` directive.
    - The explicit barrier ensures that all threads reach the barrier point before any thread is allowed to proceed further.
    - Unlike the implicit barrier, the explicit barrier is explicitly added by the programmer in the code.
      The explicit barrier provides finer control over synchronization points and allows the programmer to enforce synchronization at specific locations in the code where it is necessary.
    - Explicit barriers are especially useful when there is a need to coordinate data dependencies or when multiple threads need to synchronize their work at a particular point.
    - In summary, the main difference between implicit and explicit barriers in OpenMP is that implicit barriers are automatically inserted by the OpenMP runtime system at the end of certain constructs, while explicit barriers are explicitly added by the programmer using the #pragma omp barrier
    - directive. Implicit barriers simplify synchronization in most cases, while explicit barriers provide finer control over synchronization points when needed.
### More detail about implicit barrier
Let's dive into more detail about implicit barriers in OpenMP and their placement in specific constructs.

In OpenMP, an implicit barrier is automatically inserted at the end of certain constructs to ensure synchronization among participating threads. These constructs include:

- `parallel` Construct:
  The parallel construct creates a team of threads that execute the code block in parallel.
  At the end of the parallel region, an implicit barrier is inserted.
  The barrier ensures that all threads in the team complete their parallel work before moving to the sequential part of the code outside the parallel region.
  Example:

	```c#
	#pragma omp parallel
	{
	    // Parallel work
	    // ...
	    // Implicit barrier at the end of the parallel region
	}
	```
	
- `for` Construct:
  The for construct distributes loop iterations among the available threads for parallel execution.
  At the end of the for construct, an implicit barrier is inserted.
  The barrier ensures that all threads finish their assigned iterations before proceeding to the next sequential code outside the for loop.
  
	```cpp
	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
	    // Loop body
		  // ...
	    // Implicit barrier at the end of the for loop
	}
	```
	
- `sections` Construct:
  The sections construct allows multiple code sections to be executed in parallel by different threads.
  At the end of each section within the sections construct, an implicit barrier is inserted.
  The barrier ensures that all threads complete their respective section before moving to the next sequential code outside the sections construct.

	```c#
	#pragma omp parallel sections
	{
  	  #pragma omp section
    	{
      	  // Section 1
        	// ...
		      // Implicit barrier at the end of Section 1
    	}
   	 #pragma omp section
    	{
      	  // Section 2
        	// ...
      	  // Implicit barrier at the end of Section 2
    	}
    	// ...
    	// Implicit barrier at the end of the sections construct
  }
  ```
  It's important to note that these constructs provide implicit barriers by default to synchronize the participating threads. If explicit control over synchronization is needed at different points, the programmer can add explicit barriers using the #pragma omp barrier directive at desired locations in the code.

	Implicit barriers help maintain the expected execution order and provide synchronization points without the need for explicit barrier directives at the end of each construct. However, explicit barriers may still be necessary in certain situations to enforce additional synchronization requirements or to coordinate data dependencies.
	
### Avoid the extra cost of implicit barrier

```c#
#pragma omp parallel shared(A, B, C) private(id)
{
  id = omp_get_thread_num();
  A[id] = big_calc(id);
#pragma omp barrier // explicit barrier
#pragma omp for
  for(i = 0; i < N; i++)
  {
    C[i] = big_calc3(i, A);
  } //implicit barrier here
#pragma omp for nowait // break the implicit barrier of next for
  for(i = 0; i < N; i++)
  {
    B[i] = big_calc2(C, i);
  }
  A[id] = big_calc4(id); // not use the above B[i] calculated
}
```

**Note**:The implicit barrier also incurs significant overhead in terms of execution time, so we can use `#pragma omp for nowait` to remove the implicit barrier of `#pragma omp for` when there is no need for local synchronization. In this case, `A[id] = big_calc4(id);` didn't use `B[i]` calculated from above for-loop, so we can remove the implicit barrier.

## Mutual Exclusion Synchronization

Define a block of code using `#pragma omp critical` that only one thread at a time can execute (figure below).

```c#
float res;
#pragma omp parallel
{
  float B; int i, id, num;
  id = omp_get_thread_num();
  num = omp_get_num_threads();
  for (i = id; i < niters; i+ = num)
  {
    B = big_job(i); // big job is executed in for-loop
    #pragma omp critical
    res += consume(B);
  }
}
```

<img src="./figures/IMG_3.png" alt="IMG_3" style="zoom: 35%;" />

## Atomic Synchronization

It's suitable for updating simple binary values, such as incrementing, or reading and writing a temporary value for updating.

```c
#pragma omp parallel
{
	double tmp, B;
	B = DOIT();
	tmp = big_ugly(B);
	#pragma omp atomic
	X+ = tmp;
}
```

***Notes***: The statement inside the atomic *must* be one of the following forms: ```x+ = expr```, ```x- = expr```, ```x++```, ```++x```, ```x--```, ```--x```.

## Example2_revised (Mutual Exclusion Synchronization)

```cpp
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
```
# Worksharing

- Loop Construct
- Sections/Section Construct
- Matser/Single Construct
- Task Construct

## Loop Construct

```c#
#pragma omp parallel
{
#pragma omp for
	for (i = 0; i < N; i++)
	{
		NEAT_STUFF(i);
	}
}
```

- **Static schedule**: one of the options available for loop work sharing. It determines how the iterations of a loop are divided among multiple threads in a parallel region.
  
  When using a static schedule, iterations of the loop are divided into chunks, and each thread is assigned a chunk of iterations to execute. The chunk size is typically determined based on the total number of loop iterations and the number of threads available.
  
  The key characteristic of the static schedule is that the iterations are assigned to threads in a *predetermined* and *fixed* order. The order is determined by the iteration index and the chunk size. The first thread is assigned the first chunk, the second thread gets the second chunk, and so on. This ensures that the workload is evenly distributed among the threads.
  
  Here's an example of using the static schedule in OpenMP:
  
  ```c#
  #include <omp.h>
  int main() {
    int i;
    int num_iterations = 100;
    #pragma omp parallel for schedule(static)
    for (i = 0; i < num_iterations; i++) {
        // Loop body
    }
    return 0;
  }
  
- **Dynamic schedule**:  another option for loop work sharing. It determines how loop iterations are divided among threads in a parallel region, similar to the static schedule. However, unlike the static schedule, the dynamic schedule assigns iterations to threads on-demand, rather than in a predetermined order.

  When using a dynamic schedule, the loop iterations are divided into chunks, and each thread is initially assigned a chunk of iterations. *When a thread finishes executing its assigned chunk, it requests another chunk from the remaining iterations*. This dynamic assignment of iterations continues until all iterations are completed.

  The dynamic scheduling policy offers load balancing and is useful when the iterations of the loop have varying execution times. It allows the threads to work on smaller chunks of iterations, which can help distribute the workload more evenly and prevent any single thread from being idle while others are still executing.

  Here's an example of using the dynamic schedule in OpenMP:
  
  ```c#
  #include <omp.h>
  int main() {
      int i;
      int num_iterations = 100;
      #pragma omp parallel for schedule(dynamic)
      for (i = 0; i < num_iterations; i++) {
          // Loop body
      }
      return 0;
  }
  
  ```
  
  Similar to the static schedule, you can also specify the chunk size explicitly using the `chunk` clause with the `schedule(dynamic, chunk_size)` syntax. By default, if the chunk size is not specified, the OpenMP runtime system determines a suitable value based on the number of threads and the total number of iterations. Also note that combined construct: `#pragma omp parallel` + `#parallel omp for` = `#parallel omp parallel for`.
  
  **Remark**: The difference in execution time between the static and dynamic schedules can be attributed to the inherent characteristics of the scheduling options and the nature of the computation being performed.
  
  In the static schedule, the loop iterations are divided into equal-sized chunks, and each thread is assigned a contiguous block of iterations. This allows the iterations to be distributed evenly among the threads, resulting in better load balancing. Since the workload is evenly distributed, the threads can work more efficiently and complete their assigned iterations faster, resulting in a shorter execution time.
  
  On the other hand, in the dynamic schedule, the loop iterations are distributed dynamically among the threads in smaller chunks. This can introduce additional overhead due to the dynamic nature of scheduling. The overhead includes the need for thread synchronization, task distribution, and potential load imbalance. As a result, the dynamic schedule might have a longer execution time compared to the static schedule, especially when the computation is relatively small or there is a higher overhead associated with task distribution.
  
  It's important to note that the performance difference between scheduling options can vary depending on the specific workload, the number of threads, and the characteristics of the underlying hardware. It's recommended to experiment with different scheduling options and parameters to find the most suitable configuration for your specific application and hardware setup.

## Working With Loops

- Find the *compute-intensive* loops.
  
- Make the loop iterations *independent* so they can safely execute in any order without loop-carried dependencies.
  
- Place the appropriate *OpenMP* directive and test.
  
  Here is an example:
  
  ```c#
  int i, j, A[MAX];
  j = 5;
  for (i = 0; i < MAX; i++)
  {
    j += 2;
    A[i] = big[j];
  }
  ```
  
  The example above cannot be executed with OpenMP loops for the loop body has dependencies each iteration. But we can remove that dependencies by redesign it:
  
  ```c#
   int i, j, A[MAX];
   #pragma omp parallel for
   for (i = 0; i < MAX; i++)
   {
     j = 5 + 2 * (i+1);
     A[i] = big[j];
   }
  ```
  
## Reduction operator

- reduction is a technique used to perform a parallel reduction operation on a shared variable. It allows multiple threads to independently compute partial results and then combine those results into a single value.

- Reduction operations are typically used to aggregate or combine values across iterations of a loop or other parallel constructs. Common reduction operations include summing elements, finding the maximum or minimum value, performing logical operations like bitwise OR or AND, and more.

- The OpenMP `reduction` clause provides a convenient way to specify reduction operations in parallel constructs, such as parallel loops (`parallel for`) or parallel sections (`parallel sections`). The `reduction` clause automatically handles the synchronization and reduction logic required to combine the partial results.

  Here's an example of using the `reduction` clause in OpenMP:

  ```c#
  #include <omp.h>
  
  int main() {
      int i;
      int sum = 0;
      int num_iterations = 100;
      
      #pragma omp parallel for reduction(+:sum)
      for (i = 0; i < num_iterations; i++) {
          sum += i;
      }
      printf("Sum: %d\n", sum);
      return 0;
  }
  ```

- In the example above, the `reduction(+:sum)` clause is added to the `parallel for` construct. It specifies that the reduction operation is addition (`+`) and the shared variable to be reduced is `sum`. Each thread works on a subset of loop iterations and accumulates its partial sum into the shared `sum` variable using the reduction operation (`+=`). After the parallel execution, the value of `sum` represents the sum of all loop iterations.
- OpenMP supports a variety of reduction operations, including arithmetic operations (`+`, `-`, `*`, `/`, etc.), logical operations (`&&`, `||`, `&`, `|`, `^`), and some predefined functions (`min`, `max`, etc.). You can choose the appropriate reduction operation based on the desired aggregation behavior.
- Using the `reduction` clause simplifies the process of parallelizing reduction operations, as it handles the necessary synchronization and reduction logic automatically. It allows you to write cleaner and more concise parallel code while achieving parallelism and aggregating results efficiently.
- **OpenMP: Reduction operands/initial-values**

| Operator |   Initial Value   |Operator|Initial Value|
| :---: | :---: | :---: | :---: |
|  +   |   0   | & |~0|
|   *   |   1   | \| |0|
|  -  |   0  | ^ |0|
| Min | *Largest pos num* | && |1|
| Max | *Most neg num* | \|\| |0|

The initial value of variable to be reduced is very important.

## Discussion

**In OpenMP, the scope and visibility of variables can have different effects depending on whether they are declared inside or outside of a parallel region.**

1. Variables declared outside the `#pragma omp parallel` directive:
   - These variables have a scope that extends beyond the parallel region.
   - They are shared among all threads in the parallel region, meaning that each thread can read and write to the same memory location.
   - Changes made to these variables by any thread are visible to all other threads.
   - It's important to note that proper synchronization mechanisms, such as using OpenMP constructs like `reduction` or `critical`, should be employed to avoid data races or inconsistent results when multiple threads access and modify shared variables simultaneously.
2. Variables declared inside the `#pragma omp parallel` directive:
   - These variables have a scope limited to the parallel region.
   - They are private to each thread, meaning that each thread has its own copy of the variable.
   - Each thread can independently read from and write to its private copy of the variable without affecting the copies of other threads.
   - The initial value of a private variable is undefined. If you want to initialize private variables, you should do so explicitly within the parallel region.

Here's an example illustrating the difference:

```c#
#include <omp.h>
#include <stdio.h>

int global_variable = 0;

int main() {
    #pragma omp parallel
    {
        int private_variable = omp_get_thread_num();
        global_variable += private_variable;

        printf("Thread %d: private_variable = %d, global_variable = %d\n",
               omp_get_thread_num(), private_variable, global_variable);
    }

    printf("After parallel region: global_variable = %d\n", global_variable);

    return 0;
}
```

In the example above, `global_variable` is declared outside the parallel region and is shared among all threads. Each thread has its private copy of `private_variable` declared inside the parallel region. The code increments `global_variable` by the value of `private_variable` for each thread.

After the parallel region, the value of `global_variable` will reflect the accumulated contributions from all threads. However, note that accessing `global_variable` without proper synchronization mechanisms (e.g., `atomic`, `critical`, or `reduction`) may lead to race conditions and incorrect results.

Understanding the scope and visibility of variables in OpenMP is crucial for proper parallel programming. Proper usage of shared and private variables helps ensure correctness and avoid race conditions when multiple threads are involved.

### How to avoid data races and false sharing 

- **Private Variables**: Declare variables as private to each thread to ensure that each thread operates on its own private copy of the variable.

```c#
#include <omp.h>
#include <stdio.h>

#define ARRAY_SIZE 1000

int main() {
    int shared_array[ARRAY_SIZE];
    int private_sum = 0;
    #pragma omp parallel private(private_sum)
    {
        int thread_id = omp_get_thread_num();
        // Each thread operates on its own private_sum variable
        for (int i = thread_id; i < ARRAY_SIZE; i += omp_get_num_threads()) {
            private_sum += shared_array[i];
        }
        printf("Thread %d: private_sum = %d\n", thread_id, private_sum);
    }
    return 0;
}
```

Note that the `private_sum` outside of the parallel region didn't get updated.

- **Reduction Operations**: Use reduction operations to perform parallel reductions on shared variables. OpenMP automatically handles synchronization and combination of partial results into a single value.

```c#
#include <omp.h>
#include <stdio.h>

int main() {
    int sum = 0;
    int num_iterations = 100;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < num_iterations; i++) {
        sum += i;
    }
    printf("Sum: %d\n", sum);
    return 0;
}
```

- **Atomic Operations**: Use atomic operations to ensure that specific operations are executed atomically without data races. Atomic operations are performed on shared variables and guarantee that the operation completes without interruption.

```c#
#include <omp.h>
#include <stdio.h>

int main() {
    int shared_var = 0;
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        // Perform atomic increment operation on shared_var
        #pragma omp atomic
        shared_var += thread_id;
        printf("Thread %d: shared_var = %d\n", thread_id, shared_var);
    }
    printf("After parallel region: shared_var = %d\n", shared_var);
    return 0;
}
```

- **Critical Sections**: Use critical sections to protect shared variables from simultaneous access. Only one thread can execute the critical section at a time, preventing data races.

```c#
#include <omp.h>
#include <stdio.h>
  
int main() {
    int shared_var = 0;
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        // Perform increment operation within a critical section
        #pragma omp critical
        shared_var += thread_id;
        printf("Thread %d: shared_var = %d\n", thread_id, shared_var);
    }
    printf("After parallel region: shared_var = %d\n", shared_var);
    return 0;
}
```

By employing private variables, reduction operations, atomic operations, and critical sections as necessary, you can mitigate data races and false sharing problems in OpenMP parallel regions, ensuring the correctness and performance of your parallel code.

## Example4

```c#
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
#pragma omp parallel for reduction(+:sum)
  for (i = 0; i < num_steps; i++)
  {
    long double x = (i + 0.5) * dx;
    sum += 4.0 / (1.0 + x*x);
  }
  pi = sum * dx;
  end_time = omp_get_wtime();
  printf("The numerical integration of pi is %.10Lf.\n", pi);
  printf("Execution time: %Lf seconds.", end_time - start_time);
  return 0;
} 
```

Note that using ```#pragma omp parallel for reduction(+:sum)``` is a good solution of summation, also the variable ```long double x```  we declared inside the parallel region is praivate for each thread.

## More about for-loop schedule

### omp_set_schedule

The function `omp_set_schedule` is a runtime library routine in OpenMP that allows you to change the scheduling behavior of parallel loops at runtime. It is used to set the scheduling kind and chunk size for subsequent parallel loops.

The function `omp_set_schedule` has the following syntax:

```c#
void omp_set_schedule(omp_sched_t kind, int chunk_size);
```

The parameters of `omp_set_schedule` are as follows:

- `kind`: An OpenMP enumeration type (`omp_sched_t`) that specifies the scheduling kind. It can take one of the following values:
  - `omp_sched_static`: Static scheduling.
  - `omp_sched_dynamic`: Dynamic scheduling.
  - `omp_sched_guided`: Guided scheduling.
  - `omp_sched_auto`: Compiler-selected default scheduling.
- `chunk_size`: An integer that specifies the chunk size for the specified scheduling kind. This parameter is optional and depends on the scheduling kind. It is typically used for static and guided scheduling.

Once you call `omp_set_schedule` with the desired scheduling kind and chunk size, the specified scheduling behavior will be in effect for subsequent parallel loops until it is changed again.

Here's an example demonstrating the usage of `omp_set_schedule`:

```c#
#include <omp.h>
#include <stdio.h>

int main() {
    int chunk_size = 10;
    omp_set_schedule(omp_sched_dynamic, chunk_size);

    #pragma omp parallel for
    for (int i = 0; i < 100; i++) {
        // Loop body
        printf("Thread %d executes iteration %d\n", omp_get_thread_num(), i);
    }

    return 0;
}
```

In this example, `omp_set_schedule` is used to set dynamic scheduling with a chunk size of 10. The subsequent parallel loop will follow the specified scheduling behavior.

Please note that the actual effect of `omp_set_schedule` may depend on the OpenMP implementation and runtime environment. It is recommended to consult the documentation of your specific OpenMP implementation for detailed behavior and any limitations associated with `omp_set_schedule`.



### omp_get_schedule

The function `omp_get_schedule` is a runtime library routine in OpenMP that allows you to retrieve the current scheduling information for a parallel loop. It can be used to query the scheduling kind and chunk size that are in effect for the current execution context.

The function `omp_get_schedule` has the following syntax:

```c#
void omp_get_schedule(omp_sched_t *kind, int *chunk_size);
```

The parameters of `omp_get_schedule` are as follows:

- `kind`: A pointer to an `omp_sched_t` variable where the scheduling kind will be stored. This value will be one of the following:
  - `omp_sched_static`: Static scheduling.
  - `omp_sched_dynamic`: Dynamic scheduling.
  - `omp_sched_guided`: Guided scheduling.
  - `omp_sched_auto`: Compiler-selected default scheduling.
- `chunk_size`: A pointer to an integer variable where the chunk size will be stored. This value represents the chunk size used for static and guided scheduling. For other scheduling kinds, the chunk size value will be ignored.

After calling `omp_get_schedule`, the values of `kind` and `chunk_size` will reflect the current scheduling information.

Here's an example demonstrating the usage of `omp_get_schedule`:

```c#
#include <omp.h>
#include <stdio.h>

int main() {
    omp_sched_t kind;
    int chunk_size;
    omp_get_schedule(&kind, &chunk_size);

    printf("Current Schedule: Kind = %d, Chunk Size = %d\n", kind, chunk_size);

    #pragma omp parallel for schedule(static, 10)
    for (int i = 0; i < 100; i++) {
        // Loop body
    }

    omp_get_schedule(&kind, &chunk_size);

    printf("Updated Schedule: Kind = %d, Chunk Size = %d\n", kind, chunk_size);

    return 0;
}
```

In this example, `omp_get_schedule` is used to retrieve the current scheduling information before and after a parallel loop. The obtained values are then printed to show the current schedule kind and chunk size.

Please note that `omp_get_schedule` will return the scheduling information for the current execution context, which may be affected by the OpenMP directives or API calls preceding it. It provides a way to programmatically access the scheduling information for diagnostic or reporting purposes.

## Section Construct

### Introduction

The `sections` construct in OpenMP is a directive that allows you to divide a parallel region into multiple named sections, which can be executed by different threads concurrently. It provides a way to assign different tasks or code sections to specific threads within the parallel region.

The syntax of the `sections` construct is as follows:

```c#
#pragma omp parallel sections [clause]
{
    #pragma omp section
    {
        // Code for section 1
    }
    #pragma omp section
    {
        // Code for section 2
    }
    // ...
}
```

Key points about the `sections` construct:

- The `sections` construct divides the code block following it into individual sections, each marked with the `section` directive.
- Each `section` block represents a specific task or code section that can be executed concurrently by different threads.
- The number of sections can be determined statically or dynamically at runtime.
- The threads within the parallel region execute the sections in parallel, and the specific assignment of sections to threads may vary across different runs or implementations.
- The `sections` construct introduces an implicit barrier at the end, ensuring that all threads synchronize before continuing beyond the construct.

Example usage of the `sections` construct:

```c#
#pragma omp parallel sections
{
    #pragma omp section
    {
        // Code for section 1
    }
    #pragma omp section
    {
        // Code for section 2
    }
    // Code executed concurrently by other threads
    #pragma omp section
    {
        // Code for section 3
    }
    // ...
}
```

In this example, the `sections` construct divides the code block following it into three sections. Each section represents a specific task or code segment that can be executed concurrently by different threads within the parallel region. The exact assignment of sections to threads is determined at runtime.

The `sections` construct is useful for parallelizing independent tasks or code sections, enabling different threads to work on different parts of the problem simultaneously. It provides a structured and controlled way to assign specific sections of code to individual threads, facilitating parallel execution within the parallel region.

### Advantages

1. **Parallelization of independent tasks**: The `sections` construct enables you to divide the work within a parallel region into independent sections. Each section can be executed concurrently by different threads, maximizing parallelism and speeding up the overall execution.
2. **Explicit assignment of tasks**: By using the `section` directive within the `sections` construct, you can explicitly assign specific tasks or code sections to different threads. This control allows you to distribute the workload efficiently among the available threads based on the characteristics of the tasks.
3. **Structured approach**: The `sections` construct provides a structured and organized way to parallelize and manage different tasks or code sections. It promotes code readability and maintainability by clearly indicating the parallel sections and their respective assignments.

### Disadvantages

1. **Limited for irregular workloads**: The `sections` construct is most suitable for cases where the workload can be divided into known and evenly distributed sections. It may not be well-suited for irregular workloads or cases where dynamic load balancing is required.
2. **Limited scalability**: The scalability of the `sections` construct depends on the number of available threads and the granularity of the sections. If the number of sections is significantly smaller than the number of threads, it may result in underutilization of resources and lower parallel efficiency.

### Tips

1. **Identify independent sections**: Determine which parts of the code can be executed independently and parallelized. Ensure that there are no data dependencies or shared resources that could cause conflicts among sections.
2. **Balance the workload**: Try to divide the workload into sections of roughly equal size to achieve load balance among the participating threads. This helps ensure that all threads make similar progress and maximize parallel efficiency.
3. **Avoid excessive synchronization**: Minimize the need for explicit synchronization within sections, as it can introduce additional overhead and hinder parallel performance. Synchronization should only be used when necessary to maintain correctness or enforce dependencies.
4. **Consider the granularity**: Choose an appropriate granularity for the sections based on the nature of the tasks and the available hardware resources. Fine-grained sections may provide more parallelism but could result in increased synchronization overhead.
5. **Profile and tune**: Measure the performance of your parallel code using the `sections` construct and profile it to identify any bottlenecks or areas for optimization. Experiment with different section assignments and parallelization strategies to achieve the best performance for your specific application.

Remember that the effectiveness and performance of the `sections` construct depend on the characteristics of your code, the workload distribution, and the available hardware resources. It is recommended to experiment, profile, and benchmark your code to find the optimal use of the `sections` construct in your specific scenario.

## Master Construct

The `master` construct in OpenMP is a directive used to specify a block of code that should be executed by only the master thread in a parallel region. The master thread is typically the thread with thread ID 0, although it can be explicitly set using the `omp_set_num_threads` function.

The syntax of the `master` construct is as follows:

```c#
#pragma omp master
{
    // Code executed by the master thread
}
```

Key points about the `master` construct:

- Only the master thread executes the code block within the `master` construct. Other threads in the team skip this code block.
- The `master` construct provides a way to designate a specific task or section of code that should be executed by a single thread, typically for tasks that should be performed by only one thread, such as initialization or I/O operations.
- The master thread completes the execution of the `master` construct before any other threads proceed to the next parallel region or construct.
- There exists implicit barrier at the end of `master` construct.
- The `master` construct does not introduce an implicit barrier. If synchronization is required after the `master` construct, an explicit barrier should be used.

Example usage of the `master` construct:

```c#
#pragma omp parallel
{
    // Code executed by all threads
    #pragma omp master
    {
        // Code executed only by the master thread
        // Typically used for initialization or serial operations
    }
    // Code executed by all threads
    #pragma omp barrier
    // Explicit barrier to synchronize all threads
}
```

In the example, the code within the `master` construct is executed only by the master thread, while other threads skip that code block. The explicit barrier after the `master` construct ensures synchronization among all threads before continuing to the next parallel region or construct.

The `master` construct is useful for situations where certain tasks should be performed by a single thread, while allowing other threads to continue with their work.

## Single Construct

The `single` construct in OpenMP is a directive used to specify a block of code that should be executed by only one thread in a parallel region. Unlike the `master` construct, which designates the master thread as the exclusive executor, the `single` construct allows any thread to execute the enclosed code block.

The syntax of the `single` construct is as follows:

```c#
#pragma omp single [clause]
{
    // Code executed by only one thread
}
```

Key points about the `single` construct:

- The `single` construct designates a block of code that should be executed by only one thread. It does not specify which specific thread will execute the code block.
- When multiple `single` constructs are encountered, each construct allows only one thread to execute its corresponding code block. However, different threads may execute different `single` constructs concurrently.
- If multiple threads encounter a `single` construct simultaneously, only one thread will enter the construct while other threads skip it. Which thread executes the construct is implementation-dependent.
- There exists implicit barrier at the end of ` single` construct.
- The `single` construct can have an optional clause that provides additional control over how the construct behaves. Common clauses used with the `single` construct include `nowait` and `copyprivate`.

Example usage of the `single` construct:

```c#
#pragma omp parallel
{
    // Code executed by all threads
    #pragma omp single
    {
        // Code executed by only one thread
        // The specific thread executing this block is implementation-dependent
    }
    // Code executed by all threads
}
```

In the example, the code within the `single` construct is executed by only one thread, while other threads skip that code block. The specific thread executing the block is not predetermined and may vary across different runs or implementations.

The `single` construct is useful when there is a section of code that should be executed by only one thread, without specifying which thread specifically. It is commonly used for tasks such as file I/O, allocating shared resources, or performing operations that require exclusive access.

## Lock Routines

### Introduction
Lock routines in OpenMP are a set of functions provided by the OpenMP API for thread synchronization and mutual exclusion. They allow you to protect critical sections of code from simultaneous access by multiple threads, ensuring that only one thread can access the protected region at a time. Locks are useful when you need to prevent race conditions and maintain data integrity in multithreaded programs.

OpenMP provides several lock routines that you can use to implement different locking mechanisms:

1. `omp_init_lock()`:
   - This function initializes a simple lock, which is a basic binary lock with two states: locked and unlocked.
   - It is typically called before the parallel region to initialize the lock.
2. `omp_destroy_lock()`:
   - This function destroys a lock, releasing any resources associated with it.
   - It should be called after the parallel region or when the lock is no longer needed.
3. `omp_set_lock()`:
   - This function acquires a lock, blocking other threads from entering the protected region until the lock is released.
   - If the lock is already locked by another thread, the calling thread will wait until the lock becomes available.
4. `omp_unset_lock()`:
   - This function releases a lock, allowing other threads to acquire it and enter the protected region.
   - It should always be called after the critical section of code has been executed.
5. `omp_test_lock()`:
   - This function attempts to acquire a lock, but it returns immediately instead of waiting if the lock is already locked by another thread.
   - It returns a boolean value indicating whether the lock was successfully acquired.

Lock routines can be used to protect critical sections of code where data integrity needs to be ensured. 

### Example

Here's an example demonstrating the use of lock routines in OpenMP to protect a critical section of code that updates a shared counter variable:

```c#
#include <omp.h>
#include <stdio.h>

int main() {
    omp_lock_t lock;
    int shared_counter = 0;
    // Initialize the lock
    omp_init_lock(&lock);
    #pragma omp parallel num_threads(4)
    {
        int private_counter = 0;
        // Each thread updates its private counter
        for (int i = 0; i < 100; i++) {
            private_counter++;
        }
        // Acquire the lock to update the shared counter
        omp_set_lock(&lock);
        shared_counter += private_counter;
        omp_unset_lock(&lock);
    }
    // Destroy the lock
    omp_destroy_lock(&lock);
    printf("Shared counter: %d\n", shared_counter);
    return 0;
}
```

In this example, `omp_lock_t` is a data type used to declare a lock variable. The program creates a lock using `omp_init_lock(&lock)` and initializes a shared counter variable (`shared_counter`). Each thread has its own private counter (`private_counter`) to accumulate individual counts. The private counters are then added to the shared counter within a critical section protected by lock routines (`omp_set_lock(&lock)` and `omp_unset_lock(&lock)`).

By acquiring the lock before updating the shared counter, the critical section ensures that only one thread at a time can modify the shared variable, preventing data races and ensuring the integrity of the counter. After all threads have finished updating their private counters and released the lock, the shared counter is printed.

Please note that in this example, the shared counter is updated sequentially within the critical section due to the lock. If the updates to the shared counter do not have strict ordering requirements, you might consider other synchronization mechanisms like atomic operations or OpenMP reduction clauses to improve parallelism and avoid potential contention.

### Lock Array

In OpenMP, a lock array refers to an array of lock variables (`omp_lock_t`) used to protect multiple shared resources or critical sections in a parallel program. Each lock in the array corresponds to a specific resource or critical section and is used to provide mutual exclusion among threads accessing that resource.

Using a lock array can be beneficial when you have multiple shared resources or critical sections that need to be protected independently. By assigning a lock to each resource, you can ensure that only one thread at a time can access a particular resource, preventing data races and maintaining data integrity.

Here's an example of using a lock array to protect multiple critical sections:

```c#
#include <omp.h>
#include <stdio.h>

#define NUM_SECTIONS 5

int main() {
    omp_lock_t lock_array[NUM_SECTIONS];
    int shared_data[NUM_SECTIONS] = {0};
    // Initialize the lock array
    for (int i = 0; i < NUM_SECTIONS; i++) {
        omp_init_lock(&lock_array[i]);
    }

    #pragma omp parallel num_threads(4)
    {
        int thread_id = omp_get_thread_num();
        // Each thread updates its assigned section
        for (int i = 0; i < NUM_SECTIONS; i++) {
            if (i % 4 == thread_id) {
                // Acquire the lock for the current section
                omp_set_lock(&lock_array[i]);
                
                // Critical section: Update shared data for the current section
                shared_data[i] += thread_id + 1;

                // Release the lock for the current section
                omp_unset_lock(&lock_array[i]);
            }
        }
    }
    // Destroy the lock array
    for (int i = 0; i < NUM_SECTIONS; i++) {
        omp_destroy_lock(&lock_array[i]);
    }

    // Print the updated shared data
    for (int i = 0; i < NUM_SECTIONS; i++) {
        printf("Section %d: %d\n", i, shared_data[i]);
    }
    return 0;
}
```

In this example, we have an array of locks (`lock_array`) corresponding to an array of shared data (`shared_data`) divided into sections. Each thread updates its assigned section by acquiring the corresponding lock, performing the necessary operations within the critical section, and releasing the lock. By assigning a lock to each section, we ensure that only one thread can access a section at a time, preventing conflicts and ensuring the integrity of the shared data.

Remember to initialize the lock array using `omp_init_lock()` before using it, destroy the locks using `omp_destroy_lock()` when they are no longer needed, and appropriately acquire and release the locks using `omp_set_lock()` and `omp_unset_lock()` respectively within critical sections of code.

Using a lock array provides a structured and organized approach to protect multiple critical sections or shared resources, allowing for concurrent access while maintaining data consistency and avoiding data races.

## Runtime Library Routines

### `omp_in_parallel()`

The function `omp_in_parallel()` is an OpenMP runtime library routine used to determine if the caller thread is currently executing in a parallel region. It returns a boolean value indicating whether the thread is inside a parallel region or not.

The syntax for `omp_in_parallel()` is as follows:

```c#
int omp_in_parallel();
```

The return value of `omp_in_parallel()` is `1` if the caller thread is executing in a parallel region, and `0` if it is executing in a serial region (i.e., not in a parallel region).

Here's an example illustrating the usage of `omp_in_parallel()`:

```c#
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
```

In this example, the code is executed within a parallel region defined by the `#pragma omp parallel` directive. Each thread checks its execution status using `omp_in_parallel()` and prints a corresponding message indicating whether it is executing in a parallel or serial region.

The output of the program would be something like:

```bash
Thread is executing in a parallel region.
Thread is executing in a parallel region.
Thread is executing in a parallel region.
...
```

The `omp_in_parallel()` function is useful when you need to conditionally execute certain code based on whether the thread is inside a parallel region or a serial region. It allows you to control the behavior of the code based on the parallelism context.

**Remark**: If you use `omp_get_num_threads()` outside the parallel region, the answer will be one.

### `omp_get_dynamic()`

**Dynamic mode**: Dynamic mode allows the OpenMP implementation to adapt the number of threads dynamically to optimize resource utilization and load balancing. It can be particularly useful in scenarios where the workload varies over time or when the number of available processors or threads changes dynamically.

The function `omp_get_dynamic()` is an OpenMP runtime library routine used to query the dynamic adjustment of the number of threads in a parallel region. It allows you to check whether the dynamic adjustment of the number of threads is enabled or disabled.

The syntax for `omp_get_dynamic()` is as follows:

```c#
int omp_get_dynamic();
```

The return value of `omp_get_dynamic()` is an integer that indicates the current status of dynamic thread adjustment. Here's what the return value represents:

- If the return value is `0`, it means that dynamic adjustment is disabled. The number of threads used in a parallel region will be determined by the value set using `omp_set_num_threads()` or the default number of threads set by the OpenMP implementation.
- If the return value is `1`, it means that dynamic adjustment is enabled. The number of threads used in a parallel region may be adjusted dynamically during runtime based on the environment or the `OMP_NUM_THREADS` environment variable.

Here's an example demonstrating the usage of `omp_get_dynamic()`:

```c#
#include <omp.h>
#include <stdio.h>

int main() {
    int dynamic = omp_get_dynamic();

    if (dynamic) {
        printf("Dynamic adjustment of threads is enabled.\n");
    } else {
        printf("Dynamic adjustment of threads is disabled.\n");
    }

    return 0;
}
```

In this example, the code queries the status of dynamic thread adjustment using `omp_get_dynamic()`. It prints a message indicating whether dynamic adjustment is enabled or disabled based on the return value of the function.

The output of the program would be something like:

```bash
Dynamic adjustment of threads is enabled.
```

The status of dynamic thread adjustment can be useful to know when designing and debugging your parallel programs. It allows you to understand whether the number of threads in a parallel region can be adjusted dynamically during runtime. Depending on your requirements, you can use this information to make decisions about thread management and workload distribution.

#### `omp_set_dynamic()`

The function `omp_set_dynamic()` is an OpenMP runtime library routine used to enable or disable dynamic adjustment of the number of threads in a parallel region. It allows you to control whether the number of threads can be adjusted dynamically during runtime.

The syntax for `omp_set_dynamic()` is as follows:

```c#
void omp_set_dynamic(int dynamic_threads);
```

The `dynamic_threads` argument is an integer value that determines whether dynamic thread adjustment is enabled or disabled. Here's how the argument can be set:

- If `dynamic_threads` is non-zero, dynamic adjustment is enabled. The number of threads used in a parallel region may be adjusted dynamically during runtime.
- If `dynamic_threads` is zero, dynamic adjustment is disabled. The number of threads used in a parallel region will be determined by the value set using `omp_set_num_threads()` or the default number of threads set by the OpenMP implementation.
##### Test and Output
```c#
#include <omp.h>
#include <stdio.h>

int main() {
    omp_set_dynamic(1);  // Enable dynamic mode
    #pragma omp parallel
    {
        // Code inside the parallel region
        if (omp_get_dynamic()) {
            printf("Thread %d: Dynamic mode is enabled.\n", omp_get_thread_num());
        } else {
            printf("Thread %d: Dynamic mode is disabled.\n", omp_get_thread_num());
        }
    }

    return 0;
}
```

```bash
Thread 0: Dynamic mode is enabled.
Thread 1: Dynamic mode is enabled.
Thread 2: Dynamic mode is enabled.
...
```

### `omp_get_num_procs()`

The function `omp_num_procs()` is an OpenMP runtime library routine used to retrieve the number of available processors or processing units in the system. It provides information about the total number of processors that can potentially be used for parallel execution.

The syntax for `omp_num_procs()` is as follows:

```c#
int omp_get_num_procs();
```

The return value of `omp_num_procs()` is an integer that represents the number of available processors or processing units in the system.

Here's an example demonstrating the usage of `omp_num_procs()`:

```c#
#include <omp.h>
#include <stdio.h>

int main() {
    int num_procs = omp_get_num_procs();
    printf("Number of available processors: %d\n", num_procs);
    return 0;
}
```

In this example, the code calls `omp_num_procs()` to retrieve the number of available processors. It then prints the value to the console.

The output of the program would be something like:

```bash
Number of available processors: 8
```

The value returned by `omp_num_procs()` can be useful in various scenarios, such as determining the maximum number of threads to use in a parallel region or assessing the available resources for load balancing. It provides an insight into the system's processing capacity and can help in optimizing the performance of parallel programs.

Please note that `omp_num_procs()` reports the number of available processors as seen by the OpenMP implementation and may not necessarily reflect the physical number of cores or processors on the machine. The actual behavior can vary depending on the system and OpenMP implementation.

## Example5

```c#
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
```

# Environment Variables

Environment variables in OpenMP are global variables that can be set in the execution environment to control the behavior of OpenMP programs. These variables provide a way to configure and customize the execution of OpenMP programs without modifying the source code. All environment variables are documented at `~/.bashrc` or `~/.zshrc`.

-  `export OMP_NUM_THREADS=<desired_number_of_threads>`: set the number of threads.

- `export OMP_DYNAMIC=<TRUE/FALSE>`: enable or disable the dynamic mode.

- `export OMP_NESTED=<TRUE/FALSE>`: enable of disable nested parallelism.

- `export OMP_SCHEDULE="<type>,<chunk_size>"`: specify the default scheduling policy for *parallel loop constructs* when the `schedule` clause is not explicitly specified.

- `export OMP_PROC_BIND=<option>`: control the binding of OpenMP threads to physical processing resources such as processor cores. It specifies the desired thread-to-core binding policy for parallel execution.

- `export OMP_STACKSIZE=<storage>`: control the size of the stack allocated for each thread in an OpenMP parallel region. It allows you to specify the stack size to accommodate the memory requirements of your parallel program.

- `export OMP_WAIT_POLICY=<PASSIVE/POSITIVE>`: control the behavior of threads waiting for synchronization points, such as barriers or critical sections. It determines the waiting policy used by the OpenMP runtime when a thread is waiting for other threads to complete their work.

  ```bash
  # Environment Variables for OpenMP
  export OMP_NUM_THREADS=24
  export OMP_DYNAMIC=FALSE
  export OMP_NESTED=TRUE
  export OMP_SCHEDULE="static,4"
  export OMP_PROC_BIND=master
  export OMP_STACKSIZE=256M
  export OMP_WAIT_POLICY=POSITIVE
  ```

  ## Data Environment

  
