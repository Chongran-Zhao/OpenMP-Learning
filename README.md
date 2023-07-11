# OpenMP Learning Notes
- I use this repo to update my progress in learning OpenMP.
- The code demo will follow content in lecture slides from UT at Austin.
- The video on Youtube is "Introduction to OpenMP"
- website: https://www.youtube.com/@OpenMPARB/videos

## Basic Concepts

- Understand the concept of shared memory programming and parallelization.
- Learn about the fork-join model of OpenMP and the concept of threads.

## Directive-based Parallelism

- Explore the OpenMP directives, such as `parallel`, `for`, `sections`, `single`, and `task`, to specify parallel regions and work-sharing constructs.
- Understand how to use these directives to parallelize loops, divide work into sections, and perform tasks efficiently.

## Runtime Routines
   - Familiarize yourself with common OpenMP runtime routines for controlling parallel execution, thread management, and synchronization.
   - Study routines like `omp_get_thread_num()`, `omp_get_num_threads()`, `omp_get_max_threads()`, `omp_set_lock()`, `omp_unset_lock()`, and `omp_barrier()`.
   - Understand how these routines can be used to obtain thread and team information, manage locks for critical sections, and synchronize threads.
## Thread-Safety and Synchronization
   - Learn about race conditions and the need for synchronization in multithreaded programs.
   - Explore the concepts of mutual exclusion and critical sections.
   - Understand how to use lock routines (`omp_init_lock()`, `omp_destroy_lock()`, etc.) to protect shared resources and prevent data races.

## More Advanced Topics

1. **Nested Parallelism:** Learn about the concept of nested parallelism in OpenMP, which allows parallel regions to be nested within other parallel regions. Understand the implications, benefits, and considerations when using nested parallelism, and explore functions like `omp_get_nested()` and `omp_set_nested()`.
2. **Dynamic Thread Adjustment:** Explore the dynamic adjustment of the number of threads during runtime using functions like `omp_get_dynamic()` and `omp_set_dynamic()`. Understand how these functions can be used to enable or disable dynamic adjustment of the number of threads in a parallel region.
3. **Tasking:** Gain knowledge of task-based parallelism in OpenMP, which allows you to define and execute tasks independently. Learn about functions like `omp_task()` and `omp_taskwait()` to create and synchronize tasks, and understand the concepts of task dependencies and task scheduling.
4. **Atomic Operations:** Familiarize yourself with atomic operations provided by OpenMP, such as `omp_atomic`, which ensure atomicity of memory operations. Understand how to use atomic operations to protect critical sections without the need for locks.
5. **Memory Model and Data Sharing:** Dive deeper into the OpenMP memory model, memory consistency, and data sharing concepts. Learn about the `shared`, `private`, `firstprivate`, `lastprivate`, and `default` clauses used in OpenMP directives to control data sharing and memory consistency.
6. **OpenMP Environment Variables:** Explore various environment variables that can be used to control the behavior of an OpenMP program, such as the number of threads, nested parallelism, dynamic adjustment, and more. Understand how to set and use environment variables to customize the execution of your OpenMP program.
