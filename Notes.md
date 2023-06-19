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
