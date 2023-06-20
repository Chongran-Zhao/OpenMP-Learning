cmake_minimum_required(VERSION 3.9)

project(openmp-demo)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 11)

# Look for OpenMP package
find_package(OpenMP REQUIRED)

# Set compiler flags
set(CMAKE_CXX_FLAGS "-g -Wall -O3 ${OpenMP_CXX_FLAGS} -pg")

include_directories(.)

add_executable(example1 example1.cpp)
add_executable(example2 example2.cpp)

target_link_libraries(example1 OpenMP::OpenMP_CXX)
target_link_libraries(example2 OpenMP::OpenMP_CXX)

set_target_properties(example1 PROPERTIES LINK_FLAGS "-pg")
set_target_properties(example2 PROPERTIES LINK_FLAGS "-pg")
