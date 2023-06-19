cmake_minimum_required(VERSION 3.9)

project(openmp-demo)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 11)

# Look for OpenMP package
find_package(OpenMP REQUIRED)

# Set compiler flags
set(CMAKE_CXX_FLAGS "-g -Wall -O3 ${OpenMP_CXX_FLAGS} -pg")

include_directories(.)

add_executable(hello_world hello_world.cpp)

target_link_libraries(hello_world OpenMP::OpenMP_CXX)

set_target_properties(hello_world PROPERTIES LINK_FLAGS "-pg")
