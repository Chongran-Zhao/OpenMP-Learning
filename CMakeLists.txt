cmake_minimum_required(VERSION 3.20)

project(openmp-demo)

set(CMAKE_CXX_STANDARD 11)

# Look for OpenMP package
set(OpenMP_ROOT "/opt/homebrew/Cellar/libomp/16.0.6")
find_package(OpenMP REQUIRED)

# Set compiler flags
set(CMAKE_CXX_FLAGS "-g -Wall -O3 ${OpenMP_CXX_FLAGS}")

include_directories(.)

add_executable(example1 example1.cpp)
add_executable(example2 example2.cpp)
add_executable(example3 example3.cpp)
add_executable(example2_revised example2_revised.cpp)
add_executable(example4 example4.cpp)

target_link_libraries(example1 OpenMP::OpenMP_CXX)
target_link_libraries(example2 OpenMP::OpenMP_CXX)
target_link_libraries(example3 OpenMP::OpenMP_CXX)
target_link_libraries(example2_revised OpenMP::OpenMP_CXX)
target_link_libraries(example4 OpenMP::OpenMP_CXX)
