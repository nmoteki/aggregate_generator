cmake_minimum_required (VERSION 3.10.0)
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -DNDEBUG -DEIGEN_NO_DEBUG")
include_directories(/mnt/c/Users/moteki/WSL_Share_Ubuntu1804/C++/eigen_3_3_5)

project(aggregate_gen)

add_executable(aggregate_gen_main aggregate_gen_main.cpp)

add_library(aggregate_gen_CCA aggregate_gen_CCA.cpp)
add_library(aggregate_gen_SA aggregate_gen_SA.cpp)

target_link_libraries(aggregate_gen_main aggregate_gen_CCA)
target_link_libraries(aggregate_gen_CCA  aggregate_gen_SA)
