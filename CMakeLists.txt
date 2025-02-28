cmake_minimum_required(VERSION 3.10)
project(HeatTopology_team2)

set(CMAKE_CXX_STANDARD 20)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(src/)
# The following lines define an executable and compilation flags
add_executable(main src/main.cpp)
target_compile_options(main PRIVATE -Wextra -Wall)
