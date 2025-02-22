cmake_minimum_required(VERSION 3.10)
project(HeatTopology_team2)

set(CMAKE_CXX_STANDARD 20)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# The following lines define an executable and compilation flags
add_executable(adjointMethod src/adjointMethod.cpp)
target_compile_options(adjointMethod PRIVATE -Wextra -Wall)

add_executable(test src/test.cpp)
target_compile_options(test PRIVATE -Wextra -Wall)