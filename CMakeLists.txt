cmake_minimum_required(VERSION 3.10)
project(HeatTopology_team2)

set(CMAKE_CXX_STANDARD 20)

add_executable(HELLO src/hello.cpp)
target_compile_options(HELLO PRIVATE -Wall -Wextra)