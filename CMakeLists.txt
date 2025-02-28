cmake_minimum_required(VERSION 3.10)
project(HeatTopology_team2)

set(CMAKE_CXX_STANDARD 20)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(/usr/include/eigen3)

include_directories(src/)
# The following lines define an executable and compilation flags
add_executable(main src/main.cpp)
target_compile_options(main PRIVATE -O3 -march=native -pedantic)

# Set the log directory inside the build folder
set(LOG_DIR "${CMAKE_BINARY_DIR}/output")

# Custom target to delete all .txt log files
add_custom_target(clean_logs
    COMMAND ${CMAKE_COMMAND} -E remove -f ${LOG_DIR}/*.txt
    COMMENT "Deleting all log files in ${LOG_DIR}"
)