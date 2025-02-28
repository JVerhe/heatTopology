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

set(LOG_DIR "${CMAKE_BINARY_DIR}/output")
add_custom_target(clean_logs
    COMMAND ${CMAKE_COMMAND} -E remove -f ${LOG_DIR}/*.txt
    COMMENT "Deleting all log files in ${LOG_DIR}"
)

# This flag can be set on or off, by default it is off
option(BUILD_TESTS "Build test executables" OFF)

if(BUILD_TESTS)
    set(TEST_DIR "${CMAKE_SOURCE_DIR}/test")
    file(GLOB TEST_SOURCES "${TEST_DIR}/*.cpp")

    # List to store test executables
    set(TEST_EXECUTABLES "")

    foreach(TEST_FILE ${TEST_SOURCES})
        get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
        add_executable(${TEST_NAME} ${TEST_FILE})
        list(APPEND TEST_EXECUTABLES ${TEST_NAME})  # Store test names
    endforeach()


    add_custom_target(run_tests
        ${TEST_EXECUTABLES}
        COMMENT "Executing all compiled test cases."
        DEPENDS ${TEST_EXECUTABLES}  # Ensure tests are built before running
    )
endif()
