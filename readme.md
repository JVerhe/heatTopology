## The repository has the following layout:

- etc: contains files used in development for various purposes
- proto: contains the python prototype of the project
- src: contains the files for the C++ project
- test: contains automatic unit tests

## To compile the project:

Requirements:

The Eigen c++ library for linear algebra from: [text](https://eigen.tuxfamily.org)
`sudo apt-get install libeigen3-dev`

The CMake software build system from: [text](https://cmake.org/)
`sudo apt-get install cmake`

For compiling the first time:

1. `cd build`
2. `cmake ..`
3. `cmake --build .` (everytime you want to compile the project)

Compilation flags can be changed in CMakeLists.txt

- Deleting all the log files can be done by running the command:

`cmake --build . --target clean_logs`
