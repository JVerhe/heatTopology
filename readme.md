## The repository has the following layout:

- etc: contains files used in development for various purposes
- proto: contains the python prototype of the project
- src: contains the files for the C++ project
- test: contains automatic unit tests

## To compile the project:

Requirements:

`sudo apt-get install libeigen3-dev`

`sudo apt-get install cmake`

1. cd build
2. cmake ..
3. cmake --build . (everytime you want to compile the project)

Compilation flags can be changed in CMakeLists.txt
