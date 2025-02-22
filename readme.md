## The repository has the following layout:

- etc: contains files used in development for various purposes
- proto: contains the python prototype of the project
- src: contains the files for the C++ project
- test: contains automatic unit tests

To compile the project:

1. cd build
2. sudo apt-get install libeigen3-dev
3. cmake ..
4. cmake --build . (everytime you want to compile the project)

Compilation flags can be changed in CMakeLists.txt
TODO library management
