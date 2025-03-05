## The repository has the following layout:

- etc: contains files used in development for various purposes
- proto: contains the python prototype of the project
- src: contains the files for the C++ project
- test: contains automatic unit tests

## To compile the project:

Requirements:

The Eigen c++ library for linear algebra from [link](https://eigen.tuxfamily.org)
`sudo apt-get install libeigen3-dev`

The CMake software build system from [link](https://cmake.org/)
`sudo apt-get install cmake`

For compiling the first time:

1. `cd build`
2. `cmake ..`
3. `cmake --build .` (everytime you want to compile the project)

Compilation flags can be changed in CMakeLists.txt

## Running the main.cpp

- A config file inside the folder `/build/config/` has to be present. It should look like this.

```
Number of Discretization points in one dimension (int)
Metal fraction penalty exponent (double)
Filtering technique (0 = no filtering, 1 = sensitivity filtering, 2 = density filtering)
```

For example a `config1.txt` could look like this:

```
20
2
2
```

The main file can be executed by running `./main config1`

## Running tests

- All the tests should be placed in the `/test` directory

- Compiling and running tests:

1. `cd build`
2. `cmake -DBUILD_TESTS=ON ..` (can be set to `OFF` when you don't want to compile the tests)
3. `cmake --build .`
4. `cmake --build . --target run_tests` (everytime you want to compile the tests)

## Deleting log files

- Deleting all the log files in `/build/output` can be done inside the build directory by running the command:

`cmake --build . --target clean_logs`
