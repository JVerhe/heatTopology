## Requirements

The Eigen c++ library for linear algebra from [link](https://eigen.tuxfamily.org)

`sudo apt-get install libeigen3-dev`

The CMake software build system from [link](https://cmake.org/)

`sudo apt-get install cmake`

For plotting the results of the script, some Python packages are required. These can be installed from the main repository via the following command:

`pip install -r requirements.txt`

## Compiling the project

For compiling the first time:

1. `cd build`
2. `cmake ..`
3. `cmake --build .` (everytime you want to compile the project)

Compilation flags can be changed in CMakeLists.txt

Compilation flags can be changed in CMakeLists.txt

## Running the main.cpp

- A config file inside the folder `/build/config` has to be present. It should look like this:

```
points: Number of Discretization points in one dimension (int)
penalty: Metal fraction penalty exponent (double)
filtering: The applied filtering technique (0 = no filtering, 1 = sensitivity filtering, 2 = density filtering)
output: (0 = no output, 1 = visualise output, 2 = run simple benchmark)
```

For example a `config.txt` could look like this: (by default this specific configuration is included in the project)

```
points 30
penalty 2
filtering 2
output 1
```

The main file can then be executed by running `./main config`

## Deleting log files

Deleting all the log files in `/build/output` can be done inside the build directory by running the command:

`cmake --build . --target clean_logs`

Deleting log files can also be done via the Python visualization script.

## Running tests

- All the tests should be placed in the `/test` directory

- Compiling and running tests:

1. `cd build`
2. `cmake -DBUILD_TESTS=ON ..` (can be set to `OFF` when you don't want to compile the tests)
3. `cmake --build . --target run_tests` (everytime you want to run only the tests)
