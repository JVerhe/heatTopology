#ifndef testFramework_hpp
#define testFramework_hpp

#include <exception>
#include <string>
#include <stdexcept>
#include <type_traits>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <vector>

namespace tf {


    template <typename T>
    void check(const T& cond, int line, std::string file)
    {
        if (!cond)
        {
            std::cout << file << std::endl;
            std::ostringstream oss;
            oss << "Error at line " << line << ": Assertion is not correct";
            throw std::invalid_argument(oss.str());
        }
    }

    template <typename T>
    void compare(const T& value1, const T& value2, int line, std::string file)
    {
        if (value1 != value2)
        {
            std::cout << file << std::endl;
            std::ostringstream oss;
            oss << "Error at line " << line << ": Values are not equal: " << value1 << " != " << value2;
            throw std::invalid_argument(oss.str());
        }
    }

    template <typename T>
    void compareTolerance(const T& value1, const T& value2, double tol, int line, std::string file)
    {
        if (abs(value1 - value2) > tol)
        {
            std::cout << file << std::endl;
            std::ostringstream oss;
            oss << "Error at line " << line << ": Values are not equal: " << value1 << " != " << value2 << " with tolerance = " << tol;
            throw std::invalid_argument(oss.str());
        }
    }
    
    template <typename T>
    void print_2D_vector(const std::vector<std::vector<T>>& vec) {
        for (const auto& row : vec) {
            for (int val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl; // Newline after each row
        }
    }


}


#endif