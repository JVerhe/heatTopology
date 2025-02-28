#include <iostream>
#include <stdlib.h>

#define BOOST_TEST_MODULE test1
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(dummy_test)
{
    BOOST_TEST(1 == 1);
}
// #include <iostream>
// int main() {
//     std::cout << "This works!" << std::endl;
//     return 0;
// }