#ifndef adjointMethod_hpp
#define adjointMethod_hpp

#include <cassert>
#include <Eigen/Dense>
#include <map>
#include <cmath>


using namespace Eigen;

VectorXd adjoint(VectorXd T, VectorXd v, std::map<int, std::vector<int>> corners, double p = 1) {
    int elements = v.size();
    VectorXd gradJv(elements);
    double km = 65;
    double kp = 0.2;
    Matrix4d K0{
        {2. / 3, -1. / 6, -1. / 3, -1. / 6},
        {-1. / 6, 2. / 3, -1. / 6, -1. / 3},
        {-1. / 3, -1. / 6, 2. / 3, -1. / 6},
        {-1. / 6, -1. / 3, -1. / 6, 2. / 3}
    };

    for (int el = 0; el < elements; el++) {
        std::vector<int> elementNodes = corners[el];
        Vector4d Te;
        for (int i = 0; i < 4; i++) {
            Te(i) = T(elementNodes[i]);
        }
        gradJv(el) = (Te.transpose() * K0 * Te);
        gradJv(el) *= -0.5 * p * (std::pow(v[el], p - 1)) * (km - kp);
    }

    return gradJv;
}

#endif