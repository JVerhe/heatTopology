#ifndef optHelper_hpp
#define optHelper_hpp

#include <cassert>
#include <Eigen/Dense>
#include <map>
#include <cmath>



Eigen::VectorXd fill_in_k(const Eigen::VectorXd& v, double k_max, double k_min, double p) {
    /**
     * Computes the material property values (k) based on the density vector (v).
     *
     * @param v: Vector of densities (values between 0 and 1).
     * @param k_max: Maximum stiffness value.
     * @param k_min: Minimum stiffness value.
     * @param p: Penalization factor.
     * @return Eigen::VectorXd containing computed k values.
     */

    Eigen::VectorXd k(v.size());

    for (int i = 0; i < v.size(); ++i) {
        k(i) = k_min + (k_max - k_min) * std::pow(v(i), p);
    }

    return k;
}


double objective(const Eigen::VectorXd& v, const std::vector<std::vector<int>>& rectangles,
    const Eigen::VectorXd& T, const Eigen::Matrix4d& K0,
    double k_min, double k_max, double p) {
    /**
    * Computes the objective function value for a given material distribution.
    *
    * @param v: Density vector (values between 0 and 1).
    * @param rectangles: List of elements, each defined by 4 node indices.
    * @param T: Temperature values at each node.
    * @param K0: Base stiffness matrix (4x4) stored as a sparse matrix.
    * @param k_min: Minimum stiffness value.
    * @param k_max: Maximum stiffness value.
    * @param p: Penalization factor.
    * @return Computed objective function value.
    */

    double D = 0.0;
    Eigen::VectorXd k_values = fill_in_k(v, k_max, k_min, p);

    for (size_t e = 0; e < rectangles.size(); ++e) {
        double k_e = k_values(e);

        // Extraction de T_loc (4x1)
        Eigen::VectorXd T_loc(4);
        for (int l = 0; l < 4; ++l) {
            T_loc(l) = T(rectangles[e][l]);
        }

        // Calcul du produit K * T_loc en utilisant SparseMatrix
        Eigen::VectorXd K_T_loc = k_e * K0 * T_loc;

        // Mise à jour de D avec le produit scalaire (T_loc^T * K * T_loc)
        D += T_loc.transpose() * K_T_loc;
    }

    return D;
}


Eigen::VectorXd adjoint(const Eigen::VectorXd& T, const Eigen::VectorXd& v,
    const std::vector<std::vector<int>>& corners,
    const Eigen::MatrixXd& K0, double p = 3) {
    /**
    * Compute the derivative of the cost function with respect to the fraction of metal for each element.
    *
    * @param T: The global temperature vector.
    * @param v: The global vector with fractions of metal in each element.
    * @param corners: A nested list with first index as elements and second index as the nodes of each element.
    * @param K0: The local stiffness matrix (4x4).
    * @param p: Penalization factor that pushes fractions of metal towards 0 or 1 (p > 1).
    * @return gradJv: A column vector of size(n) containing the derivative of the cost function.
    */

    int elements = v.size();
    Eigen::VectorXd gradJv(elements);
    double km = 65;
    double kp = 0.2;

    for (int el = 0; el < elements; ++el) {
        const std::vector<int>& elementNodes = corners[el];
        Eigen::VectorXd Te(4); // Assuming each element has 4 nodes
        for (int i = 0; i < 4; ++i) {
            Te(i) = T(elementNodes[i]);
        }

        // Compute the derivative (using Eigen matrix and vector operations)
        gradJv(el) = -0.5 * p * std::pow(v(el), p - 1) * (km - kp) * (Te.transpose() * K0 * Te)(0, 0);
    }

    return gradJv;
}
#endif