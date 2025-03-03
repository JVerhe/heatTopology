#ifndef optHelper_hpp
#define optHelper_hpp

#include <cassert>
#include <Eigen/Dense>
#include <map>
#include <cmath>

using namespace Eigen;

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
void sparse_H_setup(
    const int nx,const int ny, 
    const float rmin, 
    std::vector<int>& iH, 
    std::vector<int>& jH, 
    std::vector<double>& sH
){
    for (int i1 = 0; i1 < nx; ++i1) {
        for (int j1 = 0; j1 < ny; ++j1) {
            int e1 = i1 * ny + j1;
            for (int i2 = std::max(i1 - static_cast<int>(std::ceil(rmin)) + 1, 0); i2 < std::min(i1 + static_cast<int>(std::ceil(rmin)), nx); ++i2) {
                for (int j2 = std::max(j1 - static_cast<int>(std::ceil(rmin)) + 1, 0); j2 < std::min(j1 + static_cast<int>(std::ceil(rmin)), ny); ++j2) {
                    int e2 = i2 * ny + j2;
                    double weight = std::max(0.0, rmin - std::sqrt((i1 - i2) * (i1 - i2) + (j1 - j2) * (j1 - j2)));
                    iH.push_back(e1);
                    jH.push_back(e2);
                    sH.push_back(weight);
                }
            }
        }
    }
}

void create_sparse_matrix(
    const int nx, const int ny, 
    const float rmin, 
    SparseMatrix<double>& H, 
    VectorXd& Hs
){
    vector<int> iH,jH; 
    vector<double> sH;
    sparse_H_setup(nx,ny,rmin,iH,jH,sH);

    int N = nx * ny;
    vector<Triplet<double>> tripletList;

    for (size_t k = 0; k < sH.size(); ++k) {
        tripletList.emplace_back(iH[k], jH[k], sH[k]);
    }

    H.resize(N, N);
    H.setFromTriplets(tripletList.begin(), tripletList.end());

    // Compute sum of each row
    Hs = VectorXd::Zero(N);
    for (int k = 0; k < H.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(H, k); it; ++it) {
            Hs(it.row()) += it.value();
        }
    }
}

void solve_sparse_lin_sys(const SparseMatrix<double>& K, const VectorXd& F, VectorXd& U){
    ConjugateGradient<SparseMatrix<double>, Lower | Upper> solver;

    solver.compute(K);
    solver.setMaxIterations(1000);
    solver.setTolerance(1e-4);

    U = solver.solve(F);
}


void update_densities(
    const VectorXd& x,
    VectorXd& xnew,
    const VectorXd& B,
    const float move
){
    VectorXd Bx = B.array() * x.array();
    for (int e = 0; e < x.size(); ++e) {
        if (Bx[e] <= std::max(0.0, x[e] - move)) {
            xnew[e] = std::max(0.0, x[e] - move);
        }
        else if (Bx[e] >= std::min(1.0, x[e] + move)) {
            xnew[e] = std::min(1.0, x[e] + move);
        }
        else {
            xnew[e] = Bx[e];
        }
    }
}


void find_new_densities(
    const int nx, const int ny,
    const VectorXd& x,
    VectorXd& xPhys,
    VectorXd& xnew,
    const VectorXd& dc,
    const VectorXd& dv,
    const SparseMatrix<double>& H, 
    const VectorXd& Hs,
    const int ft,
    const float max_vol_frac
){
    double l1 = 0;
    double l2 = 1e9;
    double move = 0.4;
    while ((l2 - l1) / (l2 + l1) > 0.001) {

        double lmid = 0.5 * (l1 + l2);
        VectorXd B = (-dc.array() / (dv.array() * lmid)).max(0).sqrt();
        VectorXd Bx = B.array() * x.array();

        //update_x
        xnew = Eigen::VectorXd::Zero(x.size());
        for (int e = 0; e < x.size(); ++e) {
            if (Bx[e] <= std::max(0.0, x[e] - move)) {
                xnew[e] = std::max(0.0, x[e] - move);
            }
            else if (Bx[e] >= std::min(1.0, x[e] + move)) {
                xnew[e] = std::min(1.0, x[e] + move);
            }
            else {
                xnew[e] = Bx[e];
            }
        }

        update_densities(x,xnew,B,move);

        if (ft == 0 || ft == 1) {
            xPhys = xnew;
        }
        else if (ft == 2) {
            xPhys = H * xnew;
            xPhys = xPhys.cwiseQuotient(Hs);
        }

        if (xPhys.sum() > max_vol_frac * nx * nx) {
            l1 = lmid;
        }
        else {
            l2 = lmid;
        }
    }
}
#endif

