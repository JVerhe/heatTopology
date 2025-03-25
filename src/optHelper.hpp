#ifndef optHelper_hpp
#define optHelper_hpp

#include <cassert>
#include <Eigen/Sparse>
#include <map>
#include <cmath>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace Eigen;

/**
     * Computes the material property values (k) based on the density vector (v).
     *
     * @param k: contains the material property values, to be filled. 
     * @param v: Vector of densities (values between 0 and 1).
     * @param k_max: Maximum stiffness value.
     * @param k_min: Minimum stiffness value.
     * @param p: Penalization factor.
     * @return Eigen::VectorXd containing computed k values.
*/
void fill_in_k(VectorXd& k, const VectorXd& v, const double k_max, const double k_min, const double p) {

    assert(k.size()==v.size());

    for (int i = 0; i < v.size(); ++i) {
        k(i) = k_min + (k_max - k_min) * std::pow(v(i), p);
    }
}

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
double objective(const Eigen::VectorXd& v, const std::vector<std::vector<int>>& rectangles,
    const Eigen::VectorXd& T, const Eigen::Matrix4d& K0,
    double k_min, double k_max, double p) {

    double D = 0.0;
    Eigen::VectorXd k_values = VectorXd::Zero(v.size());
    fill_in_k(k_values, v, k_max, k_min, p);

    assert(v.size() == rectangles.size());

    for (size_t e = 0; e < rectangles.size(); ++e) {
        double k_e = k_values(e);

        // Extraction of T_loc (4x1)
        Eigen::VectorXd T_loc(4);
        for (int l = 0; l < 4; ++l) {
            T_loc(l) = T(rectangles[e][l]);
        }

        // Calculation of K * T_loc
        Eigen::VectorXd K_T_loc = k_e * K0 * T_loc;

        // (T_loc^T * K * T_loc)
        D += T_loc.transpose() * K_T_loc;
    }

    return D;
}

/**
    * Compute the derivative of the cost function with respect to the fraction of metal for each element.
    *
    * @param gradJv: A column vector of size(n) containing the derivative of the cost function.
    * @param T: The global temperature vector.
    * @param v: The global vector with fractions of metal in each element.
    * @param corners: A nested list with first index as elements and second index as the nodes of each element.
    * @param K0: The local stiffness matrix (4x4).
    * @param p: Penalization factor that pushes fractions of metal towards 0 or 1 (p > 1).
    * @param kmin: minimum conductivity. 
    * @param kmax: maximum conductivity. 
    * @return void
*/
void adjoint(Eigen::VectorXd& gradJv, const Eigen::VectorXd& T, const Eigen::VectorXd& v,
    const std::vector<std::vector<int>>& corners,
    const Eigen::MatrixXd& K0, double p = 3, double kmin = 0.2, double kmax = 65) {

    int elements = v.size();

    for (int el = 0; el < elements; ++el) {
        const std::vector<int>& elementNodes = corners[el];
        Eigen::VectorXd Te(4); // Assuming each element has 4 nodes
        for (int i = 0; i < 4; ++i) {
            Te(i) = T(elementNodes[i]);
        }

        // Compute the derivative (using Eigen matrix and vector operations)
        gradJv(el) = -0.5 * p * std::pow(v(el), p - 1) * (kmax - kmin) * (Te.transpose() * K0 * Te)(0, 0);
    }
}


/**
 * @brief construct sparse matrix H via vectors iH, jH and sH, s.t. H(iH[k],jH[k]) == sH[k]
 * 
 * @param nx amount of elements in x-direction
 * @param ny amount of elements in y-direction
 * @param rmin filter radius
 * @param iH vector containing row-indices of Hcorresponding to non-zero values of H
 * @param jH vector containing column-indices of H corresponding to non-zero values of H
 * @param sH vector containing weights of H, s.t. H(iH[k],jH[k]) == sH[k]
 * 
 * @return void
 */
void sparse_H_setup(
    const int nx, const int ny,
    const float rmin,
    std::vector<int>& iH,
    std::vector<int>& jH,
    std::vector<double>& sH
) {
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

/**
 * @brief set up sparse matrix H used for sensitivity or density filtering.
 * 
 * @param nx amount of elements in x-direction
 * @param ny amount of elements in y-direction
 * @param rmin filter radius
 * @param H sparse matrix to be set up
 * @param Hs row-sum vector of H
 * 
 * @return void
 */
void create_sparse_matrix(
    const int nx, const int ny,
    const float rmin,
    SparseMatrix<double>& H,
    VectorXd& Hs
) {
    std::vector<int> iH, jH;
    std::vector<double> sH;
    sparse_H_setup(nx, ny, rmin, iH, jH, sH);

    int N = nx * ny;
    std::vector<Triplet<double>> tripletList;

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

/**
 * @brief Solves the linear system \( K T = F \) using the Conjugate Gradient (CG) solver.
 * 
 * In this approach, we solve the system of equations \( K T = F \) iteratively using the Conjugate Gradient method. 
 * The Conjugate Gradient method is well-suited for large, sparse, symmetric, and positive-definite matrices such as the stiffness matrix \( K \).
 * This solver is particularly efficient when the matrix \( K \) is large, and the solution vector \( T \) is sparse. The method avoids direct matrix inversion 
 * and instead iteratively improves the solution by minimizing the residual error.
 * 
 * The boundary conditions are applied to the solution vector \( T \) by modifying the right-hand side vector \( F \) and adjusting the entries of \( T \) at the boundary nodes 
 * before running the CG solver. This approach is more flexible when the boundary conditions are complex and need to be enforced at specific grid points.
 * 
 * @param K The stiffness matrix representing the system.
 * @param F The right-hand side vector containing the forces or source terms.
 * @param T The solution vector, which will contain the temperature or displacement values after solving.
 * 
 * @return void
 */
void solve_sparse_lin_sys(const SparseMatrix<double>& K, const VectorXd& F, VectorXd& U) {
    ConjugateGradient<SparseMatrix<double>, Lower | Upper> solver;

    solver.compute(K);
    solver.setMaxIterations(100000000);
    solver.setTolerance(1e-15);

    U = solver.solve(F);
}


/**
 * @brief Solves the linear system \( K T = F \) using the LU decomposition solver.
 * 
 * This approach uses LU decomposition to solve the system \( K T = F \) directly. LU decomposition is a standard method for solving linear systems 
 * and is particularly useful when the matrix \( K \) is dense or when we have direct access to the full matrix. In this method, we decompose 
 * the matrix \( K \) into a lower triangular matrix \( L \) and an upper triangular matrix \( U \), and then solve for \( T \) using forward 
 * and backward substitution.
 * 
 * LU decomposition is more compatible with the way we impose boundary conditions in this context. The boundary conditions can be enforced 
 * more directly by modifying the matrix \( K \) itself before solving the system, ensuring that the solution vector \( T \) satisfies the boundary 
 * conditions at the prescribed locations. This is often more straightforward for systems with Dirichlet boundary conditions, where the solution 
 * at the boundaries is fixed.
 * 
 * @param K The stiffness matrix representing the system.
 * @param F The right-hand side vector containing the forces or source terms.
 * @param T The solution vector, which will contain the temperature or displacement values after solving.
 * 
 * @return void
 */
void solve_sparse_lin_sys_LU(const SparseMatrix<double>& K, const VectorXd& F, VectorXd& U) {

    SparseLU<SparseMatrix<double>> solver;

    solver.analyzePattern(K);
    solver.factorize(K);

    U = solver.solve(F);
}

/**
 * @brief Updates densities according to heuristic optimization scheme as implemented in the reference paper.
 *
 * @param x Old densities.
 * @param xnew New densities.
 * @param Bx Vector needed for updating scheme: B@x
 * @param move Maximum change per iteration per element.
 *
 * @return void
 */
void update_densities(
    const VectorXd& x,
    VectorXd& xnew,
    const VectorXd& B,
    const float move
) {
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

/**
 * @brief Finds Lagrange parameter and updates densities according to heuristic optimization scheme as implemented in the reference paper.
 * 
 * @param nx amount of elements in x-direction
 * @param ny amount of elements in y-direction
 * @param x Vector containing densities of previous solution 
 * @param xPhys Vector containing densities, updated according to filtering scheme. 
 * @param xnew Vector containing the densities after update. 
 * @param dc gradient of objective function. 
 * @param dv gradient of volume of elements, should be a vector of ones. 
 * @param H matrix for filtering. 
 * @param Hs vector for filtering. 
 * @param ft indicates which filtering is used: 0=no filtering, 1=sensitivity filtering, 2=density filtering. 
 * @param max_vol_frac maximum fraction of metal in plate. 
 */
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
) {
    double l1 = 0;
    double l2 = 1e9;
    double move = 0.1;
    while ((l2 - l1) / (l2 + l1) > 0.001) {

        double lmid = 0.5 * (l1 + l2);
        VectorXd B = (-dc.array() / (dv.array() * lmid)).max(0).sqrt();

        update_densities(x, xnew, B, move);

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

