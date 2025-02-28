#ifndef meshHelper_hpp
#define meshHelper_hpp

#include <Eigen/Sparse>

/**
 * Creates a list of rectangles, where each rectangle is defined by the indices of its four nodes.
 * The output is a vector of size (number_of_rectangles × 4), where each row represents a rectangle.
 *
 * @param number_of_points The number of mesh points in one fourth of the mesh.
 * @return A vector containing lists of four node indices forming each rectangle.
 */
std::vector<std::vector<int>> create_rectangle_and_mesh(int number_of_points) {

    Eigen::MatrixXi mesh = Eigen::MatrixXi::Zero(number_of_points, number_of_points);

    int count = 0;
    for (int i = 0; i < number_of_points; ++i) {
        for (int j = 0; j < number_of_points; ++j) {
            mesh(i, j) = count++;
        }
    }

    std::vector<std::vector<int>> rectangles;
    for (int i = 1; i < number_of_points; ++i) {
        for (int j = 0; j < number_of_points - 1; ++j) {
            std::vector<int> line = {
                mesh(i, j),
                mesh(i, j + 1),
                mesh(i - 1, j + 1),
                mesh(i - 1, j)
            };
            rectangles.push_back(line);
        }
    }

    return rectangles;
}

/**
 * Generates a list of coordinates for mesh points in one-fourth of the domain.
 * The output is a matrix of size (number_of_points² × 2), where each row represents a (x, y) coordinate.
 *
 * @param L The length of the chip.
 * @param number_of_points The number of mesh points in one fourth of the mesh.
 * @return An Eigen matrix of size (number_of_points² × 2) containing the coordinates.
 */
Eigen::MatrixXd create_coordinates(double L, int number_of_points) {
    Eigen::MatrixXd coordinates(number_of_points * number_of_points, 2);
    double h = (L / 2.0) / (number_of_points - 1);

    for (int i = 0; i < number_of_points; ++i) {
        for (int j = 0; j < number_of_points; ++j) {
            coordinates(i * number_of_points + j, 0) = j * h; // x-coordinate
            coordinates(i * number_of_points + j, 1) = i * h; // y-coordinate
        }
    }

    return coordinates;
}

std::vector<Eigen::Vector3d> filter_boundary_points_with_index(const Eigen::MatrixXd& coordinates, double L, double lower = 0.003, double upper = 0.007) {
    /**
     * Creates a list of boundary points with their indices.
     * Each boundary point is stored as a 3D vector (index, x, y).
     *
     * @param coordinates: Matrix (N x 2) containing the (x, y) coordinates.
     * @param L: The length of the chip.
     * @param lower: Lower bound for y values.
     * @param upper: Upper bound for y values.
     * @return A vector of Eigen::Vector3d where each entry contains (index, x, y).
     */

    std::vector<Eigen::Vector3d> boundary_points;

    for (int idx = 0; idx < coordinates.rows(); ++idx) {
        double x = coordinates(idx, 0);
        double y = coordinates(idx, 1);

        if ((x == 0 || x == L) && (lower <= y && y <= upper)) {
            boundary_points.emplace_back(idx, x, y);
        }
    }

    return boundary_points;
}

#endif