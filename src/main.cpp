#include <Eigen/Sparse>
#include <iostream>
#include <csignal>
#include <cstdlib>
#include "optimization.hpp"
#include "mms.cpp"
#include "meshHelper.hpp"
#include <random>
#include <chrono>



void readConfig(const std::string& filename, int& number_of_points, int& p, int& ft, bool& optimisation, bool& k_constant) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "No test file with this name found" << std::endl;
        return;
    }
    char space[100];
    file >> space >> number_of_points >> space >> p >> space >> ft >> space >> optimisation >> space>> k_constant;
    file.close();
}

void writeOutput(const std::string& filename, Eigen::VectorXd x) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (int i = 0; i < x.size(); ++i) {
            file << x(i) << std::endl;
        }
        file.close();
        std::cout << "Saved to " << filename << std::endl;
    }
    else {
        std::cerr << "Error when saving " << filename << std::endl;
    }
    return;

}

void signalHandler(int signal) {
    std::cerr << "No config file with that name was found" << std::endl;
    std::exit(signal);
}


int main(int argc, char* argv[]) {

    if (argc == 1) {
        std::cerr << "Append the name of a config file when executing main" << std::endl;
        return 1;
    }

    int number_of_points; int p; int ft; bool optimization; bool k_constant;

    std::signal(SIGSEGV, signalHandler);
    std::string file_name = argv[1];
    std::string config_file = "config/" + file_name + ".txt";

    readConfig(config_file, number_of_points, p, ft,optimization,k_constant);

    double  L = 0.01;
    double T_k = 293;
    double vol_frac = 0.4;

    int size = (number_of_points - 1) * (number_of_points - 1);
    std::vector<std::vector<int>> rectangles = create_rectangle_and_mesh(number_of_points);

    Eigen::Matrix4d local_matrix;
    local_matrix << 2.0 / 3, -1.0 / 6, -1.0 / 3, -1.0 / 6,
        -1.0 / 6, 2.0 / 3, -1.0 / 6, -1.0 / 3,
        -1.0 / 3, -1.0 / 6, 2.0 / 3, -1.0 / 6,
        -1.0 / 6, -1.0 / 3, -1.0 / 6, 2.0 / 3;

    Eigen::VectorXd x(size); 

    if(optimization){
        ///////////////////////////// INITIAL RANDOM SOLUTION //////////////////////////////////////
        std::random_device rd;                                                                    //  
        std::mt19937 gen(rd());  // Generator with a random seed                                  //               
        std::uniform_real_distribution<> dis(0.0, 1.0);                                           //
                                                                                                  //                                                                           //
        for (int i = 0; i < size; ++i) {                                                          //
            x(i) = dis(gen) * 0.5 + 0.5;  // Fill the vector with random values between 0 and 1   //
        }                                                                                         //                  
        x = (x.array() - x.mean()) + 0.4;  // Adjust the mean to 0.4                              //
        x = x.cwiseMin(1.0).cwiseMax(0.0); // Ensure all values remain in [0,1]                   //
        ////////////////////////////////////////////////////////////////////////////////////////////
        optimize(local_matrix, x, vol_frac, number_of_points - 1, number_of_points - 1, p, rectangles, L, T_k, ft);
    }else{
        mms(local_matrix, x, vol_frac, number_of_points - 1, number_of_points - 1, p, rectangles, 2*L, T_k, k_constant);
    }

    ////////////////// Call the program to plot the result //////
    std::string callPython = "python3 ../src/plot_result.py";  //
    int rc = system(callPython.c_str());                       //
    /////////////////////////////////////////////////////////////

    return 0;
}