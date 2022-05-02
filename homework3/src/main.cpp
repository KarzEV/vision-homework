#include <iostream>

#include <Eigen/Core>

#include "UniversalSolver.h"

Eigen::VectorXd my_sin(double time, const Eigen::VectorXd& val) {
    Eigen::VectorXd ret(1);
    ret(0) = -sin(time);
    return ret;
}

void test1() {
    std::cout << "Test1:" << std::endl;
    Eigen::MatrixXd butcher_matrix(5, 5);
    butcher_matrix << 0,     0,     0,     0,     0,
                      1.0/2, 1.0/2, 0,     0,     0,
                      1.0/2, 0,     1.0/2, 0,     0,
                      1.0,   0,     0,     1.0,   0,
                      0,     1.0/6, 1.0/3, 1.0/3, 1.0/6;

    homework::UniversalSolver solver(butcher_matrix, my_sin);
    Eigen::MatrixXd init_val(1, 1);
    init_val(0, 0) = 1;
    solver.init_values(init_val, 0);

    while (solver.t() < M_PI) {
        solver.calc_step();
        std::cout << "cos: " << cos(solver.t()) << "x: "<< solver.vals() << std::endl;
    }
}

void test2() {
    std::cout << "Test2:" << std::endl;
    Eigen::MatrixXd butcher_matrix(9, 8);
    butcher_matrix << 0,    0,            0,             0,            0,            0,               0,          0,
                    1.0/5,  1.0/5,        0,             0,            0,            0,               0,          0,
                    3.0/10, 3.0/40,       9.0/40,        0,            0,            0,               0,          0,
                    4.0/5,  44.0/45,      -56.0/15,      32.0/9,       0,            0,               0,          0,
                    8.0/9,  19372.0/6561, -25360.0/2187, 64448.0/6561, -212.0/729,   0,               0,          0,
                    1,      9017.0/3168,  -355.0/33,     46732.0/5247, 49.0/176,     -5103.0/18656,   0,          0,
                    1,      35.0/384,     0,             500.0/1113,   125.0/192,    -2187.0/6784,    11.0/84,    0,
                    0,      35.0/384,     0,             500.0/1113,   125.0/192,    -2187.0/6784,    11.0/84,    0,
                    0,      5179.0/57600, 0,             7571.0/16695, 393.0/640,    -92097.0/339200, 187.0/2100, 1.0/40;

    homework::UniversalSolver solver(butcher_matrix, my_sin);
    Eigen::MatrixXd init_val(1, 1);
    init_val(0, 0) = 1;
    solver.init_values(init_val, 0);

    while (solver.t() < M_PI) {
        solver.calc_step();
        std::cout << "cos: " << cos(solver.t()) << "x: "<< solver.vals() << std::endl;
    }
}


int main() {
    //test1();
    test2();
    return 0;
}
