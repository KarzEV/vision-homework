#include "test.h"

#include <iostream>

#include "UniversalSolver.h"

namespace homework {

extern const double a;
extern const double b;
extern const double d;

Eigen::VectorXd my_sin(double time, const Eigen::VectorXd& val) {
    Eigen::VectorXd ret(1);
    ret(0) = -sin(time);
    return ret;
}

void test1() {
    std::cout << "Test1:" << std::endl;

    UniversalSolver solver(get_runge(), my_sin);
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

    UniversalSolver solver(get_DP(), my_sin);
    Eigen::MatrixXd init_val(1, 1);
    init_val(0, 0) = 1;
    solver.init_values(init_val, 0);

    while (solver.t() < M_PI) {
        solver.calc_step();
        std::cout << "cos: " << cos(solver.t()) << "x: "<< solver.vals() << std::endl;
    }
}

Eigen::VectorXd simple1(double time, const Eigen::VectorXd &val) {
    Eigen::VectorXd ret(1);
    ret(0) = a * time - b * val(0, 0);
    return ret;
}

Eigen::VectorXd simple2(double time, const Eigen::VectorXd &val) {
    Eigen::VectorXd ret(2, 1);
    ret(0, 0) = 9 * val(0, 0) + 24 * val(1, 0) + 5 * cos(time) - 1.0 / 3 * sin(time);
    ret(1, 0) = -24 * val(0, 0) - 51 * val(1, 0) + 9 * cos(time) + 1.0 / 3 * sin(time);
    return ret;
}

Eigen::VectorXd result1(double time) {
    Eigen::VectorXd result(1);
    result(0) =  a / b * (time - 1 / b) + 0.57 * exp(-b * time);
    return result;
}

Eigen::VectorXd result2(double time) {
    Eigen::VectorXd result(2);
    result(0) = 2.0 * exp(-3 * time) - exp(-39 * time) + 1.0 / 3 * cos(time);
    result(1) = - exp(-3 * time) + 2.0 * exp(-39 * time) - 1.0 / 3 * cos(time);
    return result;
}
}