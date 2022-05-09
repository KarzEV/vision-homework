#include <iostream>
#include <functional>

#include <Eigen/Core>

#include "UniversalSolver.h"
#include "test.h"

void equation(std::function<Eigen::VectorXd(double, const Eigen::VectorXd &)> src_function,
              std::function<Eigen::VectorXd(double)> result_function, const Eigen::VectorXd& init_val) {
    homework::UniversalSolver solver_runge(homework::get_runge(), src_function);
    solver_runge.init_values(init_val, 0);

    double max_diff = 0;
    while (solver_runge.t() < 0.2) {
        solver_runge.calc_step();
        double current_diff = (solver_runge.vals() - result_function(solver_runge.t())).cwiseAbs().maxCoeff();

        if (current_diff > max_diff) {
            max_diff = current_diff;
        }

        if (current_diff > 0.001) {
            solver_runge.set_step(solver_runge.get_step() / 2);
            max_diff = 0;
            solver_runge.vals() = result_function(solver_runge.t());
        }
    }
    std::cout << "dif: " << max_diff << std::endl;
    std::cout << "step: " << solver_runge.get_step() << std::endl;

    homework::UniversalSolver solver_DP(homework::get_DP(), src_function);
    solver_DP.init_values(init_val, 0);

    double min_step = 0.1;
    int total_steps = 0;
    while (solver_DP.t() < 10) {
        ++total_steps;
        solver_DP.calc_step();

        if (solver_DP.get_step() < min_step) {
            min_step = solver_DP.get_step();
        }
    }

    std::cout << "min step: " << min_step << std::endl;
    std::cout << "total steps: " << total_steps << std::endl;
}

void equation1() {
    std::cout << "equation1:" << std::endl;
    Eigen::VectorXd init_value(1, 1);
    init_value(0, 0) = homework::d;
    equation(homework::simple1, homework::result1, init_value);
}

void equation2() {
    std::cout << "equation2:" << std::endl;
    Eigen::VectorXd init_value(2, 1);
    init_value(0, 0) = 4.0 / 3;
    init_value(1, 0) = 2.0 / 3;
    equation(homework::simple2, homework::result2, init_value);
}

int main() {
    equation1();
    equation2();
    return 0;
}
