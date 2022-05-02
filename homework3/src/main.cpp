#include <iostream>
#include <functional>

#include <Eigen/Core>

#include "UniversalSolver.h"
#include "test.h"

void exercise1() {
    std::cout << "exercise1:" << std::endl;

    homework::UniversalSolver solver_runge(homework::get_runge(), homework::simple1);
    Eigen::VectorXd init_val(1);
    init_val(0) = homework::d;
    solver_runge.init_values(init_val, 0);

    double max_diff = 0;
    while (solver_runge.t() < 10) {
        solver_runge.calc_step();
        double current_diff = std::abs(solver_runge.vals()(0, 0) - homework::result1(solver_runge.t()));

        if (current_diff > max_diff) {
            max_diff = current_diff;
        }

        if (current_diff > 0.001) {
            solver_runge.set_step(solver_runge.get_step() / 2);
            Eigen::VectorXd new_vals(1, 1);
            new_vals(0, 0) = homework::result1(solver_runge.t());
            max_diff = 0;
            solver_runge.vals() = new_vals;
        }
    }
    std::cout << "dif: " << max_diff << std::endl;
    std::cout << "step: " << solver_runge.get_step() << std::endl;

    homework::UniversalSolver solver_DP(homework::get_DP(), homework::simple1);
    solver_DP.init_values(init_val, 0);

    double min_step = 0.1;
    int total_steps = 0;
    while (solver_DP.t() < 10) {
        ++total_steps;
        solver_DP.calc_step();

        if(solver_DP.get_step() < min_step) {
            min_step = solver_DP.get_step();
        }
    }

    std::cout << "min step: " << min_step << std::endl;
    std::cout << "total steps: " << total_steps << std::endl;
}

int main() {
    using namespace homework;
    exercise1();
    //exercise2();
    return 0;
}
