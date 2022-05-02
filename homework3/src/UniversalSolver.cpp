#include "UniversalSolver.h"

#include <iostream>

namespace homework {

void UniversalSolver::init_values(const Eigen::VectorXd &init_vals, double init_time) {
    values_ = init_vals;
    time_ = init_time;
}

Eigen::MatrixXd Walker::gen_tmp() {
    Eigen::MatrixXd tmp(solver_.a_.cols(), solver_.values_.size());
    Eigen::VectorXd val = solver_.values_;

    tmp.row(0) = solver_.recalc_function_(solver_.time_, val);

    for(int i = 1; i < solver_.a_.rows(); ++i) {
        Eigen::VectorXd val = solver_.values_;
        for(int j = 0; j < i; ++j) {
            val += tmp.row(j) * solver_.a_(i, j) * step_;
        }
        tmp.row(i) = solver_.recalc_function_(solver_.time_ + solver_.steps_(i) * step_, val);
    }

    return tmp;
}

void DPWalker::calc_step() {
    Eigen::VectorXd x1;
    Eigen::VectorXd x2;
    Eigen::MatrixXd tmp;
    double diff = 1;

    do {
        tmp = gen_tmp();

        x1 = solver_.b_[0] * tmp;
        x2 = solver_.b_[1] * tmp;

        diff = (x1 - x2).cwiseAbs().maxCoeff();
        if(diff > max_diff_) {
            step_ /= 2;
        } else if(diff < min_diff_) {
            step_ *= 2;
        }
    } while (diff > max_diff_);

    solver_.values_ += x1 * step_;
    solver_.time_ += step_;
}

void RungeWalker::calc_step() {
    Eigen::MatrixXd tmp = gen_tmp();
    Eigen::VectorXd x = solver_.b_[0] * tmp;
    solver_.values_ = solver_.values_ + x * step_;
    solver_.time_ += step_;
}

UniversalSolver::UniversalSolver(const Eigen::MatrixXd &butcher_matrix,
                                 std::function<Eigen::VectorXd(double, const Eigen::VectorXd &)> recalc_function,
                                 double max_diff, double min_diff)
        : recalc_function_(recalc_function) {
    a_ = butcher_matrix.block(0, 1, butcher_matrix.cols() - 1, butcher_matrix.cols() - 1);
    steps_ = butcher_matrix.block(0, 0, butcher_matrix.cols() - 1, 1);

    if(butcher_matrix.rows() != butcher_matrix.cols()) {
        walker_ = std::make_unique<DPWalker>(max_diff, min_diff, *this, 0.1);
        b_.push_back(butcher_matrix.block(butcher_matrix.rows() - 2, 1, 1, butcher_matrix.cols() - 1));
        b_.push_back(butcher_matrix.block(butcher_matrix.rows() - 1, 1, 1, butcher_matrix.cols() - 1));
    } else {
        walker_ = std::make_unique<RungeWalker>(*this, 0.1);
        b_.push_back(butcher_matrix.block(butcher_matrix.rows() - 1, 1, 1, butcher_matrix.cols() - 1));
    }
}

Eigen::MatrixXd get_runge() {
    Eigen::MatrixXd butcher_matrix(5, 5);
    butcher_matrix << 0,     0,     0,     0,     0,
                      1.0/2, 1.0/2, 0,     0,     0,
                      1.0/2, 0,     1.0/2, 0,     0,
                      1.0,   0,     0,     1.0,   0,
                      0,     1.0/6, 1.0/3, 1.0/3, 1.0/6;
    return butcher_matrix;
}

Eigen::MatrixXd get_DP() {
    Eigen::MatrixXd butcher_matrix(9, 8);
    butcher_matrix << 0,      0,            0,             0,            0,            0,               0,          0,
                      1.0/5,  1.0/5,        0,             0,            0,            0,               0,          0,
                      3.0/10, 3.0/40,       9.0/40,        0,            0,            0,               0,          0,
                      4.0/5,  44.0/45,      -56.0/15,      32.0/9,       0,            0,               0,          0,
                      8.0/9,  19372.0/6561, -25360.0/2187, 64448.0/6561, -212.0/729,   0,               0,          0,
                      1,      9017.0/3168,  -355.0/33,     46732.0/5247, 49.0/176,     -5103.0/18656,   0,          0,
                      1,      35.0/384,     0,             500.0/1113,   125.0/192,    -2187.0/6784,    11.0/84,    0,
                      0,      35.0/384,     0,             500.0/1113,   125.0/192,    -2187.0/6784,    11.0/84,    0,
                      0,      5179.0/57600, 0,             7571.0/16695, 393.0/640,    -92097.0/339200, 187.0/2100, 1.0/40;
    return butcher_matrix;
}
} // end namespace homework