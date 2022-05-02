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
} // end namespace homework