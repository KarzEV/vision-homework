#pragma once

#include <functional>
#include <vector>
#include <memory>

#include <Eigen/Core>

namespace homework {

class UniversalSolver;

class Walker {
public:
    explicit Walker(UniversalSolver &solver, double step = 0.1) : solver_(solver), step_(step) {}

    virtual void calc_step() = 0;

    Eigen::MatrixXd gen_tmp();

    double get_step() const { return step_; }
    void set_step(double step) { step_ = step; }

protected:
    UniversalSolver &solver_;
    double step_;
};

class DPWalker : public Walker {
public:
    DPWalker(double max_diff, double min_diff,
             UniversalSolver &solver, double step = 0.1) : Walker(solver, step),
                                                           max_diff_(max_diff),
                                                           min_diff_(min_diff) {}

    void calc_step() override;

private:
    double max_diff_;
    double min_diff_;
};

class RungeWalker : public Walker {
public:
    using Walker::Walker;

    void calc_step() override;
};

class UniversalSolver {
public:
    friend class Walker;
    friend class DPWalker;
    friend class RungeWalker;

    UniversalSolver(const Eigen::MatrixXd &butcher_matrix,
                    std::function<Eigen::VectorXd(double, const Eigen::VectorXd &)> recalc_function,
                    double max_diff = 10e-10, double min_diff = 10e-10);

    void init_values(const Eigen::VectorXd &init_vals, double init_time);

    void calc_step() { walker_->calc_step(); }
    [[nodiscard]] double t() const { return time_; }
    [[nodiscard]] Eigen::VectorXd& vals() { return values_; }
    [[nodiscard]] double get_step() const { return walker_->get_step(); }

    void set_step(double step) { walker_->set_step(step); }

private:
    std::unique_ptr<Walker> walker_;

    std::function<Eigen::VectorXd(double, const Eigen::VectorXd &)> recalc_function_;
    Eigen::VectorXd values_;
    Eigen::VectorXd steps_;

    Eigen::MatrixXd a_;
    std::vector<Eigen::MatrixXd> b_;

    double time_ = 0;
};

[[maybe_unused]] Eigen::MatrixXd get_runge();
[[maybe_unused]] Eigen::MatrixXd get_DP();
} // end namespace homework
