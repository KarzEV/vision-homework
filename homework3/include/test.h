#pragma once

#include <Eigen/Core>

namespace homework {
// variant 7
const double a = 2.5;
const double b = 2.1;
const double d = 0;

Eigen::VectorXd my_sin(double time, const Eigen::VectorXd& val);

[[maybe_unused]] Eigen::VectorXd simple1(double time, const Eigen::VectorXd& val);
[[maybe_unused]] Eigen::VectorXd simple2(double time, const Eigen::VectorXd& val);

[[maybe_unused]] double result1(double time);

void test1();
void test2();
} // end namespace homework