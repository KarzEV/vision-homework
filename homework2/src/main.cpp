#include <fstream>
#include <iostream>
#include <math.h>

#include <Eigen/Core>
#include <Eigen/LU>

template<typename M>
M load_csv(const std::string &path) {
    std::ifstream file(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;

    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }

    return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, Eigen::RowMajor>>(
            values.data(), rows, values.size() / rows);
}

Eigen::MatrixXd mEstimator(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B) {
    using namespace Eigen;
    MatrixXd Weight = MatrixXd::Identity(A.rows(), A.rows());
    double delta = 0.01;

    MatrixXd result;
    for (uint i = 0; i < 50; ++i) {
        result = (A.transpose() * Weight * A).lu().solve(A.transpose() * Weight * B);
        VectorXd res(B - A * result);

        for (uint j = 0; j < res.rows(); j++) {
            Weight(j, j) = 1 / (std::max(delta, fabs(res(j))));
        }
    }

    return result;
}

int main() {
    const std::string path("../src_data/test_data.csv");
    auto data = load_csv<Eigen::MatrixXd>(path);

    Eigen::MatrixXd B = data.col(0);
    for (long i = 0; i < B.rows(); ++i) {
        B(i, 0) = std::log(data(i, 0));
    }
    Eigen::MatrixXd A(data.rows(), 2);

    for (long i = 0; i < A.rows(); ++i) {
        A(i, 0) = data(i, 1);
        A(i, 1) = 1;
    }

    std::cout << mEstimator(A, B) << std::endl;
}