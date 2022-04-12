#include "my_matrix.h"

#include <stdexcept>
#include <utility>

using namespace homework;

MyMatrix::MyMatrix(size_t size) {
    MyMatrix(size, size);
}

MyMatrix::MyMatrix(size_t rows, size_t cols) {
    data_.reserve(rows);

    for(auto& row: data_) {
        row.reserve(cols);
    }
}

MyMatrix MyMatrix::operator*(double val) const {
    MyMatrix result(data_);

    for(auto& one_row: result.data_) {
        for(auto& one_val: one_row) {
            one_val *= val;
        }
    }

    return result;
}

MyMatrix MyMatrix::operator*(const MyMatrix &val) const {
    if(cols_ != val.rows_) {
        throw std::invalid_argument("incorrect size matrix");
    }

    MyMatrix result(rows_, val.cols_);

    for(size_t i = 0; i < rows_; ++i) {
        for(size_t j = 0; j < val.cols_; ++j) {
            for(size_t k = 0; k < cols_; ++k) {
                result[i][j] = data_[i][k] * val[j][k];
            }
        }
    }

    return result;
}

MyMatrix MyMatrix::operator+(const MyMatrix &val) const {
    if(rows_ != val.rows_ || cols_ != val.cols_) {
        throw std::invalid_argument("incorrect size matrix");
    }

    MyMatrix result(data_);

    for(size_t i = 0; i < rows_; ++i) {
        for(size_t j = 0; j < cols_; ++j) {
            result[i][j] += val[i][j];
        }
    }

    return result;
}

MyMatrix MyMatrix::operator-(const MyMatrix &val) const {
    if(rows_ != val.rows_ || cols_ != val.cols_) {
        throw std::invalid_argument("incorrect size matrix");
    }

    MyMatrix result(data_);

    for(size_t i = 0; i < rows_; ++i) {
        for(size_t j = 0; j < cols_; ++j) {
            result[i][j] -= val[i][j];
        }
    }

    return result;

}

bool MyMatrix::operator!=(const MyMatrix &val) const {
    return !(*this == val);
}

bool MyMatrix::operator==(const MyMatrix &val) const {
    if(rows_ != val.rows_ || cols_ != val.cols_) {
        return false;
    }

    for(size_t i = 0; i < rows_; ++i) {
        for(size_t j = 0; j < cols_; ++j) {
            if(data_[i][j] != val[i][j]) {
                return false;
            }
        }
    }

    return true;
}

std::tuple<MyMatrix, MyMatrix, std::vector<int>> MyMatrix::LU() const {
    if(rows_ != cols_) {
        throw std::invalid_argument("LU transformation is incorrect for a non-square matrix");
    }

    MyMatrix U_matrix(*this);
    MyMatrix L_matrix = MyMatrix::gen_identity(rows_, cols_);
    std::vector<int> position_x;

    for(size_t i = 0; i < rows_; ++i) {
        position_x.push_back(i);
    }

    for(size_t i = 0; i < rows_; ++i) {
        int new_position = U_matrix.find_main_element(i, i);
        position_x[i] = new_position;
        position_x[new_position] = i;

        for(size_t j = i + 1; i < rows_; ++j) {
            L_matrix[j][i] = U_matrix[j][i] / U_matrix[i][i];
            for(size_t k = 0; k < cols_; ++k) {
                U_matrix[j][k] -= L_matrix[j][i] * U_matrix[i][k];
            }
        }
    }

    return std::make_tuple(U_matrix, L_matrix, position_x);
}

size_t MyMatrix::find_main_element(size_t row, size_t col) {
    double max = data_[row][col];
    size_t max_row = row;

    for(size_t j = row + 1; row < rows_; ++j) {
        if(data_[j][col] > max) {
            max = data_[j][col];
            max_row = j;
        }
    }

    if(row != max_row) {
        std::vector<double> buf(std::move(data_[row]));
        data_[row] = std::move(data_[max_row]);
        data_[max_row] = std::move(buf);
    }

    return max_row;
}

MyMatrix MyMatrix::gen_identity(size_t rows, size_t cols) {
    MyMatrix result(rows, cols);

    for(size_t i = 0; i < rows; ++i) {
        for(size_t j = 0; j < cols; ++j) {
            if(i == j) {
                result[i][j] = 1;
            } else {
                result[i][j] = 0;
            }
        }
    }

    return result;
}

double MyMatrix::det() const {

}

MyMatrix MyMatrix::inverse() const {

}
