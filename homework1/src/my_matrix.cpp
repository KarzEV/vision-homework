#include "my_matrix.h"

#include <math.h>

#include <stdexcept>
#include <utility>

namespace homework {

    MyMatrix::MyMatrix(size_t rows, size_t cols) {
        data_.resize(rows);

        for (auto &row: data_) {
            row.resize(cols, 0);
        }

        rows_ = rows;
        cols_ = cols;
    }

    MyMatrix MyMatrix::operator*(double val) const {
        MyMatrix result(data_);

        for (auto &one_row: result.data_) {
            for (auto &one_val: one_row) {
                one_val *= val;
            }
        }

        return result;
    }

    MyMatrix MyMatrix::operator*(const MyMatrix &val) const {
        if (cols_ != val.rows_) {
            throw std::invalid_argument("incorrect size matrix");
        }

        MyMatrix result(rows_, val.cols_);

        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < val.cols_; ++j) {
                for (size_t k = 0; k < cols_; ++k) {
                    result[i][j] += data_[i][k] * val[k][j];
                }
            }
        }

        return result;
    }

    MyMatrix MyMatrix::operator+(const MyMatrix &val) const {
        if (rows_ != val.rows_ || cols_ != val.cols_) {
            throw std::invalid_argument("incorrect size matrix");
        }

        MyMatrix result(data_);

        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result[i][j] += val[i][j];
            }
        }

        return result;
    }

    MyMatrix MyMatrix::operator-(const MyMatrix &val) const {
        if (rows_ != val.rows_ || cols_ != val.cols_) {
            throw std::invalid_argument("incorrect size matrix");
        }

        MyMatrix result(data_);

        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result[i][j] -= val[i][j];
            }
        }

        return result;

    }

    bool MyMatrix::operator!=(const MyMatrix &val) const {
        return !(*this == val);
    }

    bool MyMatrix::operator==(const MyMatrix &val) const {
        if (rows_ != val.rows_ || cols_ != val.cols_) {
            return false;
        }

        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                if (std::fabs(data_[i][j] - val[i][j]) > 0.001) {
                    return false;
                }
            }
        }

        return true;
    }

    std::tuple<MyMatrix, MyMatrix, std::vector<int>> MyMatrix::LU() const {
        if (rows_ != cols_) {
            throw std::invalid_argument("LU transformation is incorrect for a non-square matrix");
        }

        MyMatrix U_matrix(*this);
        MyMatrix L_matrix = MyMatrix::gen_identity(rows_, cols_);
        std::vector<int> position_x;

        for (size_t i = 0; i < rows_; ++i) {
            position_x.push_back(i);
        }

        for (size_t i = 0; i < rows_; ++i) {
            int new_position = U_matrix.find_main_element(i, i);
            if (std::fabs(U_matrix[new_position][new_position] - 0) <= std::numeric_limits<double>::epsilon()) {
                throw std::invalid_argument("degenerate matrix");
            }

            if(new_position != i) {
                int prev_position = position_x[i];
                position_x[i] = position_x[new_position];
                position_x[new_position] = prev_position;

                for(size_t k = 0; k < i; ++k) {
                    double buf = L_matrix[new_position][k];
                    L_matrix[new_position][k] = L_matrix[i][k];
                    L_matrix[i][k] = buf;
                }
            }

            for (size_t j = i + 1; j < rows_; ++j) {
                L_matrix[j][i] = U_matrix[j][i] / U_matrix[i][i];
                for (size_t k = 0; k < cols_; ++k) {
                    U_matrix[j][k] -= L_matrix[j][i] * U_matrix[i][k];
                }
            }
        }

        return std::make_tuple(U_matrix, L_matrix, position_x);
    }

    size_t MyMatrix::find_main_element(size_t row, size_t col) {
        double max = data_[row][col];
        size_t max_row = row;

        for (size_t j = row + 1; j < rows_; ++j) {
            if (std::fabs(data_[j][col]) > std::fabs(max)) {
                max = data_[j][col];
                max_row = j;
            }
        }

        std::vector<double> buf(std::move(data_[row]));
        data_[row] = std::move(data_[max_row]);
        data_[max_row] = std::move(buf);

        return max_row;
    }

    MyMatrix MyMatrix::gen_identity(size_t rows, size_t cols) {
        MyMatrix result(rows, cols);

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                if (i == j) {
                    result[i][j] = 1;
                } else {
                    result[i][j] = 0;
                }
            }
        }

        return result;
    }

    double MyMatrix::det() const {
        if (rows_ != cols_ || rows_ == 0) {
            throw std::invalid_argument("matrix is not square");
        }

        auto [U, _, __] = LU();

        double det = data_[0][0];
        for (uint32_t i = 1; i < rows_; ++i) {
            det *= U[i][i];
        }

        return det;
    }

    MyMatrix MyMatrix::inverse() const {
        MyMatrix in_matrix(rows_, cols_);

        auto [U, L, position_x] = LU();

        for(size_t i = 0; i < get_rows(); ++i) {
            MyMatrix b(get_cols(), 1);
            b[i][0] = 1;

            b = L.inverse_L() * b;
            MyMatrix x(get_rows(), 1);

            for(int i = get_rows() - 1; i >= 0; --i) {
                for(int j = get_cols() - 1; j >= i; --j) {
                    if(i != j) {
                        x[i][0] -= U[i][j] * x[j][0];
                    } else {
                        x[i][0] = (x[i][0] + b[j][0]) / U[i][j];
                    }
                }
            }

            for(size_t k = 0; k < get_rows(); ++k) {
                in_matrix[k][position_x[i]] = x[k][0];
            }
        }

        return in_matrix;
    }

    MyMatrix::MyMatrix(const std::vector<std::vector<double>> &data) : rows_(data.size()),
                                                                       data_(data) {
        if (rows_ != 0) {
            cols_ = data.at(0).size();
        } else {
            cols_ = 0;
        }
    }

    MyMatrix MyMatrix::inverse_L() {
        MyMatrix invers = MyMatrix::gen_identity(get_rows(), get_cols());

        for(int j = get_cols() - 2; j >= 0; --j) {
            MyMatrix one_M = MyMatrix::gen_identity(get_rows(), get_cols());
            for(int i = get_rows() - 1; i > j; --i) {
                one_M[i][j] = - data_[i][j];
            }
            invers = invers * one_M;
        }

        return invers;
    }

    std::ostream &operator<<(std::ostream &os, const homework::MyMatrix &matrix) {
        for (size_t i = 0; i < matrix.get_rows(); ++i) {
            for (size_t j = 0; j < matrix.get_cols(); ++j) {
                os << matrix[i][j] << ' ';
            }
            os << std::endl;
        }

        return os;
    }

    MyMatrix solve(const MyMatrix &A, const MyMatrix &b) {
        if(b.get_cols() != 1) {
            throw std::invalid_argument("incorrect matrix b");
        }

        MyMatrix result(b.get_rows(), 1);

        auto [U_matrix, L_matrix, position_x] = A.LU();

        MyMatrix right(b.get_rows(), b.get_cols());

        for(size_t i = 0; i < b.get_rows(); ++i) {
            right[i][0] = b[position_x[i]][0];
        }

        right = L_matrix.inverse_L() * right;

        for(int i = A.get_rows() - 1; i >= 0; --i) {
            for(int j = A.get_cols() - 1; j >= i; --j) {
                if(i != j) {
                    result[i][0] -= U_matrix[i][j] * result[j][0];
                } else {
                    result[i][0] = (result[i][0] + right[j][0]) / U_matrix[i][j];
                }
            }
        }

        return result;
    }

    [[maybe_unused]] std::vector<double>
    ReaderWayPointsFromCSV::parse_csv_line_(ptr begin,
                                            ptr end) {
        std::vector<double> result_parse;
        while (begin != end) {
            result_parse.push_back(bindVariable<double>(begin, end));
            begin++;
        }

        return result_parse;
    }

    MyMatrix ReaderWayPointsFromCSV::parse(std::fstream &fs) {
        if (!fs.is_open()) {
            throw std::invalid_argument("can't open input-file csv");
        }

        std::vector<std::vector<double>> values;

        std::string line;
        sep_csv sep(",");

        uint unique_id = 0;

        while (!fs.eof()) {
            std::getline(fs, line);

            if (!line.empty()) {
                tok_type tok(line, sep);

                auto new_point = parse_csv_line_(tok.begin(), tok.end());

                values.push_back(new_point);
                break;
            }
        }

        while (!fs.eof()) {
            std::getline(fs, line);

            if (line.empty()) {
                break;
            }

            tok_type tok(line, sep);

            auto new_point = parse_csv_line_(tok.begin(), tok.end());

            values.push_back(new_point);
        }

        return MyMatrix(values);
    }

    template<typename ValueT>
    ValueT ReaderWayPointsFromCSV::bindVariable(ptr pos, ptr end) {
        if (pos == end) {
            throw std::runtime_error("bad csv format");
        }

        try {
            return boost::lexical_cast<ValueT>(*pos);
        } catch (const boost::bad_lexical_cast &) {
            throw std::runtime_error("invalid convert from csv");
        }
    }

} // end homework
