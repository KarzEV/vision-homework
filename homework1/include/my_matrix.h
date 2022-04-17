#pragma once

#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

namespace homework {

    class MyMatrix {
    private:
        class ProxyMatrix {
        public:
            ProxyMatrix(std::vector<double> &col) : col_(col) {}

            double &operator[](size_t j) { return col_[j]; }

        private:
            std::vector<double> &col_;
        };

        class ConstProxyMatrix {
        public:
            ConstProxyMatrix(const std::vector<double> &col) : col_(col) {}

            double operator[](size_t j) const { return col_[j]; }

        private:
            const std::vector<double> &col_;
        };

    public:
        explicit MyMatrix(size_t size) : MyMatrix(size, size) {}

        MyMatrix(size_t rows, size_t cols);

        explicit MyMatrix(const std::vector<std::vector<double>> &data);

        MyMatrix(const MyMatrix &other) = default;

        MyMatrix &operator=(const MyMatrix &other) = default;

        MyMatrix(MyMatrix &&other) = default;

        MyMatrix &operator=(MyMatrix &&other) = default;

    public:
        size_t get_rows() const { return rows_; }

        size_t get_cols() const { return cols_; }

        MyMatrix operator*(double val) const;

        MyMatrix operator*(const MyMatrix &val) const;

        MyMatrix operator+(const MyMatrix &val) const;

        MyMatrix operator-(const MyMatrix &val) const;

        bool operator!=(const MyMatrix &val) const;

        bool operator==(const MyMatrix &val) const;

        ProxyMatrix operator[](size_t i) { return ProxyMatrix(data_[i]); }

        ConstProxyMatrix operator[](size_t i) const { return ConstProxyMatrix(data_[i]); }

        std::tuple<MyMatrix, MyMatrix, std::vector<int>> LU() const;

        double det() const;

        MyMatrix inverse() const;

        static MyMatrix gen_identity(size_t rows, size_t cols);

    private:
        size_t find_main_element(size_t row, size_t col);

    private:
        size_t rows_;
        size_t cols_;
        std::vector<std::vector<double>> data_;
    };

    class ReaderWayPointsFromCSV {
    public:
        typedef boost::char_separator<char> sep_csv;
        typedef boost::tokenizer<sep_csv> tok_type;

        [[nodiscard]] MyMatrix parse(std::fstream &fs);

    private:
        typedef boost::tokenizer<ReaderWayPointsFromCSV::sep_csv>::iterator ptr;

    private:
        std::vector<double> parse_csv_line_(ptr begin, ptr end);

        template<typename ValueT>
        ValueT bindVariable(ptr pos, ptr end);
    };

    std::ostream &operator<<(std::ostream &os, const homework::MyMatrix &matrix);

    MyMatrix solve(const MyMatrix &A, const MyMatrix &b);
} // end namespace homework
