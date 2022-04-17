#include "my_matrix.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>

#define assertm(exp, msg) assert(((void)msg, exp))

std::vector<homework::MyMatrix> get_test_data(std::string path) {
    homework::ReaderWayPointsFromCSV reader;
    std::vector<homework::MyMatrix> test_data;

    std::fstream fs(path);

    for (size_t i = 0; i < 4; ++i) {
        test_data.push_back(reader.parse(fs));
    }

    return test_data;
}

void test(std::string path) {
    std::cout << "start test: " << path << std::endl;

    using namespace std::string_literals;
    std::string error_message("failed test: "s + path);
    std::vector<homework::MyMatrix> test_data = get_test_data(path);

    assertm(homework::solve(test_data[0], test_data[1]) == test_data[3], error_message);
    assertm(test_data[0].inverse() == test_data[2], error_message);
}

int main() {
//    test("../test_data/test_1.csv");
//    test("../test_data/test_2.csv");
    //test("../test_data/test_3.csv");
    test("../test_data/test_4.csv");
    test("../test_data/test_5.csv");
    test("../test_data/test_6.csv");

    std::cout << "tests passed" << std::endl;
    return 0;
}

