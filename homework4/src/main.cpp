#include <string>
#include <iostream>

#include "opencv2/features2d.hpp"
#include "opencv2/xfeatures2d.hpp"

#include "frame_matching.h"

const std::string PATH_TO_VIDEO("../data/sample_mpg.avi");

enum Algorithm {
    SURF = 1,
    SIFT = 2,
    BRISK = 3,
};

Algorithm start_menu() {
    std::cout << "Вы можете выбрать следующие алгоритмы:" << std::endl
              << "1 - SURF" << std::endl
              << "2 - SIFT" << std::endl
              << "3 - BRISK" << std::endl;

    int number_algorithm = 0;
    std::cin >> number_algorithm;

    return static_cast<Algorithm>(number_algorithm);
}

bool menu_adaptive_alignment() {
    std::cout << "Использовать адаптивное выравнивание?" << std::endl
              << "1 - Да" << std::endl
              << "2 - Нет" << std::endl;
    int use = 0;
    std::cin >> use;

    if(use == 1) {
        return true;
    } else if(use == 2) {
        return false;
    } else {
        std::cout << "Выравнивание не используется" << std::endl;
        return false;
    }
}

int main() {
    auto number_algorithm = start_menu();
    bool using_adaptive_alignment = menu_adaptive_alignment();

    switch (number_algorithm) {
        case Algorithm::BRISK:
            frame_matching<cv::BRISK>(PATH_TO_VIDEO, using_adaptive_alignment);
            break;
        case Algorithm::SIFT:
            frame_matching<cv::xfeatures2d::SIFT>(PATH_TO_VIDEO, using_adaptive_alignment);
            break;
        case Algorithm::SURF:
            frame_matching<cv::xfeatures2d::SURF>(PATH_TO_VIDEO, using_adaptive_alignment);
            break;
        default:
            throw std::invalid_argument("incorrect number of algorithm: " + std::to_string(number_algorithm));
    }

    return 0;
}