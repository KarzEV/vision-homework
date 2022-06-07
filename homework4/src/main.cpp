#include <string>
#include <iostream>

#include "opencv2/features2d.hpp"
#include "opencv2/xfeatures2d.hpp"

#include "frame_matching.h"

const std::string PATH_TO_VIDEO("../data/sample_mpg.avi");

int main() {
    frame_matching<cv::xfeatures2d::SURF>(PATH_TO_VIDEO);
    frame_matching<cv::xfeatures2d::SIFT>(PATH_TO_VIDEO);
    frame_matching<cv::BRISK>(PATH_TO_VIDEO);
    return 0;
}