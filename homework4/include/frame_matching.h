#pragma once

#include <string>
#include <vector>

#include "opencv2/opencv.hpp"
#include "opencv2/features2d.hpp"

inline bool break_for_key(int delay) {
    int key = cv::waitKey(delay);

    switch (key) {
        case 27:
            return true;
        case 32:
            while(true) {
                if(cv::waitKey(delay) == 32) {
                    break;
                }
            }
        default:
            return false;
    }
}

template<class Detector>
void frame_matching(const std::string& path) {
    using namespace cv;
    cv::Ptr<cv::FeatureDetector> detector = Detector::create();
    cv::Ptr<cv::DescriptorExtractor> extractor = Detector::create();
    cv::BFMatcher matcher;

    std::vector<cv::KeyPoint> keys;
    std::vector<cv::KeyPoint> prev_keys;

    cv::Mat src;
    cv::Mat prev_src;
    cv::Mat prev_deser;
    cv::Mat deser;
    cv::Mat image_with_points;

    cv::VideoCapture cap(path);
    if(!cap.isOpened() || !cap.grab()) {
        throw std::invalid_argument("incorrect path: " + path);
    }

    int delay = static_cast<int>(1000 / cap.get(cv::CAP_PROP_FPS));

    cap.read(prev_src);
    detector->detect(prev_src, prev_keys);
    extractor->compute(prev_src, prev_keys, prev_deser);

    while(true) {
        std::vector<cv::DMatch> matches;

        cap.read(src);

        if(src.empty()) {
            break;
        }

        detector->detect(src, keys);
        extractor->compute(src, keys, deser);

        matcher.match(deser, prev_deser, matches);
        cv::drawMatches(src, keys, prev_src, prev_keys, matches, image_with_points);

        cv::imshow(typeid(Detector).name(), image_with_points);

        prev_src = std::move(src);
        prev_keys = std::move(keys);
        prev_deser = std::move(deser);

        if(break_for_key(delay)) {
            break;
        }
    }
    cap.release();
    cvDestroyWindow(typeid(Detector).name());
}
