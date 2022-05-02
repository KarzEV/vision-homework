#include <string>
#include <fstream>

#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

void detect_text(string input){
    Mat original = imread(input);
    imshow("Source", original);

    Mat compressed;
    pyrDown(original, compressed);
    Mat small;
    cvtColor(compressed, small, CV_BGR2GRAY);

    Mat grad;
    Mat morphKernel = getStructuringElement(MORPH_ELLIPSE, Size(3, 3));
    morphologyEx(small, grad, MORPH_GRADIENT, morphKernel);

    Mat bw;
    threshold(grad, bw, 0.0, 255.0, THRESH_BINARY | THRESH_OTSU);

    Mat connected;
    morphKernel = getStructuringElement(MORPH_RECT, Size(9, 1));
    morphologyEx(bw, connected, MORPH_CLOSE, morphKernel);

    Mat mask = Mat::zeros(bw.size(), CV_8UC1);
    vector<vector<Point>> contours;
    vector<Vec4i> hierarchy;
    findContours(connected, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

    for(int idx = 0; idx >= 0; idx = hierarchy[idx][0]){
        Rect rect = boundingRect(contours[idx]);
        Mat maskROI(mask, rect);

        drawContours(mask, contours, idx, Scalar(255, 255, 255), CV_FILLED);

        double r = (double)countNonZero(maskROI) / (rect.width * rect.height);

        if (r > 0.45 && (rect.height > 8 && rect.width > 8)){
            rectangle(compressed, rect, Scalar(0, 255, 0), 2);
        }
    }

    imshow("text", compressed);
    waitKey(0);
}

int main(int argc, char* argv[]){
    detect_text(string("../image/test_marker.jpg"));
    return 0;
}