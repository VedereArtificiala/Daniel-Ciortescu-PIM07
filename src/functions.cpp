/**
 * @file functions.cpp
 * @brief Acest fisier contine definitiile pentru functii.
 */

#include "header.hpp"

unsigned char* BinaryImage(unsigned char* img, int w, int h, double th, double maxVal)
{
    unsigned char* result = new unsigned char[w * h];
    Mat inMat(h, w, CV_8UC1, img);
    Mat binaryMat(h, w, CV_8UC1, result);
    threshold(inMat, binaryMat, th, maxVal, cv::THRESH_BINARY);
    return result;
}

void readImage(imageAnalyze& data, const char * path)
{
    data.color     = imread(path, IMREAD_COLOR);
    data.grayScale = imread(path, IMREAD_GRAYSCALE);
}

bool verify(const imageAnalyze& data)
{
    if (data.color.empty() || data.grayScale.empty()) 
    {
        std::cout << "Eroare: Nu s-a putut citi imaginea." << std::endl;
        return false;
    }
    else
    {
        return true;
    }
}

void process(const imageAnalyze& data)
{
    // Detectarea venelor
    unsigned char* veinsBinary = BinaryImage(data.grayScale.data, data.grayScale.cols, data.grayScale.rows, data.prag, data.valoareMaxima);

    // Crearea imaginii intermediare cu venele
    Mat veinsImage(data.grayScale.rows, data.grayScale.cols, CV_8UC1, veinsBinary);

    // Colorarea venelor in imaginea color
    Mat coloredImage = data.color.clone();

    for (int i = 0; i < data.grayScale.rows; i++) 
    {
        for (int j = 0; j < data.grayScale.cols; j++) 
        {
            unsigned char pixel = veinsBinary[i * data.grayScale.cols + j];
            if (pixel == 255) 
            {
                coloredImage.at<Vec3b>(i, j) = Vec3b(0, 0, 255); // Rosu
            }
        }
    }

    // Afisarea imaginii originale
    imshow("Original Image", data.grayScale);

    // Afisarea imaginii intermediare cu venele
    imshow("Vene Intermediare", veinsImage);
    
    // Afisarea imaginii cu venele evidentiate
    imshow("Vene Evidentiate", coloredImage);
}
