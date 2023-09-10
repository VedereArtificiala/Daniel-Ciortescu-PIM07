/**
 * @file header.hpp
 * @brief Acest fisier contine declaratii pentru functii si structuri utilizate in procesarea imaginilor.
 */

#ifndef _HEADER_H_
#define _HEADER_H_

#include <opencv2/opencv.hpp>
#include <iostream>
#include <stack>

using namespace cv;
using namespace std;

/**
 * @brief Structura pentru a stoca datele imaginii si parametrii de analiza.
 */
struct imageAnalyze
{
    Mat color;                  ///< Imaginea colora.
    Mat grayScale;              ///< Imaginea alb-negru (grayscale).
};


Mat detectareaVenelorAdaptiv(const Mat& inputImage);
Mat detectareaVenelorOld(const Mat& inputImage, imageAnalyze const&);

/**
 * @brief Functie pentru a citi o imagine si a stoca datele in structura data.
 *
 * @param data Structura in care se vor stoca datele imaginii.
 * @param path Calea catre imaginea de citit.
 */
void readImage(imageAnalyze& data, const char * path);

/**
 * @brief Functie pentru a verifica daca datele imaginii sunt valide.
 *
 * @param data Structura care contine datele imaginii.
 * @return true daca datele sunt valide, false in caz contrar.
 */
bool verify(const imageAnalyze& data);

/**
 * @brief Functie pentru a procesa imaginea in functie de datele stocate in structura data.
 *
 * @param data Structura care contine datele imaginii.
 */
void process(imageAnalyze& data);

double calcul_gauss(int s, int t, double sigma);

Mat Filtru_Gauss(const Mat& inputImage);

Mat colorareaVenelor(const imageAnalyze& data, unsigned char* veinsImageRaw);

Mat filtrumedian(const Mat& inputImage);

Mat apelGaussianBlurCV(const Mat& inputImage);

Mat apelMedianBlurCV(const Mat& inputImage);

#endif // _HEADER_H_
