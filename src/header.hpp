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

/**
 * @brief Calculeaza valorile kernelului pentru filtrul Gaussian.
 *
 * @param s Variabila spatiala.
 * @param t Variabila temporala.
 * @param sigma Valoarea deviatiei standard.
 * @return Valoarea calculata pentru kernelul Gaussian.
 */
double calcul_gauss(int s, int t, double sigma);

/**
 * @brief Aplica filtrul Gaussian pe o imagine folosind metoda personalizata.
 *
 * @param inputImage Imaginea de intrare.
 * @return Imaginea rezultata dupa aplicarea filtrului Gaussian.
 */
Mat Filtru_Gauss(const Mat& inputImage);

/**
 * @brief Creeaza o imagine colorata pentru vene pe baza datelor și imaginii brute.
 *
 * @param data Structura care contine datele imaginii.
 * @param veinsImageRaw Imaginea bruta pentru vene.
 * @return Imaginea colorata pentru vene.
 */
Mat colorareaVenelor(const imageAnalyze& data, unsigned char* veinsImageRaw);

/**
 * @brief Aplica filtrul median pe o imagine folosind metoda personalizata.
 *
 * @param inputImage Imaginea de intrare.
 * @return Imaginea rezultata dupa aplicarea filtrului median.
 */
Mat filtrumedian(const Mat& inputImage);

/**
 * @brief Aplica filtrul Gaussian pe o imagine folosind OpenCV.
 *
 * @param inputImage Imaginea de intrare.
 * @return Imaginea rezultata dupa aplicarea filtrului Gaussian folosind OpenCV.
 */
Mat apelGaussianBlurCV(const Mat& inputImage);

/**
 * @brief Aplica filtrul median pe o imagine folosind OpenCV.
 *
 * @param inputImage Imaginea de intrare.
 * @return Imaginea rezultata dupa aplicarea filtrului median folosind OpenCV.
 */
Mat apelMedianBlurCV(const Mat& inputImage);

/**
 * @brief Aplica operatia de dilatare pe o imagine folosind OpenCV.
 *
 * @param inputImage Imaginea de intrare.
 * @return Imaginea rezultata dupa aplicarea dilatarii folosind OpenCV.
 */
Mat apelDilatareCV(const Mat& inputImage);

/**
 * @brief Aplica operatia de dilatare pe o imagine folosind metoda personalizata.
 *
 * @param inputImage Imaginea de intrare.
 * @return Imaginea rezultata dupa aplicarea dilatarii folosind metoda personalizata.
 */
Mat filtruDilatare(const Mat& inputImage);

/**
 * @brief Aplica operatia de eroziune pe o imagine folosind OpenCV.
 *
 * @param inputImage Imaginea de intrare.
 * @return Imaginea rezultata dupa aplicarea eroziunii folosind OpenCV.
 */
Mat apelEroziuneCV(const Mat& inputImage);

/**
 * @brief Aplica operatia de eroziune pe o imagine folosind metoda personalizata.
 *
 * @param inputImage Imaginea de intrare.
 * @return Imaginea rezultata dupa aplicarea eroziunii folosind metoda personalizată.
 */
Mat filtruEroziune(const Mat& inputImage);

#endif // _HEADER_H_
