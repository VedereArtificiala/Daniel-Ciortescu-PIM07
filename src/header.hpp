/**
 * @file header.hpp
 * @brief Acest fisier contine declaratii pentru functii si structuri utilizate in procesarea imaginilor.
 */

#ifndef _HEADER_H_
#define _HEADER_H_

#include <opencv2/opencv.hpp>
#include <iostream>

using namespace cv;

/**
 * @brief Structura pentru a stoca datele imaginii si parametrii de analiza.
 */
struct imageAnalyze
{
    Mat color;                  ///< Imaginea colora.
    Mat grayScale;              ///< Imaginea alb-negru (grayscale).

    double prag = 150;          ///< Pragul pentru procesul de binarizare (default: 150).
    double valoareMaxima = 255; ///< Valoarea maxima pentru procesul de binarizare (default: 255).
};

/**
 * @brief Functie pentru a efectua binarizarea imaginii.
 *
 * @param img Pointer catre imaginea sursa.
 * @param w Latimea imaginii.
 * @param h inaltimea imaginii.
 * @param th Pragul pentru binarizare.
 * @param maxVal Valoarea maxima pentru binarizare.
 * @return Pointer catre imaginea binarizata.
 */
unsigned char* BinaryImage(unsigned char* img, int w, int h, double th, double maxVal);

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
void process(const imageAnalyze& data);

#endif // _HEADER_H_
