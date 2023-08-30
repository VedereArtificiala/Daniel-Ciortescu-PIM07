/**
 * @file header.hpp
 * @brief Acest fisier apelul functiilor pentru procesarea imaginilor.
 */


#include "header.hpp"

/**
 * @brief Functia principala a programului.
 * 
 * Aceasta functie este punctul de intrare in program.
 */
int main() 
{
    // Initializam o structura imageAnalyze pentru stocarea datelor imaginii.
    imageAnalyze data;

    // Citim imaginea de la calea specificata.
    readImage(data, "/home/dani/Documents/imgproj/Daniel-Ciortescu-PIM07/img/03.jpg");
    
    // Verificam daca citirea imaginii a fost reusita.
    if (verify(data))
    {
        // Procesam imaginea pentru a evidentia venele si afisam rezultatul.
        process(data);
    }

    // Asteptam pana cand utilizatorul apasa o tasta.
    waitKey(0);

    // Programul se incheie cu succes.
    return 0;
}
