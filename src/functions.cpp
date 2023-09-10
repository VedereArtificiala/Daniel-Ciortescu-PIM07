/**
 * @file functions.cpp
 * @brief Acest fisier contine definitiile pentru functii.
 */

#include "header.hpp"

void process(imageAnalyze& data)
{
    imshow("Imaginea originala", data.grayScale);

    // Gaus pe original blureaza imaginea
    // Mat gaussMat = apelGaussianBlurCV(data.grayScale);
	Mat gaussMat = Filtru_Gauss(data.grayScale);
    imshow("Dupa aplicarea unui filtru Gauss", gaussMat);

    // Detectarea venelor adaptiveThresh
    Mat veinsImageOLD = detectareaVenelorOld(gaussMat, data);
	Mat veinsImage = detectareaVenelorAdaptiv(gaussMat);
    imshow("Vene Intermediare", veinsImage);

	Mat veneIntermediareMediane = apelMedianBlurCV(veinsImage);
    imshow("veneIntermediareMediane", veneIntermediareMediane);

	veinsImage = veneIntermediareMediane;


	{
		Mat veinsImageMedian = filtrumedian(veinsImage);

		// Afisarea imagnii in urma aplicarii unui filtru median
		imshow("Vene IntermediareMedian", veinsImageMedian);
	}

	Mat coloredImage = colorareaVenelor(data, veinsImage.data);
    // Afisarea imaginii cu venele evidentiate
    imshow("Vene Evidentiate", coloredImage);
    
    Mat coloredImageGaussinBlured = apelGaussianBlurCV(coloredImage);
	// Afisarea imaginii cu venele evidentiate
    imshow("Vene Evidentiate blured", coloredImageGaussinBlured);

	Mat coloredImageMedianBlured = apelMedianBlurCV(coloredImage);
	// Afisarea imaginii cu venele evidentiate
    imshow("Vene Evidentiate Median", coloredImageMedianBlured);
}

Mat detectareaVenelorAdaptiv(const Mat& inputImage) 
{
	unsigned char* img = inputImage.data;
	int w = inputImage.cols;
	int h = inputImage.rows;
	double maxVal = 255;
	int blockSize = 3;
	double C = 1.5;

	unsigned char* result = new unsigned char[w * h];

	cv::Mat inMat(h, w, CV_8UC1, img);
	cv::Mat binaryMat(h, w, CV_8UC1, result);

	cv::adaptiveThreshold(inMat, binaryMat, maxVal, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, blockSize, C); 


	Mat veinsImage(h, w, CV_8UC1, result);
	return veinsImage;
}

Mat detectareaVenelorOld(const Mat& inputImage, const imageAnalyze& data)
{
	unsigned char* img = inputImage.data;
	int w = inputImage.cols;
	int h = inputImage.rows;
	double th = 140;
	double maxVal = 255;
	
	unsigned char* result = new unsigned char[w * h];
    Mat inMat(h, w, CV_8UC1, img);
    Mat binaryMat(h, w, CV_8UC1, result);
    threshold(inMat, binaryMat, th, maxVal, cv::THRESH_BINARY);
    
	Mat veinsImage(h, w, CV_8UC1, result);
	unsigned char * veinsBinary = veinsImage.data;

	Mat coloredImage = data.color.clone();

	for (int i = 0; i < inputImage.rows; i++) 
    {
        for (int j = 0; j < inputImage.cols; j++) 
        {
            unsigned char pixel = veinsBinary[i * inputImage.cols + j];
            if (pixel == 255) 
            {
                coloredImage.at<Vec3b>(i, j) = Vec3b(0, 0, 255); // Rosu
            }
        }
    }

	return coloredImage;
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

Mat colorareaVenelor(const imageAnalyze& data, unsigned char* veinsImageRaw)
{
	// Colorarea venelor in imaginea color
    Mat coloredImage = data.color.clone();

    for (int i = 0; i < data.grayScale.rows; i++) 
    {
        for (int j = 0; j < data.grayScale.cols; j++) 
        {
            unsigned char center = veinsImageRaw[i * data.grayScale.cols + j];
            
            if (center == 0) 
            {
                int closeNeighbors = 0;
                
                for (int k = -1; k <= 1; k++) 
                {
                    for (int l = -1; l <= 1; l++) 
                    {
                        if (k == 0 && l == 0) 
                        {
                            continue;
                        }
                        int ii = i + k;
                        int jj = j + l;
                        
                        if (ii >= 0 && ii < data.grayScale.rows && jj >= 0 && jj < data.grayScale.cols) {

                            unsigned char neighb = veinsImageRaw[ii * data.grayScale.cols + jj];
                            
                            if (abs(neighb - center) <= 5) 
                            {
                                closeNeighbors++;
                            }
                        }
                    }
                }
                
                if (closeNeighbors >= 3) 
                {
                    coloredImage.at<Vec3b>(i, j) = Vec3b(145, 12, 255);
                }
            
            }
        }
    }

	return coloredImage;
}

Mat apelGaussianBlurCV(const Mat& inputImage)
{
 	cv::Size ksize(21, 21); // Dimensiunea kernel-ului Gaussian
    double sigma = 1;    	// Deviația standard
    cv::Mat outputImage;
    cv::GaussianBlur(inputImage, outputImage, ksize, sigma);

	return outputImage;
}

Mat apelMedianBlurCV(const Mat& inputImage)
{
	// Aplicați filtrul median
    int kernelSize = 3; 	// Dimensiunea kernel-ului median
    cv::Mat outputImage;
    cv::medianBlur(inputImage, outputImage, kernelSize);

	return outputImage;
}

Mat filtrumedian(const Mat& inputImage)
{
	unsigned char *img = inputImage.data;
	int w = inputImage.cols;
	int h = inputImage.rows;
	int fw = 20;
	int fh = 2;

	unsigned char *result = new unsigned char[w*h];
	unsigned char *fereastra = new unsigned char[fw*fh];
	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			int counter = -1;
			for (int i = 0; i < fw; i++)
			{
				for (int j = 0; j < fh; j++)
				{
					if ((y - fh / 2 + i >= 0) && (y - fh / 2 + i < h) && (x - fw / 2 + j) >= 0 && (x - fw / 2 + j < w))
					{
						fereastra[++counter] = img[(y - fh / 2 + i)*w + x - fw / 2 + j];
					}
				}
			}
			sort(fereastra, fereastra + fw * fh);
			result[y*w + x] = fereastra[counter / 2];
		}
	}
	delete[]fereastra;

	Mat veinsImageMedian(h, w, CV_8UC1, result);

	return veinsImageMedian;
}

double calcul_gauss(int s, int t, double sigma)
{
	double prod;
	int k = 2 * M_PI*pow(sigma, 2);
	double factor = (double)1. / k;
	prod = factor * exp(-(pow(s, 2) + pow(t, 2)) / (2 * pow(sigma, 2)));
	return prod;
}

Mat Filtru_Gauss(const Mat& inputImage)
{
	unsigned char* img = inputImage.data;
	int w =  inputImage.cols;
	int h = inputImage.rows;
	int dim_filtru_w = 3;
	int dim_filtru_h = 3;
	double sigma = 1.2;

    double filter5[3][3] = { 1, 2, 1,
                             2, 4, 2,
                             1, 2, 1 };
	double* temp5[3];
	for (int i = 0; i < 3; ++i)
    {
        temp5[i] = filter5[i];
    }
	double **filter_sol5 = temp5;

	unsigned char* result = new unsigned char[w*h];

	for (int i = 0; i < w * h; i++)
	{	
		result[i] = 0;
	}

	int a = (dim_filtru_w - 1) / 2;
	int b = (dim_filtru_h - 1) / 2;

	for (int y = a; y < h-a; y++)
	{
		for (int x = b; x < w-b; x++)
		{
			for (int s = -a; s <= a; s++)
			{
				for (int t = -b; t <= b; t++)
				{
					result[y*w + x] = result[y*w + x] + (calcul_gauss(s, t, sigma)* img[(y + s)*w + x + t]);
				}
			}
			if (result[y*w + x] > 255)
			{
				result[y*w + x] = 255;
			}
			if (result[y*w + x] < 0)
			{
				result[y*w + x] = 0;
			}
		}
	}

	Mat gaussMat(h, w, CV_8UC1, result);

	return gaussMat;
}
