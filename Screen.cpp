/*
    Screen.cpp
*/

#include <iostream>
#include <FreeImage.h>
#include "Screen.h"

Screen::Screen()
{
}

Screen::~Screen()
{
    delete[] pixelArray;
}

void Screen::init(const int width, const int height)
{
    this->width = width;
    this->height = height;
    pixelArray = new Byte[width*height*3];
}

void Screen::writePixel(const glm::dvec3& color, const int x, const int y)
{
    glm::dvec3 newColor = glm::clamp(color, 0.0, 1.0);
    for (int i = 0; i < 3; i++) {
        pixelArray[(width*y+x)*3 + 2 - i] = (Byte) (newColor[i]*255.0);
    }
}

void Screen::saveScreenshot(const char* outputFilename)
{
    FreeImage_Initialise();
    FIBITMAP *img = FreeImage_ConvertFromRawBits(pixelArray, width, height, width * 3, 
                                                 24, 0xFF0000, 0x00FF00, 0x0000FF, false);
    std::cout << "Saving screenshot: " << outputFilename << "\n";
    FreeImage_Save(FIF_PNG, img, outputFilename, 0);
    delete img;
    FreeImage_DeInitialise();
}

int Screen::getWidth()
{
    return width;
}

int Screen::getHeight()
{
    return height;
}
