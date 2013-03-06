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

void Screen::init(int width, int height)
{
    this->width = width;
    this->height = height;
    pixelArray = new Byte[width*height*3];
}

void Screen::writePixel(glm::vec3& color, int x, int y)
{
    // Note: 2 - i for color location because FreeImage expects BGR despite mine giving it 
    // different bitmask data.
    for (int i = 0; i < 3; i++) {
        pixelArray[(width*y+x)*3 + 2 - i] = (Byte) (color[i]*255.0);
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
