/*
    Screen.h

    Responsible for holding pixel color data and displaying the image.
*/

#ifndef I_SCREEN
#define I_SCREEN

#include <glm/glm.hpp>

typedef unsigned char Byte;

class Screen {
    public:
        Screen();
        ~Screen();
        void init(const int width, const int height);
        void writePixel(const glm::dvec3& color, const int x, const int y);
        void saveScreenshot(const char* outputFilename);

        int getWidth();
        int getHeight();
    private:
        int width;
        int height;
        Byte* pixelArray;
};

#endif
