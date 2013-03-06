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
        void init(int width, int height);
        void writePixel(glm::vec3& color, int x, int y);
        void saveScreenshot(const char* outputFilename);

        int getWidth();
        int getHeight();
    private:
        int width;
        int height;
        Byte* pixelArray;
};

#endif
