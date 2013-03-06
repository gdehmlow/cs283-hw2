/*
    Sampler.h

    Sampler generates a float pair (x,y) to be used to generate rays from the camera.
    There is the option for it to be random or not. Will generate multiple samples per pixel
    (necessary for Monte Carlo path tracing) if samplesPerPixel > 1.

    Restriction: if not random, samplesPerPixel must be a square number (sPP=x^2 for some x). 
*/

#ifndef I_SAMPLER
#define I_SAMPLER

typedef struct _Sample {
    float x;
    float y;
} Sample;

enum SampleType { REGULAR, JITTERED, RANDOM };

class Sampler {
    public:
        void getSample(Sample* sample);
        Sampler(int width, int height, int samplesPerPixel, SampleType t);
        ~Sampler();

    private:
        float screenWidth; 
        float screenHeight;
        int currentX;
        int currentY;
        int samplesPerPixel;
        bool isRandom;
};

#endif
