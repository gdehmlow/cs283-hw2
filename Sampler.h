/*
    Sampler.h

    Sampler generates a float pair (x,y) to be used to generate rays from the camera.
    There is the option for it to be random or not. Will generate multiple samples per pixel
    (necessary for Monte Carlo path tracing) if samplesPerPixel > 1.

    Restriction: if samplerType is not random, samplesPerPixel must be a square number (sPP=x^2 
    for some x), or it will be rounded to nearest square number.
*/

#ifndef I_SAMPLER
#define I_SAMPLER

typedef struct _Sample {
    float x;
    float y;
} Sample;

enum SamplerType { REGULAR, JITTERED, RANDOM };

class Sampler {
    public:
        bool hasSamples();
        Sample getSample();
        Sampler(int x, int y, int samplesPerPixel, SamplerType t);
        ~Sampler();

    private:
        bool hasSamps;
        int x;
        int y;
        int currentSubPixelCount;
        int sqrtSamplesPerPixel;
        int samplesPerPixel;
        float subPixelOffset;
        SamplerType samplerType;
};

#endif
