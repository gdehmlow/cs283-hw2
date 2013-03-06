/*
    Sampler.cpp
*/

#include "Sampler.h"

Sampler::Sampler(int width, int height, int samplesPerPixel, bool isRandom)
{
    this->screenWidth       = (float) width;  // Making them floats now so I don't have to later when
    this->screenHeight      = (float) height; // doing the sampling math.
    this->samplesPerPixel   = samplesPerPixel;
    this->isRandom          = isRandom;
    this->currentX          = 0;
    this->currentY          = 0;
}

Sampler::~Sampler()
{
}

void Sampler::getSample(Sample* sample)
{
    
}