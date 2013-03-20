/*
    Sampler.cpp
*/

#include <cmath>
#include <ctime>
#include <stdlib.h>
#include <iostream>
#include "Sampler.h"
using namespace std;

Sampler::Sampler(int x, int y, int samplesPerPixel, SamplerType t)
{
    this->samplerType           = t;
    this->x                     = x;
    this->y                     = y;
    this->hasSamps              = true;
    this->currentSubPixelCount  = 0;
    this->sqrtSamplesPerPixel   = static_cast<int>(sqrt(static_cast<float>(samplesPerPixel)) + .5f);
    this->samplesPerPixel       = samplesPerPixel;
    this->subPixelOffset        = 0.5/((float) this->sqrtSamplesPerPixel);
}

Sampler::~Sampler()
{
}

bool Sampler::hasSamples()
{
    return hasSamps;
}

Sample Sampler::getSample()
{
    Sample sample;

    if (samplerType == REGULAR) {
        float xInSubPixel = (float)(currentSubPixelCount % sqrtSamplesPerPixel) / sqrtSamplesPerPixel;
        float yInSubPixel = (float)(currentSubPixelCount / sqrtSamplesPerPixel) / sqrtSamplesPerPixel;
        sample.x = (float)x + xInSubPixel + subPixelOffset;
        sample.y = (float)y + yInSubPixel + subPixelOffset;
    } else if (samplerType == JITTERED) {
        float xInSubPixel = (float)(currentSubPixelCount % sqrtSamplesPerPixel) / sqrtSamplesPerPixel;
        float yInSubPixel = (float)(currentSubPixelCount / sqrtSamplesPerPixel) / sqrtSamplesPerPixel;
        sample.x = (float)x + xInSubPixel + subPixelOffset * 2.0 * rand() / double(RAND_MAX);
        sample.y = (float)y + yInSubPixel + subPixelOffset * 2.0 * rand() / double(RAND_MAX);
    } else if (samplerType == RANDOM) {
        sample.x = (float)x + rand() / double(RAND_MAX);
        sample.y = (float)y + rand() / double(RAND_MAX);
    }

    currentSubPixelCount++;

    if (currentSubPixelCount > samplesPerPixel - 2) {
        currentSubPixelCount = 0;
        hasSamps = false;
    }

    return sample;
}
