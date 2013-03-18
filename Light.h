/*
    Light.h
*/

#ifndef I_LIGHT
#define I_LIGHT

#include <glm/glm.hpp>

enum LightType { POINT, DIRECTIONAL, AREA };
enum AreaType { RECT };

typedef struct _Light {
    glm::vec4 posdir;
    glm::vec3 color;
    LightType type;

    // For rectangular area light sources
    //
    // ^ = = = = = = o 
    // | <-upVec     |
    // |             |
    // *------------>o
    // ^        ^
    // posDir   rightVec

    AreaType areaType;
    int numSamples; // Per square side! So numSamples = 3 -> 9 samples taken. 
    glm::vec4 upStep;
    glm::vec4 rightStep;
} Light;

#endif