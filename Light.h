/*
    Light.h
*/

#ifndef I_LIGHT
#define I_LIGHT

#include <glm/glm.hpp>

typedef struct _Light {
    glm::vec4 posdir;
    glm::vec3 color;
    bool point;
} Light;

#endif