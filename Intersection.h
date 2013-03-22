/*
    Intersection.h
    
    Struct to hold intersection values;
*/

#ifndef I_INTRS
#define I_INTRS

#include <glm/glm.hpp>

typedef struct _Intersection {
    glm::dvec4 position;
    glm::dvec4 normal;
    float t;
} Intersection;

#endif
