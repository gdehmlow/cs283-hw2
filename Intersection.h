/*
    Intersection.h
    
    Struct to hold intersection values;
*/

#ifndef I_INTRS
#define I_INTRS

#include <glm/glm.hpp>

typedef struct _Intersection {
    glm::vec4 position;
    glm::vec4 normal;
    float t;
} Intersection;

#endif
