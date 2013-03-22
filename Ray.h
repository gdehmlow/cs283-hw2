/*
	Ray.h
*/

#ifndef I_RAY
#define I_RAY

#include <glm/glm.hpp>

typedef struct _Ray {
    glm::dvec4 position;
    glm::dvec4 direction;
    double t;
} Ray;

#endif
