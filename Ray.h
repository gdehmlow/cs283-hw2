/*
	Ray.h
*/

#ifndef I_RAY
#define I_RAY

#include <glm/glm.hpp>

typedef struct _Ray {
    glm::vec4 position;
    glm::vec4 direction;
} Ray;

#endif
