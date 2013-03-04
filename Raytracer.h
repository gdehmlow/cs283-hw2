/*
    Raytracer.h

    Static class. Traces a ray.
*/

#ifndef I_RAYTR
#define I_RAYTR

#include <glm/glm.hpp>
class Scene;
#include "Scene.h"

class Raytracer {
    public:
        static int maxDepth;
        static int traceRay(Scene* scene, Ray* ray, int depth, glm::vec3& color, float rayRIndex);
};

#endif
