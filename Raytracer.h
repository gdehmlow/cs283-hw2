/*
    Raytracer.h

    Static class. Traces a ray.
*/

#ifndef I_RAYTR
#define I_RAYTR

#include <glm/glm.hpp>
class Scene;
#include "Scene.h"
#include "Intersection.h"
#include "Primitive.h"

class Raytracer {
    public:
        Raytracer(Scene* scene);
        ~Raytracer();

        int traceRay(Ray* ray, int depth, glm::vec3& color, float rayRIndex);

    private:
        int maxDepth;
        glm::vec3 directLighting(Ray* ray, Intersection* lightIntersect, Primitive* intersectedObject);
        Scene* scene;
};

#endif
