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

        void traceRay(Ray* ray, int depth, glm::vec3& color, float weight, int bounce, float rayRIndex);
        glm::vec3 pathTraceRay(Ray& ray, int depth, glm::vec3& color, float weight, int bounce);

    private:
        int maxDepth;
        bool findClosestIntersection(Ray& ray, int& closestIntersectionIndex, Intersection& closestIntersection);
        glm::vec3 directLighting(Ray* ray, Intersection* lightIntersect, Primitive* intersectedObject);
        Scene* scene;
};

#endif
