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

        glm::vec3 pathTraceRay(const Ray& ray, int depth, float weight, int bounce);

    private:
        int maxDepth;
        bool findClosestIntersection(const Ray& ray, int& closestIntersectionIndex, Intersection& closestIntersection);
        glm::vec3 directLighting(const Ray& ray, const Intersection& surfaceIntersection, const Primitive& intersectedObject);
        glm::vec3 indirectDiffuseLighting(const Intersection& surfaceIntersection, const Primitive& intersectedObject);
        Scene* scene;
};

#endif
