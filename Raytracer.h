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

        glm::dvec3 pathTraceRay(const Ray& ray, int depth, double weight, int bounce);

    private:
        int maxDepth;
        bool findClosestIntersection(const Ray& ray, int& closestIntersectionIndex, Intersection& closestIntersection);
        glm::dvec3 directLighting(const Ray& ray, const Intersection& surfaceIntersection, 
                                  const Primitive& intersectedObject);
        glm::dvec3 indirectDiffuseLighting(const Intersection& surfaceIntersection, const Primitive& intersectedObject);
        void sampleUniformHemisphere(Ray& ray, const Intersection& surfaceIntersection); 
        void sampleCosineWeightedHemisphere(Ray& ray, const Intersection& surfaceIntersection);
        void sampleSpecularLobe(Ray& ray, const Intersection& surfaceIntersection, const glm::dvec3& reflection, 
                                const double shininess);
        Scene* scene;
};

#endif
