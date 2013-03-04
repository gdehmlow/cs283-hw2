/*
    Sphere.h

    Defines intersections for sphere objects
*/

#ifndef I_SPHERE
#define I_SPHERE

#include <glm/glm.hpp>
#include "Shape.h"
#include "Intersection.h"
#include "AABB.h"

class Sphere : public Shape {
    public:
        Sphere(float x, float y, float z, float radius);
        bool doesIntersect(Ray* ray, float tmax);
        int intersectionPoint(Ray* ray, Intersection* intersection);
        void createAABB(glm::mat3& transformation);
        AABB* getAABB();
    private:
        AABB* aabb;
        glm::vec3 posit;
        float radius;
};

#endif
