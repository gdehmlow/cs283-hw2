/*
    Shape.h

    Abstract class which provides intersection methods.
*/

#ifndef I_SHAPE
#define I_SHAPE

#include <glm/glm.hpp>
#include "Ray.h"
#include "Intersection.h"
#include "AABB.h"

class Shape {
    public:
        virtual int intersectionPoint(Ray* ray, Intersection* intersection) = 0;
        virtual int getIntersectionPoint(const Ray& ray, Intersection& intersection) = 0;
        virtual bool doesRayIntersect(const Ray& ray, const float tmax) = 0;
        virtual bool doesIntersect(Ray* ray, float tmax) = 0;
        virtual void createAABB(glm::mat3& transformation) = 0;
        virtual AABB* getAABB() = 0;
};

#endif
