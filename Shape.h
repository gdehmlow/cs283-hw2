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
        virtual int getIntersectionPoint(const Ray& ray, Intersection& intersection, const double dt) = 0;
        virtual bool doesRayIntersect(const Ray& ray, const double tmax, const double dt) = 0;

        virtual bool doesIntersect(Ray* ray, double tmax) = 0;
        virtual int intersectionPoint(Ray* ray, Intersection* intersection) = 0;
};

#endif
