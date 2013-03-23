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
        Sphere(double x, double y, double z, double radius);
        ~Sphere();

        void setEndPosition(double x, double y, double z);
        int getIntersectionPoint(const Ray& ray, Intersection& intersection, const double dt);
        bool doesRayIntersect(const Ray& ray, const double tmax, const double dt);

    private:
        glm::dvec3 posit;
        glm::dvec3 endposit;
        double radius;
        bool isMoving;
};

#endif
