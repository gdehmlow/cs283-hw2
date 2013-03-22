/*
    Triangle.h

    Defines intersections for triangle objects.
*/

#ifndef I_TRNGL
#define I_TRNGL

#include <glm/glm.hpp>
#include "Shape.h"
#include "Intersection.h"
#include "AABB.h"

class Triangle : public Shape {
    public:
        Triangle(glm::dvec3& vertex0, glm::dvec3& vertex1, glm::dvec3& vertex2);
        Triangle(glm::dvec3& vertex0, glm::dvec3& vertex1, glm::dvec3& vertex2,
                 glm::dvec3& normal0, glm::dvec3& normal1, glm::dvec3& normal2);
        ~Triangle();

        bool doesIntersect(Ray* ray, double tmax);
        int intersectionPoint(Ray* ray, Intersection* intersection);
        
        int getIntersectionPoint(const Ray& ray, Intersection& intersection, const double dt);
        bool doesRayIntersect(const Ray& ray, const double tmax, const double dt);
    private: 
        glm::dvec3* vertexArray;
        glm::dvec3* normalArray;
        glm::dvec3 faceNormal;
        bool isFace;
};

#endif
