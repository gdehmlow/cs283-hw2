/*
    Primitive.h
    
    Public; holds shape and material data for use by the raytracer and shader.
*/

#ifndef I_PRIM
#define I_PRIM

#include <glm/glm.hpp>

#include "Shape.h"
#include "Material.h"
#include "Ray.h"
#include "Intersection.h"
#include "AABB.h"

class Primitive {
    public:
        Primitive(Shape* shape);
        ~Primitive();

        int getIntersectionPoint(const Ray& ray, Intersection& intersection, const double dt);
        bool doesRayIntersect(const Ray& ray, const float tmax, const double dt);

        // Delete soon
        int intersectionPoint(Ray* ray, Intersection* intersection);
        bool doesIntersect(Ray* ray, float tmax);

        void setAmbient(const glm::dvec3& ambient);
        void setMaterial(const Material& material);
        void setTransformation(glm::dmat4& transformation);

        Shape* shape;
        Material* material;
        glm::dvec3 ambient;
        glm::dmat4 transformation;
        glm::dmat4 inverseTransformation;
        glm::dmat4 inverseTranspose;
};

#endif
