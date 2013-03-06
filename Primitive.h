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
        int intersectionPoint(Ray* ray, Intersection* intersection);
        bool doesIntersect(Ray* ray, float tmax);
        void setAmbient(glm::vec3& ambient);
        void setMaterial(Material* material);
        void setTransformation(glm::mat4& transformation);
        AABB* getAABB();

        Shape* shape;
        Material* material;
        glm::vec3 ambient;
        glm::mat4 transformation;
        glm::mat4 inverseTransformation;
        glm::mat4 inverseTranspose;
};

#endif
