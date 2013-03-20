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
        Triangle(glm::vec3& vertex0, glm::vec3& vertex1, glm::vec3& vertex2);
        Triangle(glm::vec3& vertex0, glm::vec3& vertex1, glm::vec3& vertex2,
                 glm::vec3& normal0, glm::vec3& normal1, glm::vec3& normal2);
        ~Triangle();

        bool doesIntersect(Ray* ray, float tmax);
        int intersectionPoint(Ray* ray, Intersection* intersection);
        
        int getIntersectionPoint(const Ray& ray, Intersection& intersection);
        bool doesRayIntersect(const Ray& ray, const float tmax);
        void createAABB(glm::mat3& transformation);
        AABB* getAABB();
    private: 
        glm::vec3* vertexArray;
        glm::vec3* normalArray;
        glm::vec3 faceNormal;
        bool isFace;
        AABB* aabb;
};

#endif
