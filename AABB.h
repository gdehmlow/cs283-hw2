/*
    AABB.h

    Public bounding box for primitives
*/

#ifndef I_AABB
#define I_AABB

#include <glm/glm.hpp>

class AABB {
    public:
        AABB();
        bool doesIntersectAABB(AABB* otherBox);
        bool isWithin(glm::vec3& point);

        glm::vec3 minimum;
        glm::vec3 maximum;
};

#endif