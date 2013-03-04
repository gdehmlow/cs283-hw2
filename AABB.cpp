/*
    AABB.cpp
*/

#include "AABB.h"

typedef glm::vec3 vec3;

AABB::AABB()
{
    this->minimum = vec3(0.0,0.0,0.0);
    this->maximum = vec3(0.0,0.0,0.0);
}

bool AABB::doesIntersectAABB(AABB* otherBox)
{
    return (this->maximum.x > otherBox->minimum.x) &&
           (this->maximum.y > otherBox->minimum.y) &&
           (this->maximum.z > otherBox->minimum.z) &&
           (this->minimum.x < otherBox->maximum.x) &&
           (this->minimum.y < otherBox->maximum.y) &&
           (this->minimum.z < otherBox->maximum.z);
}

bool AABB::isWithin(glm::vec3& point)
{
    return (point.x >= this->minimum.x) &&
           (point.y >= this->minimum.y) &&
           (point.z >= this->minimum.z) &&
           (point.x <= this->maximum.x) &&
           (point.y <= this->maximum.y) &&
           (point.z <= this->maximum.z);
}
