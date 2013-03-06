/*
    Sphere.cpp
*/

#include <iostream>
#include <math.h>
#include "Sphere.h"

typedef glm::mat3 mat3;
typedef glm::mat4 mat4; 
typedef glm::vec3 vec3; 
typedef glm::vec4 vec4; 

Sphere::Sphere(float x, float y, float z, float radius)
{
    this->posit = vec3(x,y,z);
    this->radius = radius;
}

Sphere::~Sphere()
{
    //delete aabb;
}

bool Sphere::doesIntersect(Ray* ray, float tmax)
{
    vec3 dir = vec3(ray->direction[0],ray->direction[1],ray->direction[2]);
    vec3 pos = vec3(ray->position[0],ray->position[1],ray->position[2]);
    float A = glm::dot(dir,dir);
    float B = 2 * glm::dot(dir, (pos - this->posit));
    float C = glm::dot(pos - this->posit, pos - this->posit) - radius*radius;
    float disc = B*B - 4*A*C;
    float t;
    
    if (disc < 0) {
        return false;
    } else {
        float q0 = (-B + sqrt(disc)) / (2*A);
        float q1 = (-B - sqrt(disc)) / (2*A);

        if (q0 > 0 && q1 > 0) {
            t = q0 >= q1 ? q1 : q0; // pick the minimum
        } else if (q0 == q1) {
            return false;
        } else {
            t = q0 > 0 ? q0 : q1; // pick the positive one
        }

        if (t >= tmax || t < 0.0) {
            return false;
        }
    }
    return true;
}

int Sphere::intersectionPoint(Ray* ray, Intersection* intersection)
{
    vec3 dir = vec3(ray->direction[0],ray->direction[1],ray->direction[2]);
    vec3 pos = vec3(ray->position[0],ray->position[1],ray->position[2]);
    float A = glm::dot(dir,dir);
    float B = 2 * glm::dot(dir, (pos - this->posit));
    float C = glm::dot(pos - this->posit, pos - this->posit) - radius*radius;
    float disc = B*B - 4*A*C;
    int ret = 1;

    if (disc < 0) {
        return 0;    
    } else {
        float q0 = (-B + sqrt(disc)) / (2*A);
        float q1 = (-B - sqrt(disc)) / (2*A);
        float t;
        if (q0 > 0 && q1 > 0) {
            t = q0 >= q1 ? q1 : q0; // pick the minimum
        } else if (q0 == q1) {
            return 0;
        } else {
            t = q0 > 0 ? q0 : q1; // pick the positive one
            ret = -1; // inside sphere!
        }
        vec3 position = pos + dir*t;
        vec3 normal = glm::normalize(position - this->posit);

        intersection->position = vec4(position,1);
        intersection->normal = vec4(normal,0);
        intersection->t = t;
    }
    return ret;
}

void Sphere::createAABB(glm::mat3& transformation)
{
    this->aabb = new AABB();
    this->aabb->minimum = (this->posit - this->radius)*transformation;
    this->aabb->maximum = (this->posit + this->radius)*transformation;
}

AABB* Sphere::getAABB() {
    return aabb;
}
