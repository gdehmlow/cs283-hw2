/*
    Sphere.cpp
*/

#include <iostream>
#include <math.h>
#include "Sphere.h"

typedef glm::dmat3 mat3;
typedef glm::dmat4 mat4; 
typedef glm::dvec3 vec3; 
typedef glm::dvec4 vec4; 

Sphere::Sphere(double x, double y, double z, double radius)
{
    this->posit = vec3(x,y,z);
    this->radius = radius;
}

Sphere::~Sphere()
{
}

void Sphere::setEndPosition(double x, double y, double z)
{
    this->endposit = vec3(x,y,z);
}

bool Sphere::doesIntersect(Ray* ray, double tmax)
{
    vec3 dir = vec3(ray->direction[0],ray->direction[1],ray->direction[2]);
    vec3 pos = vec3(ray->position[0],ray->position[1],ray->position[2]);
    double A = glm::dot(dir,dir);
    double B = 2 * glm::dot(dir, (pos - this->posit));
    double C = glm::dot(pos - this->posit, pos - this->posit) - radius*radius;
    double disc = B*B - 4*A*C;
    double t;
    
    if (disc < 0) {
        return false;
    } else {
        double q0 = (-B + sqrt(disc)) / (2*A);
        double q1 = (-B - sqrt(disc)) / (2*A);

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

int Sphere::getIntersectionPoint(const Ray& ray, Intersection& intersection, const double dt)
{
    vec3 currentPosition;
    if (isMoving) {
        currentPosition = this->posit*(1.0 - dt) + this->endposit*(dt);
    } else {
        currentPosition = this->posit;
    }
    
    vec3 dir = vec3(ray.direction[0],ray.direction[1],ray.direction[2]);
    vec3 pos = vec3(ray.position[0],ray.position[1],ray.position[2]);
    double A = glm::dot(dir,dir);
    double B = 2 * glm::dot(dir, (pos - currentPosition));
    double C = glm::dot(pos - currentPosition, pos - currentPosition) - radius*radius;
    double disc = B*B - 4*A*C;
    int ret = 1;

    if (disc < 0) {
        return 0;    
    } else {
        double q0 = (-B + sqrt(disc)) / (2*A);
        double q1 = (-B - sqrt(disc)) / (2*A);
        double t;
        if (q0 > 0 && q1 > 0) {
            t = q0 >= q1 ? q1 : q0; // pick the minimum
        } else if (q0 == q1) {
            return 0;
        } else {
            t = q0 > 0 ? q0 : q1; // pick the positive one
            ret = -1; // inside sphere!
        }
        vec3 position = pos + dir*t;
        vec3 normal = glm::normalize(position - currentPosition);

        intersection.position = vec4(position,1);
        intersection.normal = vec4(normal,0);
        intersection.t = t;
    }
    return ret;
}

bool Sphere::doesRayIntersect(const Ray& ray, const double tmax, const double dt)
{
    vec3 currentPosition;
    if (isMoving) {
        currentPosition = this->posit*(1.0 - dt) + this->endposit*(dt);
    } else {
        currentPosition = this->posit;
    }

    vec3 dir = vec3(ray.direction[0],ray.direction[1],ray.direction[2]);
    vec3 pos = vec3(ray.position[0],ray.position[1],ray.position[2]);
    double A = glm::dot(dir,dir);
    double B = 2 * glm::dot(dir, (pos - currentPosition));
    double C = glm::dot(pos - currentPosition, pos - currentPosition) - radius*radius;
    double disc = B*B - 4*A*C;
    double t;
    
    if (disc < 0) {
        return false;
    } else {
        double q0 = (-B + sqrt(disc)) / (2*A);
        double q1 = (-B - sqrt(disc)) / (2*A);

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
    double A = glm::dot(dir,dir);
    double B = 2 * glm::dot(dir, (pos - this->posit));
    double C = glm::dot(pos - this->posit, pos - this->posit) - radius*radius;
    double disc = B*B - 4*A*C;
    int ret = 1;

    if (disc < 0) {
        return 0;    
    } else {
        double q0 = (-B + sqrt(disc)) / (2*A);
        double q1 = (-B - sqrt(disc)) / (2*A);
        double t;
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
