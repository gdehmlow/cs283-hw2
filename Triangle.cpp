/*
    Triangle.cpp
*/

#include <iostream>
#include "Triangle.h"

typedef glm::dmat3 mat3;
typedef glm::dmat4 mat4; 
typedef glm::dvec3 vec3; 
typedef glm::dvec4 vec4; 

Triangle::Triangle(vec3& vertex0, vec3& vertex1, vec3& vertex2)
{
    vertexArray = new vec3[3];
    vertexArray[0] = vertex0;
    vertexArray[1] = vertex1;
    vertexArray[2] = vertex2;
    vec3 norm = glm::normalize(glm::cross(vertex2 - vertex1, vertex0 - vertex1));
    normalArray = new vec3[3];
    normalArray[0] = norm;
    normalArray[1] = norm;
    normalArray[2] = norm;
    faceNormal = norm;
    isFace = true;
}

Triangle::Triangle(vec3& vertex0, vec3& vertex1, vec3& vertex2,
                   vec3& normal0, vec3& normal1, vec3& normal2)
{
    vertexArray = new vec3[3];
    vertexArray[0] = vertex0;
    vertexArray[1] = vertex1;
    vertexArray[2] = vertex2;
    normalArray = new vec3[3];
    normalArray[0] = normal0;
    normalArray[1] = normal1;
    normalArray[2] = normal2;
    faceNormal = glm::normalize(glm::cross(vertex2 - vertex1, vertex0 - vertex1));
    isFace = false;
}

Triangle::~Triangle()
{
    //delete aabb;
    delete[] vertexArray;
    delete[] normalArray;
}

bool Triangle::doesIntersect(Ray* ray, double tmax)
{
    vec3 dir = vec3(ray->direction[0],ray->direction[1],ray->direction[2]);
    vec3 pos = vec3(ray->position[0],ray->position[1],ray->position[2]);

    double denom = glm::dot(dir,faceNormal);

    if (denom == 0.0) {
        return false;
    } else {
        double t = (glm::dot(vertexArray[1],faceNormal) - glm::dot(pos,faceNormal)) / denom;
        vec3 p = pos + dir * t;

        vec3 u = p - vertexArray[1];             
        vec3 v = vertexArray[0] - vertexArray[1];
        vec3 w = vertexArray[2] - vertexArray[1];

        double ww = glm::dot(w,w);
        double wv = glm::dot(w,v);
        double uv = glm::dot(u,v);
        double uw = glm::dot(u,w);
        double vv = glm::dot(v,v);

        double beta = (ww*uv - wv*uw) / (vv*ww - wv*wv);
        double gamma = (vv*uw - wv*uv) / (vv*ww - wv*wv);

        if (beta >= 0 && gamma >= 0 && (beta + gamma <= 1) && t >= 0.0 && t <= tmax) {
            return true;
        } else {
            return false;
        }
    }
    return false;
}

int Triangle::getIntersectionPoint(const Ray& ray, Intersection& intersection, const double dt)
{
    vec3 dir = vec3(ray.direction[0],ray.direction[1],ray.direction[2]);
    vec3 pos = vec3(ray.position[0],ray.position[1],ray.position[2]);
    vec3 normal;
    double denom = glm::dot(dir,faceNormal);

    if (denom == 0.0) {
        return 0;
    } else {
        double t = (glm::dot(vertexArray[1],faceNormal) - glm::dot(pos,faceNormal)) / denom;

        intersection.t = t;
        vec3 p = pos + dir * t;
        intersection.position = vec4(p,1);

        // Barycentric coordinate calculations
        vec3 u = p - vertexArray[1];             
        vec3 v = vertexArray[0] - vertexArray[1];
        vec3 w = vertexArray[2] - vertexArray[1];

        double ww = glm::dot(w,w);
        double wv = glm::dot(w,v);
        double uv = glm::dot(u,v);
        double uw = glm::dot(u,w);
        double vv = glm::dot(v,v);
        double beta = (ww*uv - wv*uw) / (vv*ww - wv*wv);
        double gamma = (vv*uw - wv*uv) / (vv*ww - wv*wv);

        if (isFace) {
            normal = faceNormal;
        } else {
            normal = normalArray[1] + beta * (normalArray[0] - normalArray[1]) + 
                     gamma * (normalArray[2] - normalArray[1]);
        }
        intersection.normal = vec4(normal,0);

        if (beta >= 0 && gamma >= 0 && (beta + gamma <= 1)) {
            return 1;
        } else {
            return 0;
        }
    }
}

bool Triangle::doesRayIntersect(const Ray& ray, const double tmax, const double dt)
{
    vec3 dir = vec3(ray.direction[0],ray.direction[1],ray.direction[2]);
    vec3 pos = vec3(ray.position[0],ray.position[1],ray.position[2]);

    double denom = glm::dot(dir,faceNormal);

    if (denom == 0.0) {
        return false;
    } else {
        double t = (glm::dot(vertexArray[1],faceNormal) - glm::dot(pos,faceNormal)) / denom;
        vec3 p = pos + dir * t;

        vec3 u = p - vertexArray[1];             
        vec3 v = vertexArray[0] - vertexArray[1];
        vec3 w = vertexArray[2] - vertexArray[1];

        double ww = glm::dot(w,w);
        double wv = glm::dot(w,v);
        double uv = glm::dot(u,v);
        double uw = glm::dot(u,w);
        double vv = glm::dot(v,v);

        double beta = (ww*uv - wv*uw) / (vv*ww - wv*wv);
        double gamma = (vv*uw - wv*uv) / (vv*ww - wv*wv);

        if (beta >= 0 && gamma >= 0 && (beta + gamma <= 1) && t >= 0.0 && t <= tmax) {
            return true;
        } else {
            return false;
        }
    }
    return false;
}


int Triangle::intersectionPoint(Ray* ray, Intersection* intersection)
{
    vec3 dir = vec3(ray->direction[0],ray->direction[1],ray->direction[2]);
    vec3 pos = vec3(ray->position[0],ray->position[1],ray->position[2]);
    vec3 normal;
    double denom = glm::dot(dir,faceNormal);

    if (denom == 0.0) {
        return 0;
    } else {
        double t = (glm::dot(vertexArray[1],faceNormal) - glm::dot(pos,faceNormal)) / denom;

        intersection->t = t;
        vec3 p = pos + dir * t;
        intersection->position = vec4(p,1);

        // Barycentric coordinate calculations
        vec3 u = p - vertexArray[1];             
        vec3 v = vertexArray[0] - vertexArray[1];
        vec3 w = vertexArray[2] - vertexArray[1];

        double ww = glm::dot(w,w);
        double wv = glm::dot(w,v);
        double uv = glm::dot(u,v);
        double uw = glm::dot(u,w);
        double vv = glm::dot(v,v);
        double beta = (ww*uv - wv*uw) / (vv*ww - wv*wv);
        double gamma = (vv*uw - wv*uv) / (vv*ww - wv*wv);

        if (isFace) {
            normal = faceNormal;
        } else {
            normal = normalArray[1] + beta * (normalArray[0] - normalArray[1]) + 
                     gamma * (normalArray[2] - normalArray[1]);
        }
        intersection->normal = vec4(normal,0);

        if (beta >= 0 && gamma >= 0 && (beta + gamma <= 1)) {
            return 1;
        } else {
            return 0;
        }
    }
}
