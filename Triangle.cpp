/*
    Triangle.cpp
*/

#include <iostream>
#include "Triangle.h"

typedef glm::mat3 mat3;
typedef glm::mat4 mat4; 
typedef glm::vec3 vec3; 
typedef glm::vec4 vec4; 

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

bool Triangle::doesIntersect(Ray* ray, float tmax)
{
    vec3 dir = vec3(ray->direction[0],ray->direction[1],ray->direction[2]);
    vec3 pos = vec3(ray->position[0],ray->position[1],ray->position[2]);

    float denom = glm::dot(dir,faceNormal);

    if (denom == 0.0) {
        return false;
    } else {
        float t = (glm::dot(vertexArray[1],faceNormal) - glm::dot(pos,faceNormal)) / denom;
        vec3 p = pos + dir * t;

        vec3 u = p - vertexArray[1];             
        vec3 v = vertexArray[0] - vertexArray[1];
        vec3 w = vertexArray[2] - vertexArray[1];

        float ww = glm::dot(w,w);
        float wv = glm::dot(w,v);
        float uv = glm::dot(u,v);
        float uw = glm::dot(u,w);
        float vv = glm::dot(v,v);

        float beta = (ww*uv - wv*uw) / (vv*ww - wv*wv);
        float gamma = (vv*uw - wv*uv) / (vv*ww - wv*wv);

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
    float denom = glm::dot(dir,faceNormal);

    if (denom == 0.0) {
        return 0;
    } else {
        float t = (glm::dot(vertexArray[1],faceNormal) - glm::dot(pos,faceNormal)) / denom;

        intersection->t = t;
        vec3 p = pos + dir * t;
        intersection->position = vec4(p,1);

        // Barycentric coordinate calculations
        vec3 u = p - vertexArray[1];             
        vec3 v = vertexArray[0] - vertexArray[1];
        vec3 w = vertexArray[2] - vertexArray[1];

        float ww = glm::dot(w,w);
        float wv = glm::dot(w,v);
        float uv = glm::dot(u,v);
        float uw = glm::dot(u,w);
        float vv = glm::dot(v,v);
        float beta = (ww*uv - wv*uw) / (vv*ww - wv*wv);
        float gamma = (vv*uw - wv*uv) / (vv*ww - wv*wv);

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

void Triangle::createAABB(glm::mat3& transformation) 
{
    float min[3] = { 0.0, 0.0, 0.0 };
    float max[3] = { 0.0, 0.0, 0.0 };
    min[0] = max[0] = this->vertexArray[0].x;
    min[1] = max[1] = this->vertexArray[0].y;
    min[2] = max[2] = this->vertexArray[0].z;
    for (int i = 1; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            if (this->vertexArray[i][k] < min[k]) {
                min[k] = this->vertexArray[i][k];
            } else if (this->vertexArray[i][k] > max[k]) {
                max[k] = this->vertexArray[i][k];
            }
        }
    }
    this->aabb = new AABB();
    aabb->minimum = vec3(min[0],min[1],min[2])*transformation;
    aabb->maximum = vec3(max[0],max[1],max[2])*transformation;
}

AABB* Triangle::getAABB() {
    return aabb;
}
